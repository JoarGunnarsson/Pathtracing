#include <iostream>
#include <fstream>
#include <chrono>
#include <stdexcept>
#include "vec3.h"
#include "colors.h"
#include "denoise.h"
#include "objects.h"
#include "camera.h"
#include "utils.h"
#include "constants.h"
#include "objectunion.h"


void print_pixel_color(const vec3& rgb, std::ofstream& file){
    file << rgb[0] << ' ' << rgb[1] << ' ' << rgb[2] << '\n';
}


struct Scene{
    Object** objects;
    int number_of_objects;
    Camera* camera;
    MaterialManager* material_manager;
    Medium* medium;
};


struct PixelData{
    vec3 pixel_color = vec3(0,0,0);
    vec3 pixel_position = vec3(0,0,0);
    vec3 pixel_normal = vec3(0,0,0);
};


PixelData raytrace(Ray ray, Object** objects, const int number_of_objects, Medium* background_medium){
    MediumStack medium_stack = MediumStack();
    medium_stack.add_medium(background_medium, -1);
    PixelData data;
    vec3 color = vec3(0,0,0);
    vec3 throughput = vec3(1,1,1);
    double random_threshold = 1;
    bool allow_recursion = true;
    bool has_hit_surface = false;

    vec3 saved_point;
    double brdf_pdf;

    for (int depth = 0; depth <= constants::max_recursion_depth; depth++){
        Medium* medium = medium_stack.get_medium();
        double scatter_distance = medium -> sample_distance();

        ray.t_max = scatter_distance;
        Hit ray_hit;
        if (!find_closest_hit(ray_hit, ray, objects, number_of_objects)){
            if (scatter_distance == constants::max_ray_distance){
                break;
            }
            ray_hit.distance = constants::max_ray_distance;
        }
        
        bool scatter = scatter_distance < ray_hit.distance;
        scatter_distance = std::min(scatter_distance, ray_hit.distance);
        throughput *= medium -> sample(objects, number_of_objects, scatter_distance, scatter);
        
        if (scatter){
            vec3 scatter_point = ray.starting_position + ray.direction_vector * scatter_distance;
            vec3 scattered_direction = medium -> sample_direction(ray.direction_vector);
            if (constants::enable_next_event_estimation){
                ray_hit.intersection_point = scatter_point;
                // TODO: This is reached for scattering mediums with no scattering, is that really correct? Look over distance sampling...

                // TODO: Compute direct light using MIS for scattering. Should be identical to diffuse objects except the 
                // phase function is used instead of the brdf.
                color += compute_direct_light_scattering(ray_hit, objects, number_of_objects, medium_stack) * throughput;
                // TODO: In sample_direct, do some MIS. The other part of MIS is done if it hits an object?
                // If we do some direct light, we turn the ray diffuse. This is inherited by the next ray for n==1 transparent materials.
                // Save the phase function for the scattered ray, this will be used by MIS later.
                ray.type = DIFFUSE;
                brdf_pdf = medium -> phase_function(ray.direction_vector, scattered_direction);
            }
            ray.starting_position = scatter_point;
            ray.direction_vector = scattered_direction;
        }     
        else{
            if (!has_hit_surface){
                data.pixel_position = ray_hit.intersection_point;
                data.pixel_normal = ray_hit.normal_vector;
                has_hit_surface = true;
            }

            bool is_specular_ray = ray.type == REFLECTED || ray.type == TRANSMITTED;
            Object* hit_object = objects[ray_hit.intersected_object_index];

            //TODO: The code below should be used for scattered ray direct light too. Instead of brdf_pdf, use phase_function pdf?

            // If a light source is hit, compute the light_pdf based on the saved_point (previous hitpoint) and use MIS to add the light.
            // Could move this to a separate function, make it clearer what it is doing.
            if (hit_object -> is_light_source()){
                double weight;
                if (!constants::enable_next_event_estimation || depth == 0 || is_specular_ray){
                    weight = 1;
                }
                else{
                    double light_pdf = objects[ray_hit.intersected_object_index] -> light_pdf(ray_hit.intersection_point, saved_point, ray_hit.primitive_ID);
                    weight = mis_weight(1, brdf_pdf, 1, light_pdf);
                }
                vec3 light_emittance = hit_object -> get_light_emittance(ray_hit);

                color += weight * light_emittance * throughput * (dot_vectors(ray.direction_vector, ray_hit.normal_vector) < 0);
            }
            
            if (constants::enable_next_event_estimation){
                color += compute_direct_light(ray_hit, objects, number_of_objects, medium_stack) * throughput;//hit_object -> sample_direct(ray_hit, objects, number_of_objects, medium_stack) * throughput;
            }
            
            // Sample the material for a new direction.
            BrdfData brdf_result = hit_object -> sample(ray_hit);
            // Check if the hit material allows direct light computations to pass through the object (transparent objects with n==1).
            // If so, the outgoing ray inherits the incoming ray type. Otherwise, the brdf_pdf used in MIS is updated.
            // TODO: Rename allow_direct_light!
            // TODO: Rename is_virtual_surface variable...
            bool is_virtual_surface = hit_object -> get_material(ray_hit.primitive_ID) -> allow_direct_light();
            if (is_virtual_surface){
                brdf_result.type = ray.type;
            }
            else{ 
                brdf_pdf = brdf_result.pdf;
                // TODO: This was together with the other ray_hit stuff at the end, moved here for scattered ray MIS.
                saved_point = ray_hit.intersection_point;
            }
            throughput *= brdf_result.brdf_over_pdf;

            double incoming_dot_normal = dot_vectors(ray_hit.incident_vector, ray_hit.normal_vector);
            double outgoing_dot_normal = dot_vectors(brdf_result.outgoing_vector, ray_hit.normal_vector);
            
            bool penetrating_boundary = incoming_dot_normal * outgoing_dot_normal > 0;
            bool entering = incoming_dot_normal < 0;

            // TODO: Do the below part before sampling, so we can get the correct medium for refractive index etc?
            // TODO: Can save current_medium and next_medium, and pass that into sample and compute_direct_light.
            Medium* new_medium = hit_object -> get_material(ray_hit.primitive_ID) -> medium;
            if (penetrating_boundary && new_medium){
                // Something about this is not really working, tries to pop medium while medium is not in stack. We enter multple times too.
                // Seems to be an issue with concave objects, since the issue is not present for convex object unions (sphere etc).
                // Probably due to numeric errors. Currently relatively rare, so can be ignored, but not very good.
                if (entering){
                    medium_stack.add_medium(new_medium, ray_hit.intersected_object_index);
                }
                else{
                    medium_stack.pop_medium(ray_hit.intersected_object_index);
                }
            }
            ray.starting_position = ray_hit.intersection_point;
            ray.direction_vector = brdf_result.outgoing_vector;
            ray.type = brdf_result.type;
            // TODO: For scattered rays, this needs to be handled with care. Perhaps don't update for transparent/specular rays?...
        }

        if (depth < constants::force_tracing_limit){
            random_threshold = 1;
            allow_recursion = true;
        }
        else{
            random_threshold = std::min(throughput.max(), 0.9);
            double random_value = random_uniform(0, 1);
            allow_recursion = random_value < random_threshold;
        }
        
        if (!allow_recursion){
            break;
        }

        throughput /= random_threshold;
    }
    
    data.pixel_color = color;
    return data;
 }


PixelData compute_pixel_color(const int x, const int y, const Scene& scene){
    PixelData data;
    vec3 pixel_color = vec3(0,0,0);
    for (int i = 0; i < constants::samples_per_pixel; i++){
        Ray ray;
        ray.starting_position = scene.camera -> position;
        ray.type = TRANSMITTED;
        double new_x = x;
        double new_y = y;

        if (constants::enable_anti_aliasing){
            new_x += random_normal() / 3.0;
            new_y += random_normal() / 3.0;
        }
        
        ray.direction_vector = scene.camera -> get_starting_directions(new_x, new_y);
        PixelData sampled_data = raytrace(ray, scene.objects, scene.number_of_objects, scene.medium);
        if (i == 0){
            data.pixel_position = sampled_data.pixel_position;
            data.pixel_normal = sampled_data.pixel_normal;
        }
        pixel_color += sampled_data.pixel_color;    
    }

    vec3 average_pixel_color = pixel_color / (double) constants::samples_per_pixel;

    data.pixel_color = average_pixel_color;
    return data;
}


void print_progress(double progress){
    if (progress <= 1.0) {
        int bar_width = 60;

        std::clog << "[";
        int pos = bar_width * progress;
        for (int i = 0; i < bar_width; ++i) {
            if (i < pos) std::clog << "=";
            else if (i == pos) std::clog << ">";
            else std::clog << " ";
        }
        std::clog << "] " << int(progress * 100.0) << " %\r";
    }
}


Scene create_scene(){
    /*
    ValueMap3D* world_map = create_value_map_3D("./maps/world.map");
    ValueMap3D* sakura_map = create_value_map_3D("./maps/sakura.map");
    ValueMap3D* temple_map = create_value_map_3D("./maps/temple.map");
    ValueMap3D* cobble_map = create_value_map_3D("./maps/cobblestone.map");
    ValueMap1D* world_roughnessMap = create_value_map_1D("./maps/world_roughness.map");
    */

    MaterialManager* manager = new MaterialManager();
    MaterialData white_data;
    white_data.albedo_map = new ValueMap3D(colors::WHITE * 0.8);
    DiffuseMaterial* white_diffuse_material = new DiffuseMaterial(white_data);
    manager -> add_material(white_diffuse_material);
    
    MaterialData white_reflective_data;
    white_reflective_data.albedo_map = new ValueMap3D(colors::WHITE * 0.8);
    ReflectiveMaterial* white_reflective_material = new ReflectiveMaterial(white_reflective_data);   
    manager -> add_material(white_reflective_material);

    MaterialData red_material_data;
    red_material_data.albedo_map = new ValueMap3D(colors::RED);
    DiffuseMaterial* red_diffuse_material = new DiffuseMaterial(red_material_data);
    manager -> add_material(red_diffuse_material);

    MaterialData green_material_data;
    green_material_data.albedo_map = new ValueMap3D(colors::GREEN);
    DiffuseMaterial* green_diffuse_material = new DiffuseMaterial(green_material_data);
    manager -> add_material(green_diffuse_material);


    /*
    MaterialData hieroglyph_data;
    hieroglyph_data.albedo_map = create_value_map_3D("./maps/hieroglyph_wall.map", 1, 1);
    DiffuseMaterial* hieroglyph_material = new DiffuseMaterial(hieroglyph_data);
    manager -> add_material(hieroglyph_material);

    MaterialData sandstone_data;
    sandstone_data.albedo_map = create_value_map_3D("./maps/sandstone_floor.map");
    DiffuseMaterial* sandstone_material = new DiffuseMaterial(sandstone_data);
    manager -> add_material(sandstone_material);
    */

    MaterialData gold_data;
    gold_data.albedo_map = new ValueMap3D(colors::GOLD);
    gold_data.roughness_map = new ValueMap1D(0.1);
    gold_data.refractive_index = 0.277;
    gold_data.extinction_coefficient = 2.92;
    gold_data.is_dielectric = false;
    MicrofacetMaterial* gold_material = new MicrofacetMaterial(gold_data);
    manager -> add_material(gold_material);

    MaterialData light_material_data;
    light_material_data.albedo_map =  new ValueMap3D(colors::WHITE * 0.8);
    light_material_data.emission_color_map =  new ValueMap3D(colors::WARM_WHITE);
    light_material_data.light_intensity_map = new ValueMap1D(15.0);
    light_material_data.is_light_source = true;
    DiffuseMaterial* light_source_material = new DiffuseMaterial(light_material_data);
    manager -> add_material(light_source_material);

    MaterialData glass_data;
    glass_data.refractive_index = 1;
    BeersLawMedium* glass_medium = new BeersLawMedium(vec3(0), (vec3(1,1,1) - colors::BLUE) * 1.0);
    glass_data.medium = glass_medium;
    TransparentMaterial* glass_material = new TransparentMaterial(glass_data);
    manager -> add_material(glass_material);

    MaterialData scattering_glass_data;
    scattering_glass_data.refractive_index = 1;
    ScatteringMediumHomogenous* scattering_glass_medium = new ScatteringMediumHomogenous(vec3(4), (vec3(1,1,1)) * 1);
    scattering_glass_data.medium = scattering_glass_medium;
    TransparentMaterial* scattering_glass_material = new TransparentMaterial(scattering_glass_data);
    manager -> add_material(scattering_glass_material);

    MaterialData mirror_data;
    ReflectiveMaterial* mirror_material = new ReflectiveMaterial(mirror_data);

    /*
    MaterialData glass_data;
    glass_data.refractive_index = 1.5;
    TransparentMaterial* pane1_material = new TransparentMaterial(glass_data);

    MaterialData frosty_glass_data;
    frosty_glassData.albedo_map = pure_white_map;
    frosty_glassData.refractive_index = 1.5;
    frosty_glassData.roughness_map = new ValueMap1D(0.1);
    frosty_glassData.percentage_diffuse_map = new ValueMap1D(0.0);
    MicrofacetMaterial* _material = new MicrofacetMaterial(frosty_glassData);

    MaterialData sakura_data;
    sakura_data.albedo_map = sakura_map;
    DiffuseMaterial* sakura_material = new DiffuseMaterial(sakura_data);

    MaterialData temple_data;
    temple_data.albedo_map = temple_map;
    DiffuseMaterial* temple_material = new DiffuseMaterial(temple_data);

    MaterialData cobble_data;
    cobble_data.albedo_map = cobble_map;
    DiffuseMaterial* cobble_material = new DiffuseMaterial(cobble_data);
    

    MaterialData model_data;
    model_data.albedo_map = create_value_map_3D("./maps/bunny.map", 1, -1);
    DiffuseMaterial* model_material = new DiffuseMaterial(model_data);
    manager -> add_material(model_material);
    */

    //TODO: Fix back wall    
    Plane* this_floor = new Plane(vec3(0,0,0), vec3(1,0,0), vec3(0,0,-1), white_diffuse_material);
    Rectangle* front_wall = new Rectangle(vec3(0,1.55,-0.35), vec3(1,0,0), vec3(0,1,0), 2, 1.55*2, white_diffuse_material);
    Rectangle* left_wall = new Rectangle(vec3(-1,1.55,1.575), vec3(0,0,-1), vec3(0,1,0), 3.85, 1.55*2, red_diffuse_material);
    Rectangle* right_wall = new Rectangle(vec3(1,1.55,1.575), vec3(0,0,1), vec3(0,1,0), 3.85, 1.55*2, green_diffuse_material);
    Plane* roof = new Plane(vec3(0,2.2,0), vec3(1,0,0), vec3(0,0,1), white_diffuse_material);
    Rectangle* back_wall = new Rectangle(vec3(0,1.55,3.5), vec3(0,1,0), vec3(1,0,0), 3.85, 1.55*2, white_diffuse_material);

    /*
    Rectangle* front_pane1 = new Rectangle(vec3(-0.25,0.5,1.2), vec3(1,0,0), vec3(0,1,0), 0.5, 0.5, pane1_material);
    Rectangle* back_pane1 = new Rectangle(vec3(-0.25,0.5,1.15), vec3(-1,0,0), vec3(0,1,0), 0.5, 0.5, pane1_material);
    Object** pane1_objects = new Object*[2]{front_pane1, back_pane1};
    ObjectUnion* pane1 = new ObjectUnion(pane1_objects, 2);

    Rectangle* front_pane2 = new Rectangle(vec3(0.25,0.5,1.2), vec3(1,0,0), vec3(0,1,0), 0.5, 0.5, pane2_material);
    Rectangle* back_pane2 = new Rectangle(vec3(0.25,0.5,1.15), vec3(-1,0,0), vec3(0,1,0), 0.5, 0.5, pane2_material);
    Object** pane2_objects = new Object*[2]{front_pane2, back_pane2};
    ObjectUnion* pane2 = new ObjectUnion(pane2_objects, 2);
    */

    Sphere* ball1 = new Sphere(vec3(0, 0.8, 1), 0.35, glass_material);
    Sphere* ball2 = new Sphere(vec3(0.45, 0.5, 0.6), 0.35, scattering_glass_material);

    Rectangle* light_source = new Rectangle(vec3(0, 2.199, 1), vec3(0,0,-1), vec3(1,0,0), 1, 1, light_source_material);
    
    double desired_size = 0.5;
    vec3 desired_center = vec3(0, 0.8, 1);
    bool smooth_shade = false;
    //ObjectUnion* loaded_model = load_object_model("./models/bunny.obj", scattering_glass_material, smooth_shade, desired_center, desired_size);

    int number_of_objects = 8;
    Object** objects = new Object*[number_of_objects]{this_floor, front_wall, left_wall, right_wall, roof, back_wall, light_source, ball2};

    ScatteringMediumHomogenous* background_medium = new ScatteringMediumHomogenous(vec3(0.), (colors::WHITE) * 0.0);

    vec3 camera_position = vec3(-1, 1, 2.2);
    vec3 viewing_direction = vec3(0.8, -0.3, -1);
    vec3 screen_y_vector = vec3(0, 1, 0);
    Camera* camera = new Camera(camera_position, viewing_direction, screen_y_vector);

    Scene scene;
    scene.objects = objects;
    scene.camera = camera;
    scene.number_of_objects = number_of_objects;
    scene.material_manager = manager;
    scene.medium = background_medium;
    return scene;
}


void clear_scene(Scene& scene){
    for (int i = 0; i < scene.number_of_objects; i++){
        delete scene.objects[i];
    }

    delete[] scene.objects;
    delete scene.material_manager;
    delete scene.camera;
    delete scene.medium;
}


int main() {
    std::ofstream raw_data_file;
    raw_data_file.open(constants::raw_output_file_name);

    if (constants::enable_denoising){
            raw_data_file << "1\n";
    }
    else{
            raw_data_file << "0\n";
    }

    raw_data_file << "SIZE:" << constants::WIDTH << ' ' << constants::HEIGHT << "\n";
    std::chrono::steady_clock::time_point begin_build = std::chrono::steady_clock::now();

    Scene scene = create_scene();

    std::chrono::steady_clock::time_point end_build = std::chrono::steady_clock::now();
    std::clog << "Time taken to build scene: " << std::chrono::duration_cast<std::chrono::seconds>(end_build - begin_build).count() << "[s]" << std::endl;
    
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    vec3* pixel_buffer = new vec3[constants::WIDTH * constants::HEIGHT];
    vec3* position_buffer = new vec3[constants::WIDTH * constants::HEIGHT];
    vec3* normal_buffer = new vec3[constants::WIDTH * constants::HEIGHT];

    for (int y = constants::HEIGHT - 1; y >= 0; y--) {
        double progress = double(constants::HEIGHT - 1 - y) / (double) constants::HEIGHT;
        print_progress(progress);
        for (int x = 0; x < constants::WIDTH; x++) {
            PixelData data = compute_pixel_color(x, y, scene);
            int idx = ((constants::HEIGHT - 1 - y) * constants::WIDTH) + x;
            print_pixel_color(tone_map(data.pixel_color), raw_data_file);
            pixel_buffer[idx] = data.pixel_color;
            position_buffer[idx] = data.pixel_position;
            normal_buffer[idx] = data.pixel_normal;
        }
    }
    print_progress(1);

    raw_data_file.close();
 
    if (constants::enable_denoising){
        std::ofstream denoised_data_file;
        denoised_data_file.open(constants::denoised_output_file_name);
        
        denoised_data_file << "SIZE:" << constants::WIDTH << ' ' << constants::HEIGHT << "\n";

        denoise(pixel_buffer, position_buffer, normal_buffer);

        for (int i=0; i < constants::WIDTH * constants::HEIGHT; i++){
            print_pixel_color(tone_map(pixel_buffer[i]), denoised_data_file);
        }

        denoised_data_file.close();
    }

    std::clog << std::endl;
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    std::clog << "Time taken: " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]" << std::endl;

    clear_scene(scene);
    return 0;
}