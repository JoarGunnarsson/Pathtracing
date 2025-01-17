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


int sample_random_light(Object** objects, const int number_of_objects, int& number_of_light_sources){
    int light_source_idx_array[number_of_objects];

    number_of_light_sources = 0;
    
    for (int i = 0; i < number_of_objects; i++){
        if (objects[i] -> is_light_source()){
            light_source_idx_array[number_of_light_sources] = i;
            number_of_light_sources++;
        }
    }
    if (number_of_light_sources == 0){
        return -1;
    }
    
    int random_index = random_int(0, number_of_light_sources);
    int light_index = light_source_idx_array[random_index];
    return light_index;
}

vec3 direct_lighting(const Hit& hit, Object** objects, const int number_of_objects){
    int number_of_light_sources;
    int light_index = sample_random_light(objects, number_of_objects, number_of_light_sources);
    if (light_index == -1){
        return colors::BLACK;
    }

    double inverse_PDF;
    vec3 random_point = objects[light_index] -> random_light_point(hit.intersection_point, inverse_PDF);

    vec3 vector_towards_light = random_point - hit.intersection_point;
    double distance_to_light = vector_towards_light.length();
    vector_towards_light = normalize_vector(vector_towards_light);
    
    Ray light_ray;
    light_ray.starting_position = hit.intersection_point;
    light_ray.direction_vector =  vector_towards_light;
    light_ray.prepare();

    // This does not need to be a full find_closest hit, we can return early. Should allow parameter ignore_index.
    Hit light_hit = find_closest_hit(light_ray, objects, number_of_objects);
    bool in_shadow = light_hit.intersected_object_index != light_index;
    bool same_distance = std::abs(distance_to_light - light_hit.distance) <= constants::EPSILON;
    //bool hit_from_behind = dot_vectors(vector_towards_light, hit.normal_vector) < 0.0;
    bool inside_object = dot_vectors(hit.incident_vector, hit.normal_vector) > 0.0;

    if ( in_shadow || !same_distance || inside_object){
        return colors::BLACK;
    }

    vec3 brdf_multiplier = objects[hit.intersected_object_index] -> eval(hit);
    vec3 light_emitance = objects[light_index] -> get_light_emittance(light_hit);

    double cosine = dot_vectors(hit.normal_vector, vector_towards_light);
    cosine = std::max(0.0, cosine);

    return brdf_multiplier * cosine * light_emitance * inverse_PDF * (double) number_of_light_sources;
}


PixelData raytrace(Ray ray, Object** objects, const int number_of_objects, Medium* background_medium){
    MediumStack* medium_stack = new MediumStack();
    //medium_stack -> add_medium(background_medium, -1);
    PixelData data;
    vec3 color = vec3(0,0,0);
    vec3 throughput = vec3(1,1,1);
    double random_threshold = 1;
    bool allow_recursion = true;

    for (int depth = 0; depth <= constants::max_recursion_depth; depth++){
        ray.prepare();
        Hit ray_hit = find_closest_hit(ray, objects, number_of_objects);
        if (ray_hit.distance <= constants::EPSILON){
            break;
        }

        if (depth == 0){
            data.pixel_position = ray_hit.intersection_point;
            data.pixel_normal = ray_hit.normal_vector;
        }

        Object* hit_object = objects[ray_hit.intersected_object_index];
        /*
        double scatter_distance = medium -> sample_distance();
        // Add attenuation here. Take the minimum of scatter_distance and ray_hit.distance.
        Material* hit_material = objects[ray_hit.intersected_object_index] -> get_material(ray_hit.primitive_ID);
        double actual_distance;
        if (scatter_distance == -1){
            actual_distance = ray_hit.distance;
        }
        else{
            actual_distance = std::min(scatter_distance, ray_hit.distance);
        }
        vec3 log_attenuation = medium -> absorption_albedo * medium -> attenuation_coefficient * (-actual_distance);
        vec3 attenuation_color = exp_vector(log_attenuation);

        bool scatter = scatter_distance < ray_hit.distance && scatter_distance != -1;
        if (scatter){
            // Scatter
            ray.starting_position = ray.starting_position + ray.direction_vector * scatter_distance;
            ray.direction_vector = medium -> sample_direction(ray.direction_vector);
        }     
        */
        bool is_specular_ray = ray.type == REFLECTED || ray.type == TRANSMITTED;

        if (!constants::enable_next_event_estimation || depth == 0 || is_specular_ray){
            vec3 light_emitance = hit_object -> get_light_emittance(ray_hit);
            color += light_emitance * throughput * (dot_vectors(ray.direction_vector, ray_hit.normal_vector) < 0);
        }

        // Next event estimation does not work when inside a medium. Needs fixing.
        if (constants::enable_next_event_estimation){
            vec3 direct = direct_lighting(ray_hit, objects, number_of_objects);
            color += direct * throughput / random_threshold;
        }

        BrdfData brdf_result = hit_object -> sample(ray_hit);

        throughput *= brdf_result.brdf_multiplier;

        double incoming_dot_normal = dot_vectors(ray_hit.incident_vector, ray_hit.normal_vector);
        double outgoing_dot_normal = dot_vectors(brdf_result.outgoing_vector, ray_hit.normal_vector);
        
        bool penetrating_boundary = incoming_dot_normal * outgoing_dot_normal > 0;
        bool entering = incoming_dot_normal < 0;
        if (penetrating_boundary){
            // Something about this is not really working, tries to pop medium while medium is not in stack. We enter multple times too.
            // Seems to be an issue with concave objects, since the issue is not present for convex object unions (sphere etc).
            // Somehow, we can get total internal reflection from the inside even if we have entered and then left. 
            // Seems we think we leave but we actually dont.
            // I think this effect can be traced to sampling the outgoing vectors.
            // Indicates the issue is with convex shapes, which I suspected 6 hours ago.
            // Could also be non-water-tightness, would explain strange vertical line when close to object.
            if (entering){
                medium_stack -> add_medium(hit_object -> get_material(ray_hit.primitive_ID) -> medium, ray_hit.intersected_object_index);
            }
            else{
                medium_stack -> pop_medium(ray_hit.intersected_object_index);
            }
        }

        Ray new_ray;
        new_ray.starting_position = ray_hit.intersection_point;
        new_ray.direction_vector = brdf_result.outgoing_vector;
        new_ray.type = brdf_result.type;
        new_ray.prepare();

        Medium* medium = medium_stack -> get_medium();
        if (medium){
            vec3 Lv;
            vec3 transmittance; 
            vec3 weight;
            medium -> Integrate(objects, number_of_objects, new_ray, Lv, transmittance, weight, new_ray);
            color += Lv * weight * throughput;
            throughput *= transmittance;
        }

        // Russian Roulette
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
        ray = new_ray;
    }
    
    data.pixel_color = color;
    delete medium_stack;
    return data;
 }


PixelData compute_pixel_color(const int x, const int y, const Scene& scene){
    PixelData data;
    vec3 pixel_color = vec3(0,0,0);
    Ray ray;
    ray.starting_position = scene.camera -> position;
    for (int i = 0; i < constants::samples_per_pixel; i++){
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
    if (progress < 1.0) {
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

    MaterialData hieroglyph_data;
    hieroglyph_data.albedo_map = create_value_map_3D("./maps/hieroglyph_wall.map", 1, 1);
    DiffuseMaterial* hieroglyph_material = new DiffuseMaterial(hieroglyph_data);
    manager -> add_material(hieroglyph_material);

    MaterialData sandstone_data;
    sandstone_data.albedo_map = create_value_map_3D("./maps/sandstone_floor.map");
    DiffuseMaterial* sandstone_material = new DiffuseMaterial(sandstone_data);
    manager -> add_material(sandstone_material);

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
    light_material_data.light_intensity_map = new ValueMap1D(15.0 / 4.0);
    light_material_data.is_light_source = true;
    DiffuseMaterial* light_source_material = new DiffuseMaterial(light_material_data);
    manager -> add_material(light_source_material);

    MaterialData glass_data;
    glass_data.refractive_index = 1.5;
    glass_data.absorption_albedo = vec3(1,1,1) - colors::BLUE;
    glass_data.attenuation_coefficient = 1;
    TransparentMaterial* glass_material = new TransparentMaterial(glass_data);
    manager -> add_material(glass_material);

    MaterialData scattering_glass_data;
    scattering_glass_data.refractive_index = 1;
    scattering_glass_data.scattering_coefficient = 10;
    scattering_glass_data.attenuation_coefficient = 1;
    scattering_glass_data.absorption_albedo = vec3(1,1,1) - colors::BLUE;
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
    */

    MaterialData model_data;
    model_data.albedo_map = create_value_map_3D("./maps/bunny.map", 1, -1);
    DiffuseMaterial* model_material = new DiffuseMaterial(model_data);
    manager -> add_material(model_material);

    
    Plane* this_floor = new Plane(vec3(0,0,0), vec3(1,0,0), vec3(0,0,-1), white_diffuse_material);
    Rectangle* front_wall = new Rectangle(vec3(0,1.55,-0.35), vec3(1,0,0), vec3(0,1,0), 2, 1.55*2, white_diffuse_material);
    Rectangle* left_wall = new Rectangle(vec3(-1,1.55,1.575), vec3(0,0,-1), vec3(0,1,0), 3.85, 1.55*2, red_diffuse_material);
    Rectangle* right_wall = new Rectangle(vec3(1,1.55,1.575), vec3(0,0,1), vec3(0,1,0), 3.85, 1.55*2, green_diffuse_material);
    Plane* roof = new Plane(vec3(0,2.2,0), vec3(1,0,0), vec3(0,0,1), white_diffuse_material);
    Rectangle* back_wall = new Rectangle(vec3(0,0.425,3.5), vec3(0,1,0), vec3(1,0,0), 2, 3.85/2.0, white_diffuse_material);

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

    Sphere* ball1 = new Sphere(vec3(-0.35, 0.5, 0), 0.35, green_diffuse_material);
    Sphere* ball2 = new Sphere(vec3(0.45, 0.5, 0.6), 0.35, glass_material);

    Rectangle* light_source = new Rectangle(vec3(0, 2.199, 1), vec3(0,0,-1), vec3(1,0,0), 2, 2, light_source_material);
    
    double desired_size = 0.5;
    vec3 desired_center = vec3(0, 0.8, 1);
    bool smooth_shade = false;
    ObjectUnion* loaded_model = load_object_model("./models/dragon.obj", glass_material, smooth_shade, desired_center, desired_size);

    int number_of_objects = 8;
    Object** objects = new Object*[number_of_objects]{this_floor, front_wall, left_wall, right_wall, roof, back_wall, light_source, loaded_model};

    Medium* background_medium = new Medium(0, 0, colors::WHITE);

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

    raw_data_file << "SIZE:" << constants::WIDTH << ' ' << constants::HEIGHT << "\n";

    Scene scene = create_scene();
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