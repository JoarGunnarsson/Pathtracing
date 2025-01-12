#include <iostream>
#include <fstream>
#include <chrono>
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
    Camera* camera;
    int number_of_objects;
    MaterialManager* material_manager;
};


struct PixelData{
    vec3 pixel_color = vec3(0,0,0);
    vec3 pixel_position = vec3(0,0,0);
    vec3 pixel_normal = vec3(0,0,0);
};


vec3 direct_lighting(Hit& hit, Object** objects, const int number_of_objects){
    int light_source_idx_array[number_of_objects];

    int number_of_light_sources = 0;
    
    for (int i = 0; i < number_of_objects; i++){
        if (objects[i] -> is_light_source()){
            light_source_idx_array[number_of_light_sources] = i;
            number_of_light_sources++;
        }
    }
    if (number_of_light_sources == 0){
        return colors::BLACK;
    }
    
    int random_index = random_int(0, number_of_light_sources);
    int light_index = light_source_idx_array[random_index];

    double inverse_PDF;
    vec3 random_point = objects[light_index] -> random_light_point(hit.intersection_point, inverse_PDF);

    vec3 vector_towards_light = random_point - hit.intersection_point;
    double distance_to_light = vector_towards_light.length();
    vector_towards_light = normalize_vector(vector_towards_light);
    hit.outgoing_vector = vector_towards_light;
    
    Ray light_ray;
    light_ray.starting_position = hit.intersection_point;
    light_ray.direction_vector =  vector_towards_light;

    Hit light_hit = find_closest_hit(light_ray, objects, number_of_objects);
    bool in_shadow = light_hit.intersected_object_index != light_index;
    bool same_distance = std::abs(distance_to_light - light_hit.distance) <= constants::EPSILON;
    bool hit_from_behind = dot_vectors(vector_towards_light, hit.normal_vector) < 0.0;
    bool inside_object = dot_vectors(hit.incoming_vector, hit.normal_vector) > 0.0;
    if ( in_shadow || !same_distance || hit_from_behind || inside_object){
        return colors::BLACK;
    }

    vec3 brdf_multiplier = objects[hit.intersected_object_index] -> eval(hit);
    vec3 light_emitance = objects[light_index] -> get_light_emittance(light_hit);

    double cosine = dot_vectors(hit.normal_vector, vector_towards_light);
    cosine = std::max(0.0, cosine);

    return brdf_multiplier * cosine * light_emitance * inverse_PDF * (double) number_of_light_sources;
}


PixelData raytrace(Ray ray, Object** objects, const int number_of_objects){
    PixelData data;
    vec3 color = vec3(0,0,0);
    vec3 brdf_accumulator = vec3(1,1,1);
    vec3 throughput = vec3(1,1,1);
    int force_recusion_limit = 3;
    double random_threshold = 1;
    for (int depth = 0; depth <= constants::max_recursion_depth; depth++){
        Hit ray_hit = find_closest_hit(ray, objects, number_of_objects);
        if (ray_hit.distance <= constants::EPSILON){
            break;
        }

        if (depth == 0){
            data.pixel_position = ray_hit.intersection_point;
            data.pixel_normal = ray_hit.normal_vector;
        }
        
        bool is_specular_ray = ray.type == REFLECTED || ray.type == TRANSMITTED;

        if (!constants::enable_next_event_estimation || depth == 0 || is_specular_ray){
            vec3 light_emitance = objects[ray_hit.intersected_object_index] -> get_light_emittance(ray_hit);
            color += light_emitance * brdf_accumulator * (dot_vectors(ray.direction_vector, ray_hit.normal_vector) < 0);
        }
        bool allow_recursion;
        double random_threshold;

        if (depth < force_recusion_limit){
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

        if (constants::enable_next_event_estimation){
            vec3 direct = direct_lighting(ray_hit, objects, number_of_objects);
            color += direct * brdf_accumulator / random_threshold;
        }

        BrdfData brdf_result = objects[ray_hit.intersected_object_index] -> sample(ray_hit, objects, number_of_objects);

        throughput *= brdf_result.brdf_multiplier / random_threshold;
        ray.starting_position = ray_hit.intersection_point;
        ray.direction_vector = brdf_result.outgoing_vector;
        ray.type = brdf_result.type;

        brdf_accumulator *= brdf_result.brdf_multiplier / random_threshold;
    }

    data.pixel_color = color;
    return data;
 }


PixelData compute_pixel_color(const int x, const int y, Scene scene, std::ofstream& file){
    PixelData data;
    vec3 pixel_color = vec3(0,0,0);
    Ray ray;
    ray.starting_position = scene.camera -> position;
    for (int i = 0; i < constants::samples_per_pixel; i++){
        double x_offset = random_normal() / 2.0;
        double y_offset = random_normal() / 2.0;
        double new_x = x + x_offset;
        double new_y = y + y_offset;
        ray.direction_vector = scene.camera -> get_starting_directions(new_x, new_y);
        PixelData sampled_data = raytrace(ray, scene.objects, scene.number_of_objects);
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
    gold_data.albedo_map = new ValueMap3D(colors::WHITE);
    gold_data.roughness_map = new ValueMap1D(0);
    gold_data.refractive_index = 0.277;
    gold_data.extinction_coefficient = 2.92;
    gold_data.is_dielectric = false;
    DiffuseMaterial* gold_material = new DiffuseMaterial(gold_data);
    manager -> add_material(gold_material);

    MaterialData light_material_data;
    light_material_data.albedo_map =  new ValueMap3D(colors::WHITE * 0.8);
    light_material_data.emission_color_map =  new ValueMap3D(colors::WARM_WHITE);
    light_material_data.light_intensity_map = new ValueMap1D(15.0);
    light_material_data.is_light_source = true;
    DiffuseMaterial* light_source_material = new DiffuseMaterial(light_material_data);
    manager -> add_material(light_source_material);

    MaterialData glass_data;
    glass_data.refractive_index = 1.5;
    TransparentMaterial* glass_material = new TransparentMaterial(glass_data);
    manager -> add_material(glass_material);

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

    
    Plane* this_floor = new Plane(vec3(0,0,0), vec3(1,0,0), vec3(0,0,-1), sandstone_material);
    Rectangle* front_wall = new Rectangle(vec3(0,1.55,-0.35), vec3(1,0,0), vec3(0,1,0), 2, 1.55*2, hieroglyph_material);
    Rectangle* left_wall = new Rectangle(vec3(-1,1.55,1.575), vec3(0,0,-1), vec3(0,1,0), 3.85, 1.55*2, hieroglyph_material);
    Rectangle* right_wall = new Rectangle(vec3(1,1.55,1.575), vec3(0,0,1), vec3(0,1,0), 3.85, 1.55*2, hieroglyph_material);
    Plane* roof = new Plane(vec3(0,2.2,0), vec3(1,0,0), vec3(0,0,1), white_diffuse_material);
    Rectangle* back_wall = new Rectangle(vec3(0,0.425,3.5), vec3(0,1,0), vec3(1,0,0), 2, 3.85/2.0, hieroglyph_material);

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

    Sphere* ball1 = new Sphere(vec3(-0.35,0,0), 0.35, green_diffuse_material);
    //Sphere* ball2 = new Sphere(vec3(0.45, 0, 0.6), 0.35, light_source_material);

    Rectangle* light_source = new Rectangle(vec3(0, 2.199, 1), vec3(0,0,-1), vec3(1,0,0), 1, 1, light_source_material);
    
    double desired_size = 0.5;
    vec3 desired_center = vec3(0, 0.8, 1);
    bool smooth_shade = true;
    ObjectUnion* loaded_model = load_object_model("./models/anubis.obj", gold_material, smooth_shade, desired_center, desired_size);

    int number_of_objects = 8;
    Object** objects = new Object*[number_of_objects]{this_floor, front_wall, left_wall, right_wall, roof, back_wall, light_source, loaded_model};

    vec3 camera_position = vec3(-1, 1, 2.2);
    vec3 viewing_direction = vec3(0.8, -0.3, -1);
    vec3 screen_y_vector = vec3(0, 1, 0);
    Camera* camera = new Camera(camera_position, viewing_direction, screen_y_vector);

    Scene scene;
    scene.objects = objects;
    scene.camera = camera;
    scene.number_of_objects = number_of_objects;
    scene.material_manager = manager;
    return scene;
}


void clear_scene(Scene& scene){
    for (int i = 0; i < scene.number_of_objects; i++){
        delete scene.objects[i];
    }

    delete[] scene.objects;
    delete scene.material_manager;
    delete scene.camera;
}


int main() {

    std::ofstream data_file;
    data_file.open("./temp/result_data.txt");
    
    data_file << "SIZE:" << constants::WIDTH << ' ' << constants::HEIGHT << "\n";

    Scene scene = create_scene();
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    vec3* pixel_buffer = new vec3[constants::WIDTH * constants::HEIGHT];
    vec3* position_buffer = new vec3[constants::WIDTH * constants::HEIGHT];
    vec3* normal_buffer = new vec3[constants::WIDTH * constants::HEIGHT];

    for (int y = 0; y < constants::HEIGHT; y++) {
        double progress = double(y) / (double) constants::HEIGHT;
        print_progress(progress);
        for (int x = 0; x < constants::WIDTH; x++) {
            PixelData data = compute_pixel_color(x, y, scene, data_file);
            int idx = ((constants::HEIGHT - 1 - y) * constants::WIDTH) + x;
            pixel_buffer[idx] = data.pixel_color;
            position_buffer[idx] = data.pixel_position;
            normal_buffer[idx] = data.pixel_normal;
        }
    }

    denoise(pixel_buffer, position_buffer, normal_buffer);

    for (int i=0; i < constants::WIDTH * constants::HEIGHT; i++){
        print_pixel_color(tone_map(pixel_buffer[i]), data_file);
    }

    data_file.close();
    std::clog << std::endl;
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    std::clog << "Time taken: " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]" << std::endl;

    clear_scene(scene);
    return 0;
}