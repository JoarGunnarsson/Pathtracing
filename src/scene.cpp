#include <nlohmann/json.hpp>
#include <fstream>

#include "scene.h"
#include "objects.h"
#include "objectunion.h"
#include "materials.h"
#include "camera.h"
#include "medium.h"

using json = nlohmann::json;

void load_settings() {
    std::ifstream settings_file("./scenes/settings.json");

    json data = json::parse(settings_file);

    if (data.contains("WIDTH")) {
        constants::WIDTH = data["WIDTH"];
    }

    if (data.contains("HEIGHT")) {
        constants::HEIGHT = data["HEIGHT"];
    }

    if (data.contains("samples_per_pixel")) {
        constants::samples_per_pixel = data["samples_per_pixel"];
    }

    if (data.contains("max_recursion_depth")) {
        constants::max_recursion_depth = data["max_recursion_depth"];
    }

    if (data.contains("force_tracing_limit")) {
        constants::force_tracing_limit = data["force_tracing_limit"];
    }

    if (data.contains("air_refractive_index")) {
        constants::air_refractive_index = data["air_refractive_index"];
    }

    if (data.contains("enable_next_event_estimation")) {
        constants::enable_next_event_estimation = data["enable_next_event_estimation"];
    }

    if (data.contains("enable_anti_aliasing")) {
        constants::enable_anti_aliasing = data["enable_anti_aliasing"];
    }

    if (data.contains("enable_denoising")) {
        constants::enable_denoising = data["enable_denoising"];
    }

    if (data.contains("denoising_iterations")) {
        constants::denoising_iterations = data["denoising_iterations"];
    }

    if (data.contains("sigma_rt")) {
        constants::denoising_iterations = data["denoising_iterations"];
    }

    if (data.contains("sigma_x")) {
        constants::denoising_iterations = data["denoising_iterations"];
    }

    if (data.contains("sigma_n")) {
        constants::denoising_iterations = data["denoising_iterations"];
    }
}

Scene create_scene_old() {
    MaterialManager* manager = new MaterialManager();
    MaterialData white_data;
    white_data.albedo_map = new ValueMap3D(colors::WHITE * 0.7);
    DiffuseMaterial* white_diffuse_material = new DiffuseMaterial(white_data);
    manager->add_material(white_diffuse_material);

    MaterialData white_reflective_data;
    white_reflective_data.albedo_map = new ValueMap3D(colors::WHITE * 0.8);
    ReflectiveMaterial* white_reflective_material = new ReflectiveMaterial(white_reflective_data);
    manager->add_material(white_reflective_material);

    MaterialData red_material_data;
    red_material_data.albedo_map = new ValueMap3D(colors::RED);
    DiffuseMaterial* red_diffuse_material = new DiffuseMaterial(red_material_data);
    manager->add_material(red_diffuse_material);

    MaterialData green_material_data;
    green_material_data.albedo_map = new ValueMap3D(colors::GREEN);
    DiffuseMaterial* green_diffuse_material = new DiffuseMaterial(green_material_data);
    manager->add_material(green_diffuse_material);

    MaterialData gold_data;
    gold_data.albedo_map = new ValueMap3D(colors::GOLD);
    gold_data.roughness_map = new ValueMap1D(0.3);
    gold_data.refractive_index = 0.277;
    gold_data.extinction_coefficient = 2.92;
    gold_data.is_dielectric = false;
    MetallicMicrofacet* gold_material = new MetallicMicrofacet(gold_data);
    manager->add_material(gold_material);

    MaterialData light_material_data;
    light_material_data.albedo_map = new ValueMap3D(colors::WHITE * 0.8);
    light_material_data.emission_color_map = new ValueMap3D(colors::WARM_WHITE);
    light_material_data.light_intensity_map = new ValueMap1D(200.0);
    light_material_data.is_light_source = true;
    DiffuseMaterial* light_source_material = new DiffuseMaterial(light_material_data);
    manager->add_material(light_source_material);

    MaterialData glass_data;
    glass_data.refractive_index = 1.33;
    BeersLawMedium* glass_medium = new BeersLawMedium(vec3(0), (vec3(1, 1, 1) - colors::BLUE) * 0.0, vec3(0));
    glass_data.medium = glass_medium;
    TransparentMaterial* glass_material = new TransparentMaterial(glass_data);
    manager->add_material(glass_material);

    MaterialData scattering_glass_data;
    scattering_glass_data.refractive_index = 1.;
    ScatteringMediumHomogenous* scattering_glass_medium =
        new ScatteringMediumHomogenous(colors::WHITE, vec3(2.7, 1., 1.1) * 0.0, vec3(0));
    scattering_glass_data.medium = scattering_glass_medium;
    TransparentMaterial* scattering_glass_material = new TransparentMaterial(scattering_glass_data);
    manager->add_material(scattering_glass_material);

    MaterialData mirror_data;
    ReflectiveMaterial* mirror_material = new ReflectiveMaterial(mirror_data);
    manager->add_material(mirror_material);

    Plane* this_floor = new Plane(vec3(0, -0.3, 0), vec3(1, 0, 0), vec3(0, 0, -1), white_diffuse_material);
    Rectangle* front_wall =
        new Rectangle(vec3(0, 1.55, -0.35), vec3(1, 0, 0), vec3(0, 1, 0), 2, 1.55 * 2 * 1.5, white_diffuse_material);
    Rectangle* left_wall =
        new Rectangle(vec3(-1, 1.55, 1.575), vec3(0, 0, -1), vec3(0, 1, 0), 3.85, 1.55 * 2 * 1.5, red_diffuse_material);
    Rectangle* right_wall =
        new Rectangle(vec3(1, 1.55, 1.575), vec3(0, 0, 1), vec3(0, 1, 0), 3.85, 1.55 * 2 * 1.5, green_diffuse_material);
    Plane* roof = new Plane(vec3(0, 2.2, 0), vec3(1, 0, 0), vec3(0, 0, 1), white_diffuse_material);
    Rectangle* back_wall =
        new Rectangle(vec3(0, 1.55, 3.5), vec3(0, 1, 0), vec3(1, 0, 0), 3.85, 1.55 * 2 * 1.5, white_diffuse_material);

    Sphere* ball1 = new Sphere(vec3(0, 0.5, 1.8), 0.3, glass_material);
    Plane* ball2 = new Plane(vec3(0, 0.19, 0), vec3(1, 0, 0), vec3(0, 0, -1), scattering_glass_material);

    Rectangle* light_source =
        new Rectangle(vec3(0, 2.199, 1.8), vec3(1, 0, 0), vec3(0, 0, 1), 0.4, 0.4, light_source_material);

    double desired_size = 0.6;
    vec3 desired_center = vec3(-0.3, 0.1, 1.3);
    bool smooth_shade = false;
    bool transform_object = true;
    // TODO: Actually, use struct called object_transform, can set it to nullptr if no transformation should be made.
    //ObjectUnion* loaded_model = load_object_model("./models/water_cube.obj", scattering_glass_material, smooth_shade, transform_object, desired_center, desired_size);

    int number_of_objects = 7;
    Object** objects =
        new Object* [number_of_objects] { this_floor, front_wall, left_wall, right_wall, roof, ball1, light_source };

    BeersLawMedium* background_medium = new BeersLawMedium(vec3(0.4), (colors::WHITE) *0.0, vec3(0));

    vec3 camera_position = vec3(0, 0.8, 4.2);
    vec3 viewing_direction = vec3(0, 0, -1);
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

Scene create_scene() {
    MaterialManager* manager = new MaterialManager();
    MaterialData white_data;
    white_data.albedo_map = new ValueMap3D(colors::WHITE * 0.7);
    DiffuseMaterial* white_diffuse_material = new DiffuseMaterial(white_data);
    manager->add_material(white_diffuse_material);

    MaterialData red_material_data;
    red_material_data.albedo_map = new ValueMap3D(colors::RED);
    DiffuseMaterial* red_diffuse_material = new DiffuseMaterial(red_material_data);
    manager->add_material(red_diffuse_material);

    MaterialData green_material_data;
    green_material_data.albedo_map = new ValueMap3D(colors::GREEN);
    ReflectiveMaterial* green_diffuse_material = new ReflectiveMaterial(green_material_data);
    manager->add_material(green_diffuse_material);

    MaterialData mirror_data;
    ReflectiveMaterial* mirror_material = new ReflectiveMaterial(mirror_data);
    manager->add_material(mirror_material);

    MaterialData light_material_data;
    light_material_data.albedo_map = new ValueMap3D(colors::WHITE * 0.8);
    light_material_data.emission_color_map = new ValueMap3D(colors::WARM_WHITE);
    light_material_data.light_intensity_map = new ValueMap1D(150.0);
    light_material_data.is_light_source = true;
    DiffuseMaterial* light_source_material = new DiffuseMaterial(light_material_data);
    manager->add_material(light_source_material);

    MaterialData glass_data;
    glass_data.refractive_index = 1.3;
    BeersLawMedium* glass_medium = new BeersLawMedium(vec3(0), (vec3(1, 1, 1) - colors::BLUE) * 0.0, vec3(0));
    glass_data.medium = glass_medium;
    TransparentMaterial* glass_material = new TransparentMaterial(glass_data);
    manager->add_material(glass_material);

    MaterialData whiskey_data;
    whiskey_data.refractive_index = 1.356;
    BeersLawMedium* whiskey_medium = new BeersLawMedium(vec3(0), (vec3(1.0) - vec3(0.1, 0.3, 0.44)) * 10, vec3(0));
    whiskey_data.medium = whiskey_medium;
    TransparentMaterial* whiskey_material = new TransparentMaterial(whiskey_data);
    manager->add_material(whiskey_material);

    MaterialData floor_data;
    floor_data.albedo_map = create_value_map_3D("./maps/wooden_floor.map", 1, -1);
    floor_data.roughness_map = new ValueMap1D(0.1);
    floor_data.refractive_index = 1.5;
    DiffuseMaterial* floor_material = new DiffuseMaterial(floor_data);
    manager->add_material(floor_material);

    MaterialData chair_data;
    chair_data.albedo_map = create_value_map_3D("./maps/wooden_chair.map", 1, -1);
    chair_data.roughness_map = new ValueMap1D(0.1);
    chair_data.refractive_index = 1.1;
    GlossyMaterial* chair_material = new GlossyMaterial(chair_data);
    manager->add_material(chair_material);

    MaterialData table_data;
    table_data.albedo_map = create_value_map_3D("./maps/wooden_table.map", 1, -1);
    table_data.roughness_map = new ValueMap1D(0.1);
    table_data.refractive_index = 1.5;
    GlossyMaterial* table_material = new GlossyMaterial(table_data);
    manager->add_material(table_material);

    MaterialData painting_ramen_data;
    painting_ramen_data.albedo_map = create_value_map_3D("./maps/ramen.map", 1, -1);
    DiffuseMaterial* painting_ramen_material = new DiffuseMaterial(painting_ramen_data);
    manager->add_material(painting_ramen_material);

    MaterialData sacco_data;
    sacco_data.albedo_map = new ValueMap3D(vec3(0.04));
    sacco_data.roughness_map = new ValueMap1D(0.25);
    sacco_data.refractive_index = 1.3;
    GlossyMaterial* sacco_material = new GlossyMaterial(sacco_data);
    manager->add_material(sacco_material);

    MaterialData brown_sacco_data;
    brown_sacco_data.albedo_map = create_value_map_3D("./maps/brown_sacco.map", 1, -1);
    brown_sacco_data.roughness_map = new ValueMap1D(0.25);
    brown_sacco_data.refractive_index = 1.3;
    GlossyMaterial* brown_sacco_material = new GlossyMaterial(brown_sacco_data);
    manager->add_material(brown_sacco_material);

    MaterialData glass_contents_data;
    glass_contents_data.refractive_index = 1.;
    BeersLawMedium* glass_contents_medium =
        new BeersLawMedium(vec3(0.2, 0.2, 0.3) * 0, vec3(2.7, 1., 1.1) * 0.03, vec3(0));
    glass_contents_data.medium = glass_contents_medium;
    TransparentMaterial* glass_contents_material = new TransparentMaterial(glass_contents_data);
    manager->add_material(glass_contents_material);

    MaterialData wall_data;
    wall_data.albedo_map = create_value_map_3D("./maps/plaster_wall.map", 1, -1);
    DiffuseMaterial* wall_material = new DiffuseMaterial(wall_data);
    manager->add_material(wall_material);

    //BeersLawMedium* glass_medium = new BeersLawMedium(vec3(0), (vec3(1,1,1) - colors::BLUE) * 0.0, vec3(0));

    // Walls
    ObjectUnion* floor_obj =
        load_object_model("./models/realistic_room/floor.obj", floor_material, false, false, vec3(0), 0);
    ObjectUnion* roof_obj =
        load_object_model("./models/realistic_room/roof.obj", wall_material, false, false, vec3(0), 0);
    ObjectUnion* wall_obj =
        load_object_model("./models/realistic_room/left_wall.obj", wall_material, false, false, vec3(0), 0);
    ObjectUnion* wall1_obj =
        load_object_model("./models/realistic_room/back_wall.obj", wall_material, false, false, vec3(0), 0);
    ObjectUnion* wall2_obj =
        load_object_model("./models/realistic_room/right_wall.obj", wall_material, false, false, vec3(0), 0);

    ObjectUnion* left_moulding_obj = load_object_model("./models/realistic_room/left_moulding.obj",
                                                       white_diffuse_material, false, false, vec3(0), 0);
    ObjectUnion* back_moulding_obj = load_object_model("./models/realistic_room/back_moulding.obj",
                                                       white_diffuse_material, false, false, vec3(0), 0);
    ObjectUnion* right_moulding_obj = load_object_model("./models/realistic_room/right_moulding.obj",
                                                        white_diffuse_material, false, false, vec3(0), 0);

    // Dining room
    /*
    ObjectUnion* table_obj = load_object_model("./models/realistic_room/table.obj", table_material, false, false, vec3(0), 0);
    ObjectUnion* chair_obj = load_object_model("./models/realistic_room/chair.obj", chair_material, false, false, vec3(0), 0);
    ObjectUnion* chair1_obj = load_object_model("./models/realistic_room/chair1.obj", chair_material, false, false, vec3(0), 0);
    ObjectUnion* chair2_obj = load_object_model("./models/realistic_room/chair2.obj", chair_material, false, false, vec3(0), 0);
    ObjectUnion* chair3_obj = load_object_model("./models/realistic_room/chair3.obj", chair_material, false, false, vec3(0), 0);
    ObjectUnion* glass_obj = load_object_model("./models/realistic_room/glass.obj", glass_material, false, false, vec3(0), 0);
    ObjectUnion* glass1_obj = load_object_model("./models/realistic_room/glass1.obj", glass_material, false, false, vec3(0), 0);
    ObjectUnion* glass2_obj = load_object_model("./models/realistic_room/glass2.obj", glass_material, false, false, vec3(0), 0);
    ObjectUnion* glass3_obj = load_object_model("./models/realistic_room/glass3.obj", glass_material, false, false, vec3(0), 0);
    ObjectUnion* painting_ramen_obj = load_object_model("./models/realistic_room/painting_ramen.obj", painting_ramen_material, false, false, vec3(0), 0);
    */

    // Living room
    ObjectUnion* sacco_obj =
        load_object_model("./models/realistic_room/sacco.obj", sacco_material, true, false, vec3(0), 0);
    ObjectUnion* sacco1_obj =
        load_object_model("./models/realistic_room/sacco1.obj", brown_sacco_material, true, false, vec3(0), 0);
    ObjectUnion* glass4_obj =
        load_object_model("./models/realistic_room/glass4.obj", glass_material, false, false, vec3(0), 0);
    ObjectUnion* glass5_obj =
        load_object_model("./models/realistic_room/glass5.obj", glass_material, false, false, vec3(0), 0);
    ObjectUnion* caraffe_obj =
        load_object_model("./models/realistic_room/water_caraffe.obj", glass_material, true, false, vec3(0), 0);
    ObjectUnion* caraffe_contents_obj =
        load_object_model("./models/realistic_room/caraffe_contents.obj", whiskey_material, true, false, vec3(0), 0);
    ObjectUnion* mini_table_obj =
        load_object_model("./models/realistic_room/mini_table.obj", table_material, false, false, vec3(0), 0);
    ObjectUnion* contents_obj =
        load_object_model("./models/realistic_room/glass_contents.obj", whiskey_material, false, false, vec3(0), 0);
    ObjectUnion* contents1_obj =
        load_object_model("./models/realistic_room/glass_contents1.obj", whiskey_material, false, false, vec3(0), 0);
    //ObjectUnion* bunny_obj = load_object_model("./models/realistic_room/bunny.obj", glass_material, false, false, vec3(0), 0);
    /*
    */

    // Lamp
    ObjectUnion* light_source =
        load_object_model("./models/realistic_room/roof_lamp.obj", light_source_material, false, false, vec3(0), 0);

    int number_of_objects = 18;
    Object** objects = new Object* [number_of_objects] {
        floor_obj, roof_obj, wall_obj, wall1_obj, wall2_obj, left_moulding_obj, back_moulding_obj, right_moulding_obj,
            sacco_obj, sacco1_obj, glass4_obj, contents_obj, glass5_obj, contents1_obj, caraffe_obj,
            caraffe_contents_obj, mini_table_obj, light_source
    };

    ScatteringMediumHomogenous* background_medium =
        new ScatteringMediumHomogenous(vec3(0.), (colors::WHITE) *0.0, vec3(0));

    vec3 camera_position = vec3(-0.7, 1.7, -0.16);
    vec3 viewing_direction = vec3(0.7, -0.3, -0.8);
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
