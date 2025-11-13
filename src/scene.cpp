#include <nlohmann/json.hpp>
#include <fstream>
#include <string>
#include <unordered_map>

#include "scene.h"
#include "objects.h"
#include "objectunion.h"
#include "materials.h"
#include "camera.h"
#include "medium.h"
#include "valuemap.h"

using json = nlohmann::json;

template<typename T> void load_one_setting(const json& data, const std::string& key, T& setting) {
    if (data.contains(key)) {
        setting = data[key];
    }
}

void require_field(const json& data, const std::string& key) {
    if (!data.contains(key)) {
        throw std::runtime_error("JSON structure missing key '" + key + "'");
    }
}

void load_settings(const std::string& file_path) {
    std::ifstream settings_file(file_path);
    json data = json::parse(settings_file);

    load_one_setting(data, "WIDTH", constants::WIDTH);
    load_one_setting(data, "HEIGHT", constants::HEIGHT);
    load_one_setting(data, "samples_per_pixel", constants::samples_per_pixel);
    load_one_setting(data, "max_recursion_depth", constants::max_recursion_depth);
    load_one_setting(data, "force_tracing_limit", constants::force_tracing_limit);
    load_one_setting(data, "air_refractive_index", constants::air_refractive_index);
    load_one_setting(data, "enable_next_event_estimation", constants::enable_next_event_estimation);
    load_one_setting(data, "enable_anti_aliasing", constants::enable_anti_aliasing);
    load_one_setting(data, "enable_denoising", constants::enable_denoising);
    load_one_setting(data, "denoising_iterations", constants::denoising_iterations);
    load_one_setting(data, "sigma_rt", constants::sigma_rt);
    load_one_setting(data, "sigma_x", constants::sigma_x);
    load_one_setting(data, "sigma_n", constants::sigma_n);
}

vec3 get_vec3_param(const json& data, const std::string& key) {
    json vector_data = data[key];
    if (vector_data.size() != 3) {
        throw std::runtime_error("Data array " + vector_data.dump() + " should contain 3 elements, but contains " +
                                 std::to_string(vector_data.size()));
    }
    return vec3(vector_data[0], vector_data[1], vector_data[2]);
}

template<typename T> bool key_in_map(const T& map, const std::string& key) {
    return map.find(key) != map.end();
}

template<typename T> void require_unique_key(const T& map, const std::string& key) {
    if (key_in_map(map, key)) {
        throw std::runtime_error("Duplicate entry '" + key + "' found");
    }
}

template<typename T> void require_key_exists(const T& map, const std::string& key) {
    if (!key_in_map(map, key)) {
        throw std::runtime_error("Entry '" + key + "' does not exist");
    }
}

// TODO: Could template this.
ValueMap1D* load_valuemap1d(const json& data) {
    require_field(data, "parameters");
    json parameters = data["parameters"];

    if (parameters.contains("data")) {
        json map_data = parameters["data"];
        if (map_data.size() != 1) {
            throw std::runtime_error("Data array " + map_data.dump() + " should contain 1 element, but contains " +
                                     std::to_string(map_data.size()));
        }
        double value = map_data[0];
        return new ValueMap1D(value);
    }
    else if (parameters.contains("file")) {
        return create_value_map_1D(parameters["file"], 1, -1);
    }
    else {
        throw std::runtime_error("ValueMap must contain either 'data' or 'file'");
    }
}

ValueMap3D* load_valuemap3d(const json& data) {
    require_field(data, "parameters");
    json parameters = data["parameters"];

    if (parameters.contains("data")) {
        vec3 colour_data = get_vec3_param(parameters, "data");
        return new ValueMap3D(colour_data);
    }
    else if (parameters.contains("file")) {
        return create_value_map_3D(parameters["file"], 1, -1);
    }
    else {
        throw std::runtime_error("ValueMap must contain either 'data' or 'file'");
    }
}

Medium* load_medium(const json& data, const SceneStore& store) {
    require_field(data, "parameters");
    json parameters = data["parameters"];

    require_field(parameters, "scattering_albedo");
    require_field(parameters, "absorption_albedo");
    require_field(parameters, "emission_coefficient");

    vec3 scattering_albedo = get_vec3_param(parameters, "scattering_albedo");
    vec3 absorption_albedo = get_vec3_param(parameters, "absorption_albedo");
    vec3 emission_coefficient = get_vec3_param(parameters, "emission_coefficient");

    require_field(data, "subtype");
    std::string medium_type = data["subtype"];
    if (medium_type == "BeersLawMedium") {
        return new BeersLawMedium(scattering_albedo, absorption_albedo, emission_coefficient);
    }
    else if (medium_type == "HomogenousScatteringMedium") {
        return new HomogenousScatteringMedium(scattering_albedo, absorption_albedo, emission_coefficient);
    }
    else {
        throw std::runtime_error(medium_type + " is not a valid medium type");
    }
}

Material* load_material(const json& data, const SceneStore& store) {
    MaterialData material_data;

    require_field(data, "parameters");
    json parameters = data["parameters"];

    if (parameters.contains("albedo_map")) {
        std::string albedo_map = parameters["albedo_map"];
        require_key_exists(store.valuemap3d_store, albedo_map);
        material_data.albedo_map = store.valuemap3d_store.find(albedo_map)->second;
    }
    if (parameters.contains("refractive_index")) {
        material_data.refractive_index = (double) parameters["refractive_index"];
    }
    if (parameters.contains("extinction_coefficient")) {
        material_data.extinction_coefficient = (double) parameters["extinction_coefficient"];
    }
    if (parameters.contains("emission_color_map")) {
        std::string emission_color_map = parameters["emission_color_map"];
        require_key_exists(store.valuemap3d_store, emission_color_map);
        material_data.emission_color_map = store.valuemap3d_store.find(emission_color_map)->second;
    }
    if (parameters.contains("light_intensity_map")) {
        std::string light_intensity_map = parameters["light_intensity_map"];
        require_key_exists(store.valuemap1d_store, light_intensity_map);
        material_data.light_intensity_map = store.valuemap1d_store.find(light_intensity_map)->second;
    }
    if (parameters.contains("is_dielectric")) {
        material_data.is_dielectric = (bool) parameters["is_dielectric"];
    }
    if (parameters.contains("roughness_map")) {
        std::string roughness_map = parameters["roughness_map"];
        require_key_exists(store.valuemap1d_store, roughness_map);
        material_data.roughness_map = store.valuemap1d_store.find(roughness_map)->second;
    }
    if (parameters.contains("is_light_source")) {
        material_data.is_light_source = (bool) parameters["is_light_source"];
    }
    if (parameters.contains("medium")) {
        std::string medium = parameters["medium"];
        require_key_exists(store.medium_store, medium);
        material_data.medium = store.medium_store.find(medium)->second;
    }

    require_field(data, "subtype");
    std::string material_type = data["subtype"];
    if (material_type == "Diffuse") {
        return new DiffuseMaterial(material_data);
    }
    else if (material_type == "Reflective") {
        return new ReflectiveMaterial(material_data);
    }
    else if (material_type == "Transparent") {
        return new TransparentMaterial(material_data);
    }
    else if (material_type == "Glossy") {
        return new GlossyMaterial(material_data);
    }
    else if (material_type == "Metallic") {
        return new MetallicMicrofacet(material_data);
    }
    else if (material_type == "TransparentMicrofacet") {
        return new TransparentMicrofacetMaterial(material_data);
    }
    else {
        throw std::runtime_error(material_type + " is not a valid material type");
    }
}

Object* load_object(const json& data, const SceneStore& store) {
    require_field(data, "subtype");
    std::string object_type = data["subtype"];

    require_field(data, "parameters");
    json parameters = data["parameters"];

    require_field(parameters, "material");
    require_key_exists(store.material_store, parameters["material"]);
    Material* material = store.material_store.find(parameters["material"])->second;

    if (object_type == "Sphere") {
        require_field(parameters, "position");
        require_field(parameters, "radius");

        vec3 position = get_vec3_param(parameters, "position");
        double radius = parameters["radius"];
        return new Sphere(position, radius, material);
    }
    else if (object_type == "Plane") {
        require_field(parameters, "position");
        require_field(parameters, "v1");
        require_field(parameters, "v2");

        vec3 position = get_vec3_param(parameters, "position");
        vec3 v1 = get_vec3_param(parameters, "v1");
        vec3 v2 = get_vec3_param(parameters, "v2");
        return new Plane(position, v1, v2, material);
    }
    else if (object_type == "Rectangle") {
        require_field(parameters, "position");
        require_field(parameters, "v1");
        require_field(parameters, "v2");
        require_field(parameters, "L1");
        require_field(parameters, "L2");

        vec3 position = get_vec3_param(parameters, "position");
        vec3 v1 = get_vec3_param(parameters, "v1");
        vec3 v2 = get_vec3_param(parameters, "v2");
        double L1 = parameters["L1"];
        double L2 = parameters["L2"];
        return new Rectangle(position, v1, v2, L1, L2, material);
    }
    else if (object_type == "Triangle") {
        require_field(parameters, "p1");
        require_field(parameters, "p2");
        require_field(parameters, "p3");

        vec3 p1 = get_vec3_param(parameters, "p1");
        vec3 p2 = get_vec3_param(parameters, "p2");
        vec3 p3 = get_vec3_param(parameters, "p3");
        return new Triangle(p1, p2, p3, material);
    }
    else if (object_type == "ObjectUnion") {
        require_field(parameters, "file");
        require_field(parameters, "enable_smooth_shading");
        require_field(parameters, "move_object");

        std::string file_name = parameters["file"];
        bool enable_smooth_shading = parameters["enable_smooth_shading"];
        bool move_object = parameters["move_object"];
        vec3 center = vec3(0);
        double size = 0;

        if (move_object) {
            require_field(parameters, "center");
            require_field(parameters, "size");

            vec3 center = get_vec3_param(parameters, "center");
            double size = parameters["size"];
        }
        return load_object_model(file_name, material, enable_smooth_shading, move_object, center, size);
    }
    else {
        throw std::runtime_error(object_type + " is not a valid object type");
    }
}

Camera* load_camera(const json& data) {
    require_field(data, "camera");
    json camera_data = data["camera"];

    require_field(camera_data, "camera_position");
    require_field(camera_data, "viewing_direction");
    require_field(camera_data, "screen_y_vector");

    vec3 camera_position = get_vec3_param(camera_data, "camera_position");
    vec3 viewing_direction = get_vec3_param(camera_data, "viewing_direction");
    vec3 screen_y_vector = get_vec3_param(camera_data, "screen_y_vector");
    return new Camera(camera_position, viewing_direction, screen_y_vector);
}

void populate_scene_store(json& scene_data, SceneStore& store) {
    // TODO: Need to iterate in a specific order: First value maps, then media, then materials, then objects.
    require_field(scene_data, "valuemaps");
    for (json& element : scene_data["valuemaps"]) {
        require_field(element, "name");
        std::string name = element["name"];

        require_field(element, "type");
        std::string type = element["type"];

        if (type == "ValueMap1D") {
            require_unique_key(store.valuemap1d_store, name);
            store.valuemap1d_store[name] = load_valuemap1d(element);
        }
        else if (type == "ValueMap3D") {
            require_unique_key(store.valuemap3d_store, name);
            store.valuemap3d_store[name] = load_valuemap3d(element);
        }
        else {
            throw std::runtime_error("Found invalid type in ValueMap array: " + type);
        }
    }
    require_field(scene_data, "media");
    for (json& element : scene_data["media"]) {
        require_field(element, "name");
        std::string name = element["name"];

        require_field(element, "type");
        std::string type = element["type"];

        if (type == "Medium") {
            require_unique_key(store.medium_store, name);
            store.medium_store[name] = load_medium(element, store);
        }
        else {
            throw std::runtime_error("Found invalid type in medium array: " + type);
        }
    }
    for (json& element : scene_data["materials"]) {
        require_field(element, "name");
        std::string name = element["name"];

        require_field(element, "type");
        std::string type = element["type"];

        if (type == "Material") {
            require_unique_key(store.material_store, name);
            store.material_store[name] = load_material(element, store);
        }
        else {
            throw std::runtime_error("Found invalid type in material array: " + type);
        }
    }
    for (json& element : scene_data["objects"]) {
        require_field(element, "name");
        std::string name = element["name"];

        require_field(element, "type");
        std::string type = element["type"];

        if (type == "Object") {
            require_unique_key(store.object_store, name);
            store.object_store[name] = load_object(element, store);
        }
        else {
            throw std::runtime_error("Found invalid type in object array: " + type);
        }
    }
}
Scene load_scene(const std::string& file_path) {
    // TODO: Add managers for all maps, so that memory cleanup works properly.
    std::ifstream scene_file(file_path);
    json scene_data = json::parse(scene_file);

    SceneStore store;
    populate_scene_store(scene_data, store);

    MaterialManager* manager = new MaterialManager();

    int number_of_objects = store.object_store.size();
    Object** objects = new Object*[number_of_objects];
    int i = 0;
    for (auto& [key, value] : store.object_store) {
        objects[i] = value;
        i++;
    }

    require_field(scene_data, "background_medium");
    require_key_exists(store.medium_store, scene_data["background_medium"]);
    Medium* background_medium = store.medium_store.find(scene_data["background_medium"])->second;

    Camera* camera = load_camera(scene_data);

    Scene scene;
    scene.objects = objects;
    scene.camera = camera;
    scene.number_of_objects = number_of_objects;
    scene.material_manager = manager;
    scene.medium = background_medium;
    return scene;
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
    MetallicMicrofacetMaterial* gold_material = new MetallicMicrofacetMaterial(gold_data);
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
    HomogenousScatteringMedium* scattering_glass_medium =
        new HomogenousScatteringMedium(colors::WHITE, vec3(2.7, 1., 1.1) * 0.0, vec3(0));
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

    HomogenousScatteringMedium* background_medium =
        new HomogenousScatteringMedium(vec3(0.), (colors::WHITE) *0.0, vec3(0));

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
