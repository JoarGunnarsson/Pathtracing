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

template<typename T>
void PointerManager::clear_pointer_array(T pointer_array) {
    for (size_t i = 0; i < pointer_array.size(); ++i) {
        delete pointer_array[i];
    }
}

PointerManager::~PointerManager() {
    clear_pointer_array(valuemaps);
    clear_pointer_array(media);
    clear_pointer_array(materials);
}

void PointerManager::add_valuemap(ValueMap* map) {
    valuemaps.push_back(map);
};

void PointerManager::add_medium(Medium* medium) {
    media.push_back(medium);
};

void PointerManager::add_material(Material* material) {
    materials.push_back(material);
};

template<typename T>
void load_one_setting(const json& data, const std::string& key, T& setting) {
    if (data.contains(key)) {
        setting = data[key];
    }
}

void load_settings(const std::string& file_path) {
    std::ifstream settings_file(file_path);
    if (!settings_file.is_open()) {
        throw std::runtime_error("Could not open file '" + file_path + "'");
    }
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

void require_field(const json& data, const std::string& key) {
    if (!data.contains(key)) {
        throw std::runtime_error("JSON structure missing key '" + key + "'");
    }
}

template<typename T>
bool key_in_map(const T& map, const std::string& key) {
    return map.find(key) != map.end();
}

template<typename T>
void require_unique_key(const T& map, const std::string& key, const std::string& key_type) {
    if (key_in_map(map, key)) {
        throw std::runtime_error("Duplicate " + key_type + " entry '" + key + "' found");
    }
}

template<typename T>
void require_key_exists(const T& map, const std::string& key, const std::string& key_type) {
    if (!key_in_map(map, key)) {
        throw std::runtime_error(key_type + " entry '" + key + "' does not exist");
    }
}

vec3 get_vec3_param(const json& data, const std::string& key) {
    json vector_data = data[key];
    if (vector_data.size() != 3) {
        throw std::runtime_error("Data array " + vector_data.dump() + " should contain 3 elements, but contains " +
                                 std::to_string(vector_data.size()));
    }
    return vec3(vector_data[0], vector_data[1], vector_data[2]);
}

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
        require_key_exists(store.valuemap3d_store, albedo_map, "ValueMap3D");
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
        require_key_exists(store.valuemap3d_store, emission_color_map, "ValueMap3D");
        material_data.emission_color_map = store.valuemap3d_store.find(emission_color_map)->second;
    }
    if (parameters.contains("light_intensity_map")) {
        std::string light_intensity_map = parameters["light_intensity_map"];
        require_key_exists(store.valuemap1d_store, light_intensity_map, "ValueMap1D");
        material_data.light_intensity_map = store.valuemap1d_store.find(light_intensity_map)->second;
    }
    if (parameters.contains("is_dielectric")) {
        material_data.is_dielectric = (bool) parameters["is_dielectric"];
    }
    if (parameters.contains("roughness_map")) {
        std::string roughness_map = parameters["roughness_map"];
        require_key_exists(store.valuemap1d_store, roughness_map, "ValueMap1D");
        material_data.roughness_map = store.valuemap1d_store.find(roughness_map)->second;
    }
    if (parameters.contains("is_light_source")) {
        material_data.is_light_source = (bool) parameters["is_light_source"];
    }
    if (parameters.contains("medium")) {
        std::string medium = parameters["medium"];
        require_key_exists(store.medium_store, medium, "Medium");
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
    else if (material_type == "MetallicMicrofacet") {
        return new MetallicMicrofacetMaterial(material_data);
    }
    else if (material_type == "ReflectiveMicrofacet") {
        return new ReflectiveMicrofacetMaterial(material_data);
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
    require_key_exists(store.material_store, parameters["material"], "Material");
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

Camera load_camera(const json& data) {
    require_field(data, "camera");
    json camera_data = data["camera"];

    require_field(camera_data, "camera_position");
    require_field(camera_data, "viewing_direction");
    require_field(camera_data, "screen_y_vector");

    vec3 camera_position = get_vec3_param(camera_data, "camera_position");
    vec3 viewing_direction = get_vec3_param(camera_data, "viewing_direction");
    vec3 screen_y_vector = get_vec3_param(camera_data, "screen_y_vector");
    return Camera(camera_position, viewing_direction, screen_y_vector);
}

void populate_scene_store(json& scene_data, SceneStore& store, PointerManager* manager) {
    require_field(scene_data, "valuemaps");
    for (json& element : scene_data["valuemaps"]) {
        require_field(element, "name");
        std::string name = element["name"];

        require_field(element, "type");
        std::string type = element["type"];

        if (type == "ValueMap1D") {
            require_unique_key(store.valuemap1d_store, name, "valuemap");
            ValueMap1D* map = load_valuemap1d(element);
            store.valuemap1d_store[name] = map;
            manager->add_valuemap(map);
        }
        else if (type == "ValueMap3D") {
            require_unique_key(store.valuemap3d_store, name, "valuemap");
            ValueMap3D* map = load_valuemap3d(element);
            store.valuemap3d_store[name] = map;
            manager->add_valuemap(map);
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
            require_unique_key(store.medium_store, name, "medium");
            Medium* medium = load_medium(element, store);
            store.medium_store[name] = medium;
            manager->add_medium(medium);
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
            require_unique_key(store.material_store, name, "material");
            Material* material = load_material(element, store);
            store.material_store[name] = material;
            manager->add_material(material);
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
            require_unique_key(store.object_store, name, "object");
            Object* object = load_object(element, store);
            store.object_store[name] = object;
        }
        else {
            throw std::runtime_error("Found invalid type in object array: " + type);
        }
    }
}

Scene load_scene(const std::string& file_path) {
    std::ifstream scene_file(file_path);
    if (!scene_file.is_open()) {
        throw std::runtime_error("Could not open file '" + file_path + "'");
    }
    json scene_data = json::parse(scene_file);

    PointerManager* manager = new PointerManager();
    SceneStore store;
    populate_scene_store(scene_data, store, manager);

    int number_of_objects = store.object_store.size();
    Object** objects = new Object*[number_of_objects];
    int i = 0;
    for (auto& [key, value] : store.object_store) {
        objects[i] = value;
        i++;
    }

    require_field(scene_data, "background_medium");
    require_key_exists(store.medium_store, scene_data["background_medium"], "Medium");
    Medium* background_medium = store.medium_store.find(scene_data["background_medium"])->second;

    Camera camera = load_camera(scene_data);

    Scene scene;
    scene.objects = objects;
    scene.camera = camera;
    scene.number_of_objects = number_of_objects;
    scene.pointer_manager = manager;
    scene.medium = background_medium;
    return scene;
}
