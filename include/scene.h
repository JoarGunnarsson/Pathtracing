#ifndef SCENE_H
#define SCENE_H

#include <unordered_map>

#include "objects.h"
#include "materials.h"
#include "medium.h"
#include "camera.h"

struct Scene {
    Object** objects;
    int number_of_objects;
    Camera* camera;
    MaterialManager* material_manager;
    Medium* medium;
};

struct SceneStore {
    std::unordered_map<std::string, ValueMap1D*> valuemap1d_store;
    std::unordered_map<std::string, ValueMap3D*> valuemap3d_store;
    std::unordered_map<std::string, Medium*> medium_store;
    std::unordered_map<std::string, Material*> material_store;
    std::unordered_map<std::string, Object*> object_store;
};

void load_settings(const std::string& file_path);
Scene load_scene(const std::string& file_path);
Scene create_scene();

#endif
