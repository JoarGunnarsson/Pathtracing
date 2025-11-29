#ifndef SCENE_H
#define SCENE_H

#include <unordered_map>

#include "objects.h"
#include "materials.h"
#include "medium.h"
#include "camera.h"

class PointerManager {
    std::vector<ValueMap*> valuemaps;
    std::vector<Medium*> media;
    std::vector<Material*> materials;

  public:
    PointerManager() {}
    ~PointerManager();

    template<typename T>
    void clear_pointer_array(T pointer_array);

    // PointerManager takes owenership of pointers passed to the following methods.
    void add_valuemap(ValueMap* map);
    void add_medium(Medium* medium);
    void add_material(Material* material);
};

struct Scene {
    Object** objects;
    size_t number_of_objects;
    Camera camera;
    PointerManager* pointer_manager;
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
