#ifndef SCENE_H
#define SCENE_H

#include <unordered_map>

#include "objects.h"
#include "materials.h"
#include "medium.h"
#include "camera.h"

class PointerManager {
    std::vector<ValueMap const*> valuemaps;
    std::vector<Medium const*> media;
    std::vector<Material const*> materials;

  public:
    PointerManager() {}
    ~PointerManager();

    template<typename T>
    void clear_pointer_array(T pointer_array);

    // PointerManager takes owenership of pointers passed to the following methods.
    void add_valuemap(ValueMap const* const map);
    void add_medium(Medium const* const medium);
    void add_material(Material const* const material);
};

struct Scene {
    Object const* const* objects;
    int number_of_objects;
    Camera camera;
    PointerManager* pointer_manager;
    Medium const* medium;
};

struct SceneStore {
    std::unordered_map<std::string, ValueMap1D const*> valuemap1d_store;
    std::unordered_map<std::string, ValueMap3D const*> valuemap3d_store;
    std::unordered_map<std::string, Medium const*> medium_store;
    std::unordered_map<std::string, Material const*> material_store;
    std::unordered_map<std::string, Object const*> object_store;
};

void load_settings(const std::string& file_path);
Scene load_scene(const std::string& file_path);
Scene create_scene();

#endif
