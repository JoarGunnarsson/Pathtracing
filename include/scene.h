#ifndef SCENE_H
#define SCENE_H

#include <unordered_map>
#include <variant>

#include "utils.h"
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

struct AtrousParamters {
    int iterations;
    double sigma_rt;
    double sigma_x;
    double sigma_n;
};

struct MedianFilterParamters {
    int kernel_size;
    double threshold;
};

struct DenoisingTask {
    std::string mode;
    std::variant<AtrousParamters, MedianFilterParamters> params;
};

struct Scene {
    Object const* const* objects;
    int number_of_objects;
    Camera camera;
    PointerManager* pointer_manager;
    Medium const* medium;
    vec3 background_color;
};

struct SceneStore {
    std::unordered_map<std::string, ValueMap1D const*> valuemap1d_store;
    std::unordered_map<std::string, ValueMap3D const*> valuemap3d_store;
    std::unordered_map<std::string, Medium const*> medium_store;
    std::unordered_map<std::string, Material const*> material_store;
    std::unordered_map<std::string, Object const*> object_store;
};

void load_settings(const std::string& file_path);
std::vector<DenoisingTask> load_denoising_settings(const std::string& file_path);
Scene load_scene(const std::string& file_path);

#endif
