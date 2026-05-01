#ifndef OBJECTUNION_H
#define OBJECTUNION_H

#include <fstream>
#include "constants.h"
#include "vec3.h"
#include "utils.h"
#include "materials.h"
#include "objects.h"
#include "bvh.h"

struct ObjectTransform {
    bool move_object;
    bool rotate_object;
    bool scale_object;
    vec3 center;
    vec3 orientation;
    double size;
};

class ObjectUnion : public Object {
  public:
    ObjectUnion(Object** _objects, const int _number_of_objects, const bool construct_BVH = false);
    ~ObjectUnion();

    bool is_light_source() const override;
    bool allow_direct_light(const vec3& intersection_point, const int primitive_ID) const override;
    Material const* get_material(const int primitive_ID) const override;
    vec3 eval(const Hit& hit, const vec3& outgoing_vector) const override;
    BrdfData sample(const Hit& hit) const override;
    double brdf_pdf(const vec3& outgoing_vector, const Hit& hit) const override;
    vec3 get_light_emittance(const Hit& hit) const override;
    bool find_closest_object_hit(Hit& hit, Ray& ray) const override;
    vec3 get_normal_vector(const vec3& surface_point, const int primitive_ID) const override;
    int sample_random_primitive_index() const;
    vec3 generate_random_surface_point() const override;
    double light_pdf(const vec3& surface_point, const vec3& intersection_point, const int primitive_id) const override;
    vec3 random_light_point(const vec3& intersection_point, double& inverse_PDF) const override;

  private:
    bool use_BVH;
    bool contains_light_source = false;
    int number_of_objects;
    int* light_source_conversion_indices;
    double* cumulative_area;
    size_t number_of_light_sources;
    BVH::BoundingVolumeHierarchy bvh;
    Object** objects;
};

struct DataSizes {
    size_t num_vertices = 0;
    size_t num_vertex_UVs = 0;
    size_t num_vertex_normals = 0;
    size_t num_triangles = 0;
};

int number_of_char_occurances(const std::string& line, const char character);
std::string get_nth_word(const std::string& line, const char delimiter, const int n);
DataSizes get_vertex_data_sizes(const std::string& file_name);
void populate_vertex_arrays(const std::string& file_name, vec3* vertex_array, vec3* vertex_UV_array,
                            vec3* vertex_normal_array);
double maximum_distance(const vec3& center, const vec3* vertex_array, const size_t number_of_vertices);
void change_vectors(const ObjectTransform& transform, vec3* const vertex_array, const size_t number_of_vertices);

struct PopulateVertexVectorData {
    bool v_success = false;
    bool uv_success = false;
    bool n_success = false;
    const std::string vertex_data;
    vec3 v;
    vec3 uv;
    vec3 n;
    vec3 const* vertex_array;
    vec3 const* vertex_UV_array;
    vec3 const* vertex_normal_array;

    PopulateVertexVectorData(const std::string& data, const vec3* vertex_array, const vec3* vertex_UV_array,
                             const vec3* vertex_normal_array) :
        vertex_data(data),
        vertex_array(vertex_array),
        vertex_UV_array(vertex_UV_array),
        vertex_normal_array(vertex_normal_array) {}
};

void populate_vertex_vectors(PopulateVertexVectorData& args);

struct TriangleConstructionArgs {
    const bool enable_smooth_shading;
    const int idx1;
    const int idx2;
    const int idx3;
    vec3 const* vertex_array;
    vec3 const* vertex_UV_array;
    vec3 const* vertex_normal_array;
    const std::string triangle_data;
    Material const* material;

    TriangleConstructionArgs(const std::string& data, const int idx1, const int idx2, const int idx3,
                             Material const* material, const vec3* vertex_array, const vec3* vertex_UV_array,
                             const vec3* vertex_normal_array, const bool enable_smooth_shading) :
        enable_smooth_shading(enable_smooth_shading),
        idx1(idx1),
        idx2(idx2),
        idx3(idx3),
        vertex_array(vertex_array),
        vertex_UV_array(vertex_UV_array),
        vertex_normal_array(vertex_normal_array),
        triangle_data(data),
        material(material) {}
};

struct TriangleCreationResult {
    bool success = false;
    Triangle* triangle;
};

TriangleCreationResult construct_triangle(TriangleConstructionArgs& args);
int populate_triangle_array(std::string file_name, vec3* vertex_array, vec3* vertex_UV_array, vec3* vertex_normal_array,
                            Object** triangle_array, Material const* material, const bool enable_smooth_shading);
ObjectUnion* load_object_model(std::string file_name, Material const* material, const bool enable_smooth_shading,
                               const ObjectTransform& transform);

#endif
