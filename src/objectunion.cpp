#include "objectunion.h"

ObjectUnion::ObjectUnion(Object** _objects, const int _number_of_objects, const bool construct_BVH) : Object() {
    objects = _objects;
    number_of_objects = _number_of_objects;
    number_of_light_sources = 0;
    area = 0;

    for (int i = 0; i < number_of_objects; i++) {
        area += objects[i]->area;
        if (objects[i]->is_light_source()) {
            number_of_light_sources++;
        }
    }

    size_t area_array_size = static_cast<size_t>(number_of_objects);
    cumulative_area = new double[area_array_size];

    light_source_conversion_indices = new int[number_of_light_sources];
    int j = 0;
    for (int i = 0; i < number_of_objects; i++) {
        if (!objects[i]->is_light_source()) {
            continue;
        }
        cumulative_area[j] = objects[i]->area;
        if (j != 0) {
            cumulative_area[j] += cumulative_area[j - 1];
        }
        light_source_conversion_indices[j] = i;
        j++;
    }

    use_BVH = construct_BVH;
    if (construct_BVH) {
        bvh = BVH::BoundingVolumeHierarchy(_objects, _number_of_objects, 6);
    }

    for (int i = 0; i < number_of_objects; i++) {
        objects[i]->primitive_ID = static_cast<int>(i);
        if (objects[i]->is_light_source()) {
            contains_light_source = true;
        }
    }
}

ObjectUnion::~ObjectUnion() {
    for (int i = 0; i < number_of_objects; i++) {
        delete objects[i];
    }
    delete[] objects;
    delete[] cumulative_area;
    delete[] light_source_conversion_indices;
}

Material* ObjectUnion::get_material(const int primitive_ID) const {
    return objects[primitive_ID]->material;
}

bool ObjectUnion::is_light_source() const {
    return contains_light_source;
}

vec3 ObjectUnion::eval(const Hit& hit, const vec3& outgoing_vector) const {
    return objects[hit.primitive_ID]->eval(hit, outgoing_vector);
}

BrdfData ObjectUnion::sample(const Hit& hit) const {
    return objects[hit.primitive_ID]->sample(hit);
}

double ObjectUnion::brdf_pdf(const vec3& outgoing_vector, const Hit& hit) const {
    return objects[hit.primitive_ID]->brdf_pdf(outgoing_vector, hit);
}

vec3 ObjectUnion::get_light_emittance(const Hit& hit) const {
    return objects[hit.primitive_ID]->get_light_emittance(hit);
}

bool ObjectUnion::find_closest_object_hit(Hit& hit, Ray& ray) const {
    if (use_BVH) {
        return bvh.intersect(hit, ray);
    }

    bool success = find_closest_hit(hit, ray, objects, number_of_objects);
    hit.primitive_ID = hit.intersected_object_index;
    hit.intersected_object_index = -1;
    return success;
}

vec3 ObjectUnion::get_normal_vector(const vec3& surface_point, const int primitive_ID) const {
    return objects[primitive_ID]->get_normal_vector(surface_point, primitive_ID);
}

int ObjectUnion::sample_random_primitive_index() const {
    double random_area_split = random_uniform(0, area);
    int max = static_cast<int>(number_of_light_sources) - 1;
    int min = 0;
    int index;

    if (cumulative_area[0] >= random_area_split) {
        return light_source_conversion_indices[0];
    }

    while (min <= max) {
        index = (max - min) / 2 + min;

        if (cumulative_area[index] < random_area_split) {
            min = index + 1;
        }
        else if (cumulative_area[index] == random_area_split ||
                 (cumulative_area[index] >= random_area_split && cumulative_area[index - 1] < random_area_split)) {
            break;
        }
        else {
            max = index - 1;
        }
    }

    return light_source_conversion_indices[index];
}

vec3 ObjectUnion::generate_random_surface_point() const {
    return objects[sample_random_primitive_index()]->generate_random_surface_point();
}

double ObjectUnion::light_pdf(const vec3& surface_point, const vec3& intersection_point, const int primitive_id) const {
    return 1.0 / (cumulative_area[number_of_light_sources - 1] *
                  area_to_angle_PDF_factor(surface_point, intersection_point, primitive_id));
}

vec3 ObjectUnion::random_light_point(const vec3& intersection_point, double& pdf) const {
    int random_index = sample_random_primitive_index();
    vec3 random_point = objects[random_index]->generate_random_surface_point();
    pdf = light_pdf(random_point, intersection_point, random_index);
    return random_point;
}

int number_of_char_occurances(const std::string& line, const char character) {
    int count = 0;
    for (size_t i = 0; i < line.length(); i++) {
        if (line[i] == character) {
            count++;
        }
    }
    return count;
}

std::string get_nth_word(const std::string& line, const char delimiter, const int n) {
    int number_of_words = number_of_char_occurances(line, delimiter);

    if (n > number_of_words) {
        return "";
    }

    size_t start = 0;
    size_t end = 0;
    for (int i = 0; i < n + 1; i++) {
        start = end;
        end = line.find(delimiter, start) + 1;
    }
    return line.substr(start, end - start - 1);
}

DataSizes get_vertex_data_sizes(const std::string& file_name) {
    std::ifstream model_file(file_name);
    if (!model_file.is_open()) {
        throw std::runtime_error("Could not open file '" + file_name + "'");
    }

    std::string line;
    DataSizes nums;

    while (std::getline(model_file, line)) {
        std::string first_word = get_nth_word(line, ' ', 0);

        bool is_vertex = first_word == "v";
        bool is_vertex_UV = first_word == "vt";
        bool is_vertex_normal = first_word == "vn";
        bool is_shape = first_word == "f";

        if (is_vertex) {
            nums.num_vertices++;
        }
        else if (is_vertex_UV) {
            nums.num_vertex_UVs++;
        }
        else if (is_vertex_normal) {
            nums.num_vertex_normals++;
        }
        else if (is_shape) {
            int number_of_spaces = 0;
            for (size_t i = 0; i < line.size(); i++) {
                if (line.substr(i, 1) == " ") {
                    number_of_spaces++;
                }
            }

            bool is_triangle = number_of_spaces == 3;
            bool is_quad = number_of_spaces == 4;
            if (is_triangle) {
                nums.num_triangles++;
            }

            else if (is_quad) {
                nums.num_triangles += 2;
            }
        }
    }
    return nums;
}

void populate_vertex_arrays(const std::string& file_name, vec3* vertex_array, vec3* vertex_UV_array,
                            vec3* vertex_normal_array) {
    int vertex_idx = 0;
    int vertex_UV_idx = 0;
    int vertex_normal_idx = 0;

    std::ifstream model_file(file_name);
    std::string line;
    while (std::getline(model_file, line)) {
        std::string first_word = get_nth_word(line, ' ', 0);
        bool is_vertex = first_word == "v";
        bool is_vertex_UV = first_word == "vt";
        bool is_vertex_normal = first_word == "vn";

        if (is_vertex) {
            double v1 = std::stod(get_nth_word(line, ' ', 1));
            double v2 = std::stod(get_nth_word(line, ' ', 2));
            double v3 = std::stod(get_nth_word(line, ' ', 3));
            vertex_array[vertex_idx] = vec3(v1, v2, v3);
            vertex_idx++;
        }
        else if (is_vertex_UV) {
            double u = std::stod(get_nth_word(line, ' ', 1));
            double v = std::stod(get_nth_word(line, ' ', 2));
            vertex_UV_array[vertex_UV_idx] = vec3(u, v, 0);
            vertex_UV_idx++;
        }
        else if (is_vertex_normal) {
            double n1 = std::stod(get_nth_word(line, ' ', 1));
            double n2 = std::stod(get_nth_word(line, ' ', 2));
            double n3 = std::stod(get_nth_word(line, ' ', 3));
            vertex_normal_array[vertex_normal_idx] = vec3(n1, n2, n3);
            vertex_normal_idx++;
        }
    }
}

double maximum_vector_component(vec3 const* vertex_array, size_t const number_of_vertices, int const dimension) {
    double current_max = vertex_array[0][dimension];
    for (size_t i = 0; i < number_of_vertices; i++) {
        current_max = std::max(vertex_array[i][dimension], current_max);
    }
    return current_max;
}

double minimum_vector_component(vec3 const* vertex_array, size_t const number_of_vertices, int const dimension) {
    double current_min = vertex_array[0][dimension];
    for (size_t i = 0; i < number_of_vertices; i++) {
        current_min = std::min(vertex_array[i][dimension], current_min);
    }
    return current_min;
}

vec3 compute_AABB_midpoint(vec3 const* vertex_array, size_t const number_of_vertices) {
    vec3 max_point, min_point;
    for (int i = 0; i < 3; i++) {
        max_point[i] = maximum_vector_component(vertex_array, number_of_vertices, i);
        min_point[i] = minimum_vector_component(vertex_array, number_of_vertices, i);
    }
    return (max_point + min_point) / 2.0;
}

double maximum_distance(vec3 const& center, vec3 const* vertex_array, size_t const number_of_vertices) {
    double max_distance = 0;
    for (size_t i = 0; i < number_of_vertices; i++) {
        double distance = (vertex_array[i] - center).length();
        if (distance > max_distance) {
            max_distance = distance;
        }
    }
    return max_distance;
}

void change_vectors(const vec3& desired_center, const double desired_size, vec3* vertex_array,
                    const size_t number_of_vertices) {
    vec3 current_center = compute_AABB_midpoint(vertex_array, number_of_vertices);
    double max_distance = maximum_distance(current_center, vertex_array, number_of_vertices);

    for (size_t i = 0; i < number_of_vertices; i++) {
        vertex_array[i] = ((vertex_array[i] - current_center) / max_distance) * desired_size + desired_center;
    }
}

void populate_vertex_vectors(PopulateVertexVectorData& args) {
    std::string v_idx = get_nth_word(args.vertex_data, '/', 0);
    std::string UV_idx = get_nth_word(args.vertex_data, '/', 1);
    std::string n_idx = get_nth_word(args.vertex_data, '/', 2);

    if (v_idx != "") {
        args.v = args.vertex_array[std::stoi(v_idx) - 1];
        args.v_success = true;
    }

    if (UV_idx != "") {
        args.uv = args.vertex_UV_array[std::stoi(UV_idx) - 1];
        args.uv_success = true;
    }

    if (n_idx != "") {
        args.n = args.vertex_normal_array[std::stoi(n_idx) - 1];
        args.n_success = true;
    }
}

TriangleCreationResult construct_triangle(TriangleConstructionArgs& args) {
    std::string v1_data = get_nth_word(args.triangle_data, ' ', args.idx1);
    PopulateVertexVectorData data1 =
        PopulateVertexVectorData(v1_data, args.vertex_array, args.vertex_UV_array, args.vertex_normal_array);
    populate_vertex_vectors(data1);

    std::string v2_data = get_nth_word(args.triangle_data, ' ', args.idx2);
    PopulateVertexVectorData data2 =
        PopulateVertexVectorData(v2_data, args.vertex_array, args.vertex_UV_array, args.vertex_normal_array);
    populate_vertex_vectors(data2);

    std::string v3_data = get_nth_word(args.triangle_data, ' ', args.idx3);
    PopulateVertexVectorData data3 =
        PopulateVertexVectorData(v3_data, args.vertex_array, args.vertex_UV_array, args.vertex_normal_array);
    populate_vertex_vectors(data3);

    bool loaded_vertices_successfully = data1.v_success && data2.v_success && data3.v_success;

    TriangleCreationResult result;

    if (!loaded_vertices_successfully) {
        return result;
    }

    Triangle* triangle = new Triangle(data1.v, data2.v, data3.v, args.material);

    bool loaded_UV_successfully = data1.uv_success && data2.uv_success && data3.uv_success;
    if (loaded_UV_successfully) {
        triangle->set_vertex_UV(data1.uv, data2.uv, data3.uv);
    }

    bool loaded_normals_successfully = data1.n_success && data2.n_success && data3.n_success;
    if (loaded_normals_successfully && args.enable_smooth_shading) {
        triangle->set_vertex_normals(data1.n, data2.n, data3.n);
    }
    result.success = true;
    result.triangle = triangle;
    return result;
}

int populate_triangle_array(std::string file_name, vec3* vertex_array, vec3* vertex_UV_array, vec3* vertex_normal_array,
                            Object** triangle_array, Material* material, const bool enable_smooth_shading) {
    std::ifstream model_file(file_name);
    std::string line;
    int shape_idx = 0;
    while (std::getline(model_file, line)) {
        std::string first_word = get_nth_word(line, ' ', 0);

        bool is_shape = first_word == "f";
        if (!is_shape) {
            continue;
        }
        int number_of_spaces = number_of_char_occurances(line, ' ');
        bool is_triangle = number_of_spaces == 3;
        bool is_quad = number_of_spaces == 4;
        if (is_triangle) {
            TriangleConstructionArgs args = TriangleConstructionArgs(
                line, 1, 2, 3, material, vertex_array, vertex_UV_array, vertex_normal_array, enable_smooth_shading);
            TriangleCreationResult result = construct_triangle(args);
            if (!result.success) {
                continue;
            }
            triangle_array[shape_idx] = result.triangle;
            shape_idx++;
        }

        else if (is_quad) {
            TriangleConstructionArgs args1 = TriangleConstructionArgs(
                line, 1, 2, 3, material, vertex_array, vertex_UV_array, vertex_normal_array, enable_smooth_shading);
            TriangleCreationResult result1 = construct_triangle(args1);
            if (result1.success) {
                triangle_array[shape_idx] = result1.triangle;
                shape_idx++;
            }

            TriangleConstructionArgs args2 = TriangleConstructionArgs(
                line, 1, 3, 4, material, vertex_array, vertex_UV_array, vertex_normal_array, enable_smooth_shading);
            TriangleCreationResult result2 = construct_triangle(args2);
            if (result1.success) {
                triangle_array[shape_idx] = result2.triangle;
                shape_idx++;
            }
        }
    }
    return shape_idx;
}

ObjectUnion* load_object_model(std::string file_name, Material* material, const bool enable_smooth_shading,
                               const ObjectTransform& transform) {
    DataSizes nums = get_vertex_data_sizes(file_name);

    vec3* vertex_array = new vec3[nums.num_vertices];
    vec3* vertex_UV_array = new vec3[nums.num_vertex_UVs];
    vec3* vertex_normal_array = new vec3[nums.num_vertex_normals];

    populate_vertex_arrays(file_name, vertex_array, vertex_UV_array, vertex_normal_array);
    if (transform.move_object) {
        change_vectors(transform.center, transform.size, vertex_array, nums.num_vertices);
    }

    Object** triangles = new Object*[nums.num_triangles];
    int num_valid_triangles = populate_triangle_array(file_name, vertex_array, vertex_UV_array, vertex_normal_array,
                                                      triangles, material, enable_smooth_shading);
    ObjectUnion* loaded_object = new ObjectUnion(triangles, num_valid_triangles, true);
    return loaded_object;
}
