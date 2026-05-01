#include "bvh.h"
#include <algorithm>

namespace BVH {
vec3 get_max_point(Object const* const* triangles, const int number_of_triangles) {
    if (number_of_triangles == 0) {
        return vec3(0, 0, 0);
    }
    vec3 max_point = triangles[0]->max_axis_point();
    for (int i = 1; i < number_of_triangles; i++) {
        vec3 this_max_point = triangles[i]->max_axis_point();
        for (int j = 0; j < 3; j++) {
            max_point.e[j] = std::max(this_max_point[j], max_point[j]);
        }
    }
    return max_point + vec3(constants::EPSILON);
}

vec3 get_min_point(Object const* const* triangles, const int number_of_triangles) {
    if (number_of_triangles == 0) {
        return vec3(0, 0, 0);
    }
    vec3 min_point = triangles[0]->min_axis_point();
    for (int i = 1; i < number_of_triangles; i++) {
        vec3 this_min_point = triangles[i]->min_axis_point();
        for (int j = 0; j < 3; j++) {
            min_point.e[j] = std::min(this_min_point[j], min_point[j]);
        }
    }
    return min_point - vec3(constants::EPSILON);
}

BoundingBox::BoundingBox(Object const* const* _triangles, const int number_of_triangles) {
    p1 = get_min_point(_triangles, number_of_triangles);
    p2 = get_max_point(_triangles, number_of_triangles);
    x_interval = Interval(p1[0], p2[0]);
    y_interval = Interval(p1[1], p2[1]);
    z_interval = Interval(p1[2], p2[2]);

    width = p2[0] - p1[0];
    length = p2[1] - p1[1];
    height = p2[2] - p1[2];

    axis_length[0] = width;
    axis_length[1] = length;
    axis_length[2] = height;
}

Interval BoundingBox::get_interval(const int axis) const {
    switch (axis) {
        case 0:
            return x_interval;
        case 1:
            return y_interval;
        default:
            return z_interval;
    }
}

bool BoundingBox::intersect(Ray& ray, double& distance) const {
    Interval ray_interval(0, constants::max_ray_distance);

    for (int axis = 0; axis < 3; axis++) {
        Interval ax = get_interval(axis);
        const double d_inv = 1.0 / ray.direction_vector[axis];

        const double t0 = (ax.min - ray.starting_position[axis]) * d_inv;
        const double t1 = (ax.max - ray.starting_position[axis]) * d_inv;

        const double t_min = std::min(t0, t1);
        const double t_max = std::max(t0, t1);

        ray_interval.min = std::max(ray_interval.min, t_min);
        ray_interval.max = std::min(ray_interval.max, t_max);

        if (ray_interval.max <= ray_interval.min) {
            return false;
        }
    }

    distance = std::max(ray_interval.min, constants::EPSILON);
    return true;
}

void sort_by_axis(Object** triangles, const int number_of_triangles, const int axis) {
    std::sort(triangles, triangles + number_of_triangles, [axis](Object const* obj1, Object const* obj2) {
        return (obj1->compute_centroid())[axis] < (obj2->compute_centroid())[axis];
    });
}

int surface_area_heuristic(Object** triangles, const int number_of_triangles) {
    // 'triangles' must be sorted with respect to some axis beforehand.

    if (number_of_triangles <= 1) {
        throw std::runtime_error("Cannot split triangle array of size 1 or less!");
    }
    int split_index = 0;
    int bucket_width =
        ceil_division(number_of_triangles - (constants::bvh_n_axis_splits - 1), constants::bvh_n_axis_splits);

    double LSA, RSA;
    double best_metric = -1;
    double metric;
    for (int i = bucket_width; i < number_of_triangles - bucket_width; i += bucket_width + 1) {
        vec3 left_max_point = get_max_point(triangles, i + 1);
        vec3 left_min_point = get_min_point(triangles, i + 1);
        double lw = left_max_point[0] - left_min_point[0];
        double ld = left_max_point[1] - left_min_point[1];
        double lh = left_max_point[2] - left_min_point[2];

        vec3 right_max_point = get_max_point(triangles + i, number_of_triangles - (i + 1));
        vec3 right_min_point = get_min_point(triangles + i, number_of_triangles - (i + 1));
        double rw = right_max_point[0] - right_min_point[0];
        double rd = right_max_point[1] - right_min_point[1];
        double rh = right_max_point[2] - right_min_point[2];

        LSA = 2 * lw * ld + 2 * lw * lh + 2 * ld * lh;
        RSA = 2 * rw * rd + 2 * rw * rh + 2 * rd * rh;
        metric = LSA * (i + 1) + RSA * (number_of_triangles - (i + 1));

        if (metric < best_metric || best_metric == -1) {
            best_metric = metric;
            split_index = i;
        }
    }
    return split_index;
}

Node::Node(Object** _triangles, const int _number_of_triangles, const int _leaf_size) {
    leaf_size = _leaf_size;
    bounding_box = BoundingBox(_triangles, _number_of_triangles);

    if (_number_of_triangles <= leaf_size) {
        triangles = _triangles;
        number_of_triangles = _number_of_triangles;
        is_leaf_node = true;
        return;
    }

    is_leaf_node = false;
    int axis = get_split_axis();

    sort_by_axis(_triangles, _number_of_triangles, axis);

    int split_index = surface_area_heuristic(_triangles, _number_of_triangles);

    size_t left_split_size = static_cast<size_t>(split_index + 1);
    size_t right_split_size = static_cast<size_t>(_number_of_triangles - split_index - 1);

    Object** node1_triangles = new Object*[left_split_size];
    Object** node2_triangles = new Object*[right_split_size];

    // Triangle with index 'split_index' will be included in the left node.
    for (size_t i = 0; i < left_split_size; i++) {
        node1_triangles[i] = _triangles[i];
    }
    for (size_t i = 0; i < right_split_size; i++) {
        node2_triangles[i] = _triangles[i + left_split_size];
    }

    // TODO: This should be a size_t, but it's too much work fixing it right now.
    // Might break for really large 3D models.
    node1 = new Node(node1_triangles, static_cast<int>(left_split_size), _leaf_size);
    node2 = new Node(node2_triangles, static_cast<int>(right_split_size), _leaf_size);
}

int Node::get_split_axis() {
    int axis;
    double max_length = 0;
    for (int i = 0; i < 3; i++) {
        if (bounding_box.axis_length[i] >= max_length) {
            axis = i;
            max_length = bounding_box.axis_length[i];
        }
    }
    return axis;
}

bool Node::intersect(Ray& ray, Hit& hit) {
    if (is_leaf_node) {
        if (number_of_triangles == 0) {
            return false;
        }
        Hit triangle_hit;
        bool hits_triangles = find_closest_hit(triangle_hit, ray, triangles, number_of_triangles);
        if (hits_triangles && triangle_hit.distance < hit.distance && triangle_hit.distance > constants::EPSILON) {
            hit = triangle_hit;
            return true;
        }
        return false;
    }

    double d1, d2;
    bool bvh1_hit = node1->bounding_box.intersect(ray, d1);
    bool bvh2_hit = node2->bounding_box.intersect(ray, d2);

    bool node1_success = false;
    bool node2_success = false;
    if (bvh1_hit && bvh2_hit) {
        if (d1 < d2) {
            node1_success = node1->intersect(ray, hit);
            if (d2 < hit.distance || hit.distance == -1) {
                node2_success = node2->intersect(ray, hit);
            }
        }
        else {
            node2_success = node2->intersect(ray, hit);
            if (d1 < hit.distance || hit.distance == -1) {
                node1_success = node1->intersect(ray, hit);
            }
        }
        return node1_success || node2_success;
    }
    else if (bvh1_hit) {
        return node1->intersect(ray, hit);
    }

    else if (bvh2_hit) {
        return node2->intersect(ray, hit);
    }
    return false;
}

BoundingVolumeHierarchy::BoundingVolumeHierarchy(Object** triangles, const int number_of_triangles,
                                                 const int leaf_size) {
    root_node = new Node(triangles, number_of_triangles, leaf_size);
}

bool BoundingVolumeHierarchy::intersect(Hit& hit, Ray& ray) const {
    double distance_to_bounding_box;
    if (!root_node->bounding_box.intersect(ray, distance_to_bounding_box)) {
        return false;
    }

    return root_node->intersect(ray, hit);
    ;
}
}
