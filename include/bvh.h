#ifndef BVH_H
#define BVH_H

#include "objects.h"
#include <chrono>

namespace BVH {
vec3 get_max_point(Object** triangles, int number_of_triangles);
vec3 get_min_point(Object** triangles, int number_of_triangles);

struct Interval {
    double min;
    double max;

    Interval() {}
    Interval(const double min, const double max) : min(min), max(max) {}
};

class BoundingBox {
  public:
    double axis_length[3];

    BoundingBox() {}
    BoundingBox(Object** _triangles, int number_of_triangles);

    inline bool is_within_bounds(const double x, const double lower, const double higher) const {
        return lower <= x && x <= higher;
    }

    Interval get_interval(const int axis) const;
    bool intersect(Ray& ray, double& distance) const;

  private:
    vec3 p1;
    vec3 p2;
    Interval x_interval;
    Interval y_interval;
    Interval z_interval;
    double width;
    double height;
    double length;
};

void sort_by_axis(Object** triangles, int number_of_triangles, int axis);

class Node {
  public:
    BoundingBox bounding_box;

    Node() {}
    Node(Object** _triangles, int _number_of_triangles, int _leaf_size = 12);

    int get_split_axis();
    bool intersect(Ray& ray, Hit& hit);

  private:
    int leaf_size;
    bool is_leaf_node;
    Node* node1;
    Node* node2;
    Object** triangles;
    int number_of_triangles;
};

class BoundingVolumeHierarchy {
  public:
    BoundingVolumeHierarchy() {}
    BoundingVolumeHierarchy(Object** triangles, int number_of_triangles, int leaf_size);

    bool intersect(Hit& hit, Ray& ray) const;

  private:
    Node* root_node;
};
}

#endif
