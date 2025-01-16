#ifndef BVH_H
#define BVH_H

#include "objects.h"


namespace BVH{
    vec3 get_max_point(Object** triangles, int number_of_triangles);
    vec3 get_min_point(Object** triangles, int number_of_triangles);


    class BoundingBox{
        public:
            vec3 p1;
            vec3 p2;
            double width;
            double height;
            double length;
            double axis_length[3];

            BoundingBox(){}
            BoundingBox(Object** _triangles, int number_of_triangles);
        
        inline bool is_within_bounds(const double x, const double lower, const double higher){
            return lower <= x && x <= higher;
        }

        double intersect(const Ray& ray);    
    };


    void sort_by_axis(Object** triangles, int number_of_triangles, int axis);


    class Node{
        public:
            int leaf_size;
            bool is_leaf_node;
            Node* node1;
            Node* node2;
            int depth;
            Object** triangles;
            int number_of_triangles;
            BoundingBox bounding_box;

            Node(){}
            Node(Object** _triangles, int _number_of_triangles, int _leaf_size=12, int depth=0);

        int get_split_axis();
        void intersect(const Ray& ray, Hit& hit);
    };


    class BoundingVolumeHierarchy{
        public:
            Node* root_node;
            
            BoundingVolumeHierarchy(){}
            BoundingVolumeHierarchy(Object** triangles, int number_of_triangles, int leaf_size);

            Hit intersect(const Ray& ray) const;
    };
}

#endif