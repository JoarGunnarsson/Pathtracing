#include "bvh.h"


namespace BVH{
    vec3 get_max_point(Object** triangles, int number_of_triangles){
        if (number_of_triangles == 0){
            return vec3(0,0,0);
        }
        vec3 max_point = triangles[0] -> max_axis_point();
        for (int i = 1; i < number_of_triangles; i++){
            vec3 this_max_point = triangles[i] -> max_axis_point();
            for (int j = 0; j < 3; j++){
                max_point.e[j] = std::max(this_max_point[j], max_point[j]);
            }
        }
        return max_point;
    }


    vec3 get_min_point(Object** triangles, int number_of_triangles){
        if (number_of_triangles == 0){
            return vec3(0,0,0);
        }
        vec3 min_point = triangles[0] -> min_axis_point();
        for (int i = 1; i < number_of_triangles; i++){
            vec3 this_min_point = triangles[i] -> min_axis_point();
            for (int j = 0; j < 3; j++){
                min_point.e[j] = std::min(this_min_point[j], min_point[j]);
            }
        }
        return min_point;
    }


    BoundingBox::BoundingBox(Object** _triangles, int number_of_triangles){
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

    Interval BoundingBox::get_interval(const int axis) const{
        if (axis == 0){
            return x_interval;
        }
        else if (axis == 1){
            return y_interval;
        }
        return z_interval;
    }

    bool BoundingBox::intersect(Ray& ray, double& distance) const{
        Interval ray_interval(-1, ray.t_max);

        for (int axis = 0; axis < 3; axis++){
            Interval ax = get_interval(axis);
            double d_inv = 1.0 / ray.direction_vector[axis];

            double t0 = (ax.min - ray.starting_position[axis]) * d_inv;
            double t1 = (ax.max - ray.starting_position[axis]) * d_inv;

            if (t0 < t1){
                if (t0 > ray_interval.min){
                    ray_interval.min = t0;
                }

                if (t1 < ray_interval.max){
                    ray_interval.max = t1;
                }
            }
            else{
                if (t1 > ray_interval.min){
                    ray_interval.min = t1;
                }

                if (t0 < ray_interval.max){
                    ray_interval.max = t0;
                }
            }

            if (ray_interval.max <= ray_interval.min){
                return false;
            }
        }
        if (ray_interval.min > constants::EPSILON){
            distance = ray_interval.min;
            return true;
        }
        else if (ray_interval.max > constants::EPSILON){
            distance = ray_interval.max;
            return true;
        }
        return false;
    }    


    void sort_by_axis(Object** triangles, int number_of_triangles, int axis){
        std::sort(triangles, triangles + number_of_triangles, [axis](Object* obj1, Object* obj2){ 
            return (obj1 -> compute_centroid())[axis] < (obj2 -> compute_centroid())[axis]; 
            });
    }



    Node::Node(Object** _triangles, int _number_of_triangles, int _leaf_size, int depth){
        leaf_size = _leaf_size;
        bounding_box = BoundingBox(_triangles, _number_of_triangles);
        if (_number_of_triangles <= leaf_size){
            triangles = _triangles;
            number_of_triangles = _number_of_triangles;
            is_leaf_node = true;
            return;
        }
        is_leaf_node = false;
        int axis = get_split_axis();

        sort_by_axis(_triangles, _number_of_triangles, axis);
        int split_index = _number_of_triangles / 2;

        Object** node1_triangles = new Object*[split_index];
        Object** node2_triangles = new Object*[_number_of_triangles - split_index];

        for (int i = 0; i < split_index; i++){
            node1_triangles[i] = _triangles[i];
        }
        for (int i = split_index; i < _number_of_triangles; i++){
            node2_triangles[i - split_index] = _triangles[i];
        }

        node1 = new Node(node1_triangles, split_index, _leaf_size, depth+1);
        node2 = new Node(node2_triangles, _number_of_triangles - split_index, _leaf_size, depth+1);
    }

    int Node::get_split_axis(){
        int axis;
        double max_length = 0;
        for (int i = 0; i < 3; i++){
            if (bounding_box.axis_length[i] >= max_length){
                axis = i;
                max_length = bounding_box.axis_length[i];
            }
        }
        return axis;
    }

    bool Node::intersect(Ray& ray, Hit& hit){
        if (is_leaf_node){
            if (number_of_triangles == 0){
                return false;
            }
            Hit triangle_hit;
            bool hits_triangles = find_closest_hit(triangle_hit, ray, triangles, number_of_triangles);
            if (hits_triangles && triangle_hit.distance < hit.distance){
                hit.distance = triangle_hit.distance;
                hit.primitive_ID = triangle_hit.primitive_ID;
                return true;
            }
            return false;
        }

        double d1;
        bool bvh1_hit = node1 -> bounding_box.intersect(ray, d1);
        double d2;
        bool bvh2_hit = node2 -> bounding_box.intersect(ray, d2);
        
        if (bvh1_hit && bvh2_hit){
            if (d1 > d2){
                bool node1_success = node1 -> intersect(ray, hit);
                bool node2_success;
                if (d2 < hit.distance || hit.distance == -1){
                    node2_success = node2 -> intersect(ray, hit);
                }
                return node1_success || node2_success;
            }
            else{
                bool node1_success;
                bool node2_success = node2 -> intersect(ray, hit);
                if (d1 < hit.distance || hit.distance == -1){
                    node1_success = node1 -> intersect(ray, hit);
                }
                return node1_success || node2_success;
            }

        }
        else if (bvh1_hit){
            return node1 -> intersect(ray, hit);
        }
        
        else if (bvh2_hit){
            return node2 -> intersect(ray, hit);
        }
        return false;
    }


    BoundingVolumeHierarchy::BoundingVolumeHierarchy(Object** triangles, int number_of_triangles, int leaf_size){
        root_node = new Node(triangles, number_of_triangles, leaf_size);
    }

    bool BoundingVolumeHierarchy::intersect(Hit& hit, Ray& ray) const{
        double distance_to_bounding_box;
        
        if (!root_node -> bounding_box.intersect(ray, distance_to_bounding_box)){
            return false;
        }

        //hit.distance = -1;
        //hit.primitive_ID = -1;
        return root_node -> intersect(ray, hit);
    }
}
