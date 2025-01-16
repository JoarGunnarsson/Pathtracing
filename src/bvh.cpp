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

        width = p2[0] - p1[0];
        length = p2[1] - p1[1];
        height = p2[2] - p1[2];
        axis_length[0] = width;
        axis_length[1] = length;
        axis_length[2] = height;
    }

    double BoundingBox::intersect(const Ray& ray){
        double t[6];
        bool inside_bounds[6];
        for (int i = 0; i < 3; i++){
            if (std::abs(ray.direction_vector[i]) < constants::EPSILON){
                t[i] = -1;
                inside_bounds[i] = false;
                continue;
            }
            t[i] = (p1[i] - ray.starting_position[i]) / ray.direction_vector[i];
            vec3 hit_point = ray.direction_vector * t[i] + ray.starting_position;
            vec3 difference_vector = hit_point - p1;
            if (i == 0){
                inside_bounds[i] = is_within_bounds(difference_vector[1], 0, length) && is_within_bounds(difference_vector[2], 0, height);
            }
            else if (i ==  1){
                inside_bounds[i] = is_within_bounds(difference_vector[0], 0, width) && is_within_bounds(difference_vector[2], 0, height);
            }
            else if (i == 2){
                inside_bounds[i] = is_within_bounds(difference_vector[0], 0, width) && is_within_bounds(difference_vector[1], 0, length);
            }
        }

        for (int i = 0; i < 3; i++){
            if (std::abs(ray.direction_vector[i]) < constants::EPSILON){
                t[i+3] = -1;
                inside_bounds[i+3] = false;
                continue;
            }
            t[i+3] = (p2[i] - ray.starting_position[i]) / ray.direction_vector[i];
            vec3 hit_point = ray.direction_vector * t[i+3] + ray.starting_position;
            vec3 difference_vector = hit_point - p2;
            if (i == 0){
                inside_bounds[i+3] = is_within_bounds(difference_vector[1], -length, 0) && is_within_bounds(difference_vector[2], -height, 0);
            }
            else if (i ==  1){
                inside_bounds[i+3] = is_within_bounds(difference_vector[0], -width, 0) && is_within_bounds(difference_vector[2], -height, 0);
            }
            else if (i == 2){
                inside_bounds[i+3] = is_within_bounds(difference_vector[0], -width, 0) && is_within_bounds(difference_vector[1], -length, 0);
            }
        }

        double min_t = -1;
        for (int i = 0; i < 6; i++){
            if (inside_bounds[i] && (min_t == -1 || min_t > t[i]) && t[i] > constants::EPSILON){
                min_t = t[i];
            }
        }

        return min_t;
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

    void Node::intersect(const Ray& ray, Hit& hit){
        if (is_leaf_node){
            if (number_of_triangles == 0){
                return;
            }
            Hit closest_hit = find_closest_hit(ray, triangles, number_of_triangles);
            if (closest_hit.distance > constants::EPSILON && (closest_hit.distance < hit.distance || hit.distance == -1)){
                hit.distance = closest_hit.distance;
                hit.object_ID = closest_hit.object_ID;
            }
            return;
        }

        double d1 = node1 -> bounding_box.intersect(ray);
        double d2 = node2 -> bounding_box.intersect(ray);
        
        bool node1_hit = d1 > constants::EPSILON && (d1 < hit.distance || hit.distance == -1);
        bool node2_hit = d2 > constants::EPSILON && (d2 < hit.distance || hit.distance == -1);
        
        if (node1_hit && node2_hit){
            if (d1 > d2){
                node1 -> intersect(ray, hit);
                if (d2 < hit.distance || hit.distance == -1){
                    node2 -> intersect(ray, hit);
                }
            }
            else{
                node2 -> intersect(ray, hit);
                if (d1 < hit.distance || hit.distance == -1){
                    node1 -> intersect(ray, hit);
                }
            }

        }
        else if (node1_hit){
            node1 -> intersect(ray, hit);
        }
        else if (node2_hit){
            node2 -> intersect(ray, hit);
        }
    }


    BoundingVolumeHierarchy::BoundingVolumeHierarchy(Object** triangles, int number_of_triangles, int leaf_size){
        root_node = new Node(triangles, number_of_triangles, leaf_size);
    }

    Hit BoundingVolumeHierarchy::intersect(const Ray& ray) const{
        double distance_to_bounding_box = root_node -> bounding_box.intersect(ray);

        Hit hit;
        hit.distance = -1;
        hit.object_ID = -1;
        if (distance_to_bounding_box > constants::EPSILON){
            root_node -> intersect(ray, hit);
        }
        return hit;
    }
}
