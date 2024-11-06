#ifndef CAMERA_H
#define CAMERA_H
#include "vec3.h"
#include "constants.h"


class Camera{
    public:
        vec3 position;
        vec3 viewing_direction;
        vec3 screen_y_vector;
        vec3 screen_position;
        vec3 screen_x_vector;
        double screen_width;
        double screen_height;
        Camera(){}
        Camera(vec3 _position=vec3(0,0,0), vec3 _viewing_direction=vec3(0,0,1), vec3 _y_vector=vec3(0,1,0)){
            position = _position;
            viewing_direction = normalize_vector(_viewing_direction);
            if (dot_vectors(viewing_direction, _y_vector) != 0){
                vec3 perpendicularVector = cross_vectors(viewing_direction, _y_vector);
                _y_vector = cross_vectors(perpendicularVector, viewing_direction);
            }
            screen_y_vector = normalize_vector(_y_vector);
            screen_width = 1.0;
            screen_height = screen_width * (double) constants::HEIGHT / (double) constants::WIDTH; 
            screen_x_vector = cross_vectors(viewing_direction, screen_y_vector);
            screen_position = position + viewing_direction;
        }

    vec3 index_to_position(double x, double y){
        double local_x_coordinate = x * screen_width / (double) constants::WIDTH - (double) screen_width / 2.0;
        vec3 localX = screen_x_vector * local_x_coordinate;

        double local_y_coordinate = y * screen_height / (double) constants::HEIGHT - (double) screen_height / 2.0;
        vec3 localY = screen_y_vector * local_y_coordinate;

        return localX + localY + screen_position;
    }

    vec3 get_starting_directions(double x, double y){
        vec3 pixel_vector = index_to_position(x, y);
        vec3  direction_vector = pixel_vector - position;
        return normalize_vector(direction_vector);
 }
}; 


#endif