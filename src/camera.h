#ifndef CAMERA_H
#define CAMERA_H
#include "vec3.h"
#include "constants.h"


class Camera{
    public:
        vec3 position;

        Camera(){}
        Camera(vec3 _position=vec3(0,0,0), vec3 _viewing_direction=vec3(0,0,1), vec3 _y_vector=vec3(0,1,0));

    vec3 index_to_position(double x, double y) const;
    vec3 get_starting_directions(double x, double y) const;

    private:
        vec3 viewing_direction;
        vec3 screen_x_vector;
        vec3 screen_y_vector;
        vec3 screen_position;
        double screen_width;
        double screen_height;
};


#endif
