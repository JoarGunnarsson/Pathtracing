#ifndef CAMERA_H
#define CAMERA_H
#include "vec3.h"
#include "constants.h"

class Camera {
  public:
    vec3 position;

    Camera() {}
    Camera(vec3 const& _position, const double X, const double Y, const double Z, const double _camera_width,
           const double _focal_length);

    vec3 index_to_position(const double x, const double y) const;
    vec3 get_starting_directions(const double x, const double y) const;

  private:
    vec3 viewing_direction;
    vec3 screen_x_vector;
    vec3 screen_y_vector;
    vec3 screen_position;
    double camera_width;
    double screen_height;
    double focal_length;
};

#endif
