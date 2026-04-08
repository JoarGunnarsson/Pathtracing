#ifndef CAMERA_H
#define CAMERA_H
#include "vec3.h"
#include "constants.h"
#include "utils.h"

class Camera {
  public:
    vec3 position;

    Camera() {}
    Camera(vec3 const& _position, const double X, const double Y, const double Z, const double _camera_width,
           const double _focal_length, const std::string& _depth_of_field_mode, const double _aperture_size,
           const double _focus_distance);

    vec3 index_to_position(const double x, const double y) const;
    void adjust_depth_of_field(Ray& ray) const;
    Ray make_ray(const double x, const double y) const;

  private:
    vec3 viewing_direction;
    vec3 screen_x_vector;
    vec3 screen_y_vector;
    vec3 screen_position;
    double camera_width;
    double screen_height;
    double focal_length;
    std::string depth_of_field_mode;
    double aperture_size;
    double focus_distance;
};

#endif
