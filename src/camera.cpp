#include "camera.h"
#include "utils.h"

Camera::Camera(vec3 const& _position, double X, double Y, double Z) {
    position = _position;
    vec3 forward = vec3(0, -1, 0);
    vec3 up = vec3(0, 0, -1);
    viewing_direction = rotate(forward, Y, Z, X);
    screen_y_vector = rotate(up, Y, Z, X);

    screen_width = 1.0;
    screen_height = screen_width * (double) constants::HEIGHT / (double) constants::WIDTH;
    screen_x_vector = cross_vectors(viewing_direction, screen_y_vector);
    screen_position = position + viewing_direction;
}

vec3 Camera::index_to_position(double x, double y) const {
    double local_x_coordinate = x * screen_width / (double) constants::WIDTH - (double) screen_width / 2.0;
    vec3 local_x = screen_x_vector * local_x_coordinate;

    double local_y_coordinate = y * screen_height / (double) constants::HEIGHT - (double) screen_height / 2.0;
    vec3 local_y = screen_y_vector * local_y_coordinate;

    return local_x + local_y + screen_position;
}

vec3 Camera::get_starting_directions(double x, double y) const {
    vec3 pixel_vector = index_to_position(x, y);
    vec3 direction_vector = pixel_vector - position;
    return normalize_vector(direction_vector);
}
