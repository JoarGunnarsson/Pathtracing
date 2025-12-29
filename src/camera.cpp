#include "camera.h"

vec3 rotate(vec3 const& p1, double alpha, double beta, double gamma) {
    // Rotates the given vector using YZX Tait–Bryan angles

    alpha = alpha * M_PI / 180;
    beta = beta * M_PI / 180;
    gamma = gamma * M_PI / 180;
    double x = p1[0];
    double y = p1[1];
    double z = p1[2];

    double cos_beta, sin_beta, cos_alpha, sin_alpha, cos_gamma, sin_gamma;
    double e1, e2, e3;

    cos_alpha = cos(alpha);
    sin_alpha = sin(alpha);

    cos_beta = cos(beta);
    sin_beta = sin(beta);

    cos_gamma = cos(gamma);
    sin_gamma = sin(gamma);

    e1 = x * (cos_alpha * cos_beta) + y * (sin_alpha * sin_gamma - cos_alpha * cos_gamma * sin_beta) +
         z * (cos_gamma * sin_alpha + cos_alpha * sin_beta * sin_gamma);
    e2 = x * (sin_beta) + y * (cos_beta * cos_gamma) + z * (-cos_beta * sin_gamma);
    e3 = x * (-cos_beta * sin_alpha) + y * (cos_alpha * sin_gamma + cos_gamma * sin_alpha * sin_beta) +
         z * (cos_alpha * cos_gamma - sin_alpha * sin_beta * sin_gamma);

    return vec3(e1, e2, e3);
}

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
