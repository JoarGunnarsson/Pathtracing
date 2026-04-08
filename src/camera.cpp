#include "camera.h"
#include "utils.h"

Camera::Camera(vec3 const& _position, double X, double Y, double Z, const double _camera_width,
               const double _focal_length, const std::string& _depth_of_field_mode, const double _aperture_size,
               const double _focus_distance) {
    position = _position;
    vec3 forward = vec3(0, -1, 0);
    vec3 up = vec3(0, 0, -1);
    viewing_direction = rotate(forward, Y, Z, X);
    screen_y_vector = rotate(up, Y, Z, X);

    camera_width = _camera_width;
    focal_length = _focal_length;
    depth_of_field_mode = _depth_of_field_mode;
    aperture_size = _aperture_size;
    focus_distance = _focus_distance;

    screen_height = camera_width * (double) constants::HEIGHT / (double) constants::WIDTH;
    screen_x_vector = cross_vectors(viewing_direction, screen_y_vector);
    screen_position = position + viewing_direction;
}

vec3 Camera::index_to_position(double x, double y) const {
    double local_x_coordinate = x * camera_width / (double) constants::WIDTH - (double) camera_width / 2.0;
    vec3 local_x = screen_x_vector * local_x_coordinate;

    double local_y_coordinate = y * screen_height / (double) constants::HEIGHT - (double) screen_height / 2.0;
    vec3 local_y = screen_y_vector * local_y_coordinate;

    return position + local_x + local_y + viewing_direction * focal_length;
}

void Camera::adjust_depth_of_field(Ray& ray) const {
    double t = focus_distance / dot_vectors(ray.direction_vector, viewing_direction);
    vec3 intersection_with_focal_plane = ray.starting_position + t * ray.direction_vector;

    vec3 adjust;
    if (depth_of_field_mode == "circle") {
        double r = aperture_size / 2.0 * std::sqrt(random_uniform(0, 1));
        double phi = random_uniform(0, 2 * M_PI);
        adjust = r * (std::cos(phi) * screen_x_vector + std::sin(phi) * screen_y_vector);
    }
    else if (depth_of_field_mode == "square") {
        double r1 = random_uniform(0, aperture_size) - aperture_size / 2.0;
        double r2 = random_uniform(0, aperture_size) - aperture_size / 2.0;
        adjust = (r1 * screen_y_vector + r2 * screen_x_vector);
    }
    else {
        throw std::runtime_error("Invalid 'depth_of_field_mode' " + depth_of_field_mode + " in scene definition.");
    }

    ray.starting_position += adjust;
    ray.direction_vector = normalize_vector(intersection_with_focal_plane - ray.starting_position);
}

Ray Camera::make_ray(const double x, const double y) const {
    vec3 pixel_vector = index_to_position(x, y); // Position on image plane.
    vec3 direction_vector = normalize_vector(pixel_vector - position);

    Ray ray;
    ray.starting_position = position;
    ray.direction_vector = direction_vector;
    if (depth_of_field_mode != "none") {
        adjust_depth_of_field(ray);
    }

    ray.type = TRANSMITTED;
    return ray;
}
