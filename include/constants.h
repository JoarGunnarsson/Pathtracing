#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <string>

namespace constants {
// Configurable constants
extern int WIDTH;
extern int HEIGHT;
extern int samples_per_pixel;
extern int samples_per_iteration;
extern int max_recursion_depth;
extern int min_recursion_steps;
extern int number_of_threads;

extern bool enable_next_event_estimation;
extern bool enable_anti_aliasing;

extern bool enable_atrous_filtering;
extern int denoising_iterations;
extern double sigma_rt;
extern double sigma_x;
extern double sigma_n;
extern bool enable_median_filtering;
extern int median_kernel_size;

// Fixed constants
extern const int max_number_of_threads;
extern const double EPSILON;
extern const double max_ray_distance;

extern const std::string raw_file_name;
extern const std::string raw_denoised_file_name;
}

#endif
