#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <string>

namespace constants {
extern int WIDTH;
extern int HEIGHT;
extern int samples_per_pixel;
extern int max_recursion_depth;
extern int force_tracing_limit;
extern double air_refractive_index;
extern int number_of_threads;

extern bool enable_next_event_estimation;
extern bool enable_anti_aliasing;

extern bool enable_denoising;
extern int denoising_iterations;
extern double sigma_rt;
extern double sigma_x;
extern double sigma_n;

extern const int max_number_of_threads;
extern const double EPSILON;
extern const double max_ray_distance;

extern const std::string raw_file_name;
extern const std::string raw_denoised_file_name;
}

#endif
