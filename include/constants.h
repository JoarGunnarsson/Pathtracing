#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <string>

namespace constants {
// Configurable constants
extern size_t WIDTH;
extern size_t HEIGHT;
extern size_t samples_per_pixel;
extern size_t samples_per_iteration;
extern int max_recursion_depth;
extern int min_recursion_steps;
extern size_t number_of_threads;
extern bool use_gamma_correction;

extern bool enable_next_event_estimation;
extern bool enable_anti_aliasing;

extern bool enable_atrous_filtering;
extern int denoising_iterations;
extern double sigma_rt;
extern double sigma_x;
extern double sigma_n;

extern bool enable_median_filtering;
extern int median_kernel_size;
extern double median_filter_threshold;

// Fixed constants
extern const size_t max_number_of_threads;
extern const double EPSILON;
extern const double max_ray_distance;

extern const std::string raw_pixeldata_file;
extern const std::string raw_positiondata_file;
extern const std::string raw_normaldata_file;
extern const std::string raw_denoised_file_name;
}

#endif
