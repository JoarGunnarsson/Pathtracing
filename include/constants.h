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
extern int bvh_leaf_size;

extern bool enable_next_event_estimation;
extern bool enable_anti_aliasing;

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
