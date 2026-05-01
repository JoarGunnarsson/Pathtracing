#include <string>
#include <thread>
#include "constants.h"

namespace constants {
// Default value of settings. Is overidden by the provided settings.json, if any.
size_t WIDTH = 1000;
size_t HEIGHT = 1000;
size_t samples_per_pixel = 10;
size_t samples_per_iteration = 16;
int max_recursion_depth = 100;
int min_recursion_steps = 3;
size_t number_of_threads = 1;
bool use_gamma_correction = true;

bool enable_next_event_estimation = true;
bool enable_anti_aliasing = true;

bool enable_atrous_filtering = false;
int denoising_iterations = 5;
double sigma_rt = 0.7;
double sigma_x = 0.5;
double sigma_n = 0.4;

bool enable_median_filtering = false;
int median_kernel_size = 3;
double median_filter_threshold = 0;

const size_t max_number_of_threads = std::max(std::thread::hardware_concurrency() - 1, static_cast<unsigned int>(1));
const double EPSILON = 0.000001;
const double max_ray_distance = 1.0 / 0.0;

const std::string raw_pixeldata_file = "./temp/raw_pixel.dat";
const std::string raw_positiondata_file = "./temp/raw_position.dat";
const std::string raw_normaldata_file = "./temp/raw_normal.dat";
const std::string raw_denoised_file_name = "./temp/raw_denoised.dat";
}