#include <string>
#include <thread>
#include "constants.h"

namespace constants {
// Default value of settings. Is overidden by the provided settings.json, if any.
int WIDTH = 1000;
int HEIGHT = 1000;
int samples_per_pixel = 10;
int samples_per_iteration = 16;
int max_recursion_depth = 100;
int min_recursion_steps = 3;
int number_of_threads = 1;

bool enable_next_event_estimation = true;
bool enable_anti_aliasing = true;

bool enable_atrous_filtering = false;
int denoising_iterations = 5;
double sigma_rt = 0.7;
double sigma_x = 0.5;
double sigma_n = 0.4;

bool enable_median_filtering = false;
int median_kernel_size = 3;

const int max_number_of_threads = static_cast<int>(std::thread::hardware_concurrency()) - 1;
const double EPSILON = 0.000001;
const double max_ray_distance = 1.0 / 0.0;

const std::string raw_file_name = "./temp/raw.dat";
const std::string raw_denoised_file_name = "./temp/raw_denoised.dat";
}