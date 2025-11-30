#include <string>
#include <thread>
#include "constants.h"

namespace constants {
// Default value of settings. Is overidden by the provided settings.json, if any.
int WIDTH = 1000;
int HEIGHT = 1000;
int samples_per_pixel = 10;
int max_recursion_depth = 100;
int force_tracing_limit = 3;
double air_refractive_index = 1;
int number_of_threads = 1;

bool enable_next_event_estimation = true;
bool enable_anti_aliasing = true;

bool enable_denoising = false;
int denoising_iterations = 5;
double sigma_rt = 0.7;
double sigma_x = 0.5;
double sigma_n = 0.4;

const int max_number_of_threads = std::thread::hardware_concurrency() - 1;
const double EPSILON = 0.000001;
const double max_ray_distance = 1.0 / 0.0;

const std::string raw_file_name = "./temp/raw.dat";
const std::string raw_denoised_file_name = "./temp/raw_denoised.dat";
}