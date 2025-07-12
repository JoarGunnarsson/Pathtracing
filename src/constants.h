#ifndef CONSTANTS_H
#define CONSTANTS_H


namespace constants
{
    const int WIDTH = 1000;
    const int HEIGHT = 1000;
    const int samples_per_pixel = 10;
    const int max_recursion_depth = 100;
    const int force_tracing_limit = 3;

    const double EPSILON = 0.000001;
    const double max_ray_distance = 1.0 / 0.0;
    const double air_refractive_index = 1;

    const bool enable_next_event_estimation = true;

    const bool enable_anti_aliasing = true;

    const bool enable_denoising = true;
    const int denoising_iterations = 5;
    const double sigma_rt = 1;
    const double sigma_x = 0.5;
    const double sigma_n = 0.4;

    const char* const raw_file_name = "./temp/raw.dat";
    const char* const raw_denoised_file_name = "./temp/raw_denoised.dat";
}

#endif