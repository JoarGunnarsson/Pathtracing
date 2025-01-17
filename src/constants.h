#ifndef CONSTANTS_H
#define CONSTANTS_H


namespace constants
{
    const int WIDTH = 100;
    const int HEIGHT = 100;
    const int samples_per_pixel = 100;
    const int max_recursion_depth = 100;
    const int force_tracing_limit = 3;

    const double EPSILON = 0.000001;
    const double max_ray_distance = 1.0 / 0.0;
    const double air_refractive_index = 1;

    const bool enable_next_event_estimation = true;

    const bool enable_anti_aliasing = true;

    const bool enable_denoising = false;
    const int denoising_iterations = 5;
    const double sigma_rt = 2;
    const double sigma_x = 5;
    const double sigma_n = 0.1;

    const std::string raw_output_file_name = "./temp/raw_data.txt";
    const std::string denoised_output_file_name = "./temp/denoised_data.txt";
}

#endif