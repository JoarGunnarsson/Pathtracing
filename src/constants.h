#ifndef CONSTANTS_H
#define CONSTANTS_H


namespace constants
{
    const int WIDTH = 300;
    const int HEIGHT = 300;
    const int samples_per_pixel = 5;
    const int max_recursion_depth = 100;
    const int force_tracing_limit = 3;

    const double EPSILON = 0.000001;
    const double air_refractive_index = 1;

    const bool enable_next_event_estimation = true;

    const bool enable_anti_aliasing = true;

    const bool enable_denoising = true;
    const int denoising_iterations = 5;

    const std::string raw_output_file_name = "./temp/raw_data.txt";
    const std::string denoised_output_file_name = "./temp/denoised_data.txt";
}

#endif