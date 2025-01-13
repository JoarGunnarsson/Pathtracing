#ifndef CONSTANTS_H
#define CONSTANTS_H


namespace constants
{
    const double EPSILON = 0.000001;
    const int samples_per_pixel = 10;
    const int WIDTH = 1000;
    const int HEIGHT = 1000;
    const bool enable_next_event_estimation = true;
    const double air_refractive_index = 1;
    const int max_recursion_depth = 100;
    const bool use_denoising = true;
    const int denoising_iterations = 5;
}

#endif