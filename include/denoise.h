
#ifndef DENOISE_H
#define DENOISE_H

#include "vec3.h"
#include "utils.h"
#include "constants.h"
#include "scene.h"

struct KernelData {
    double kernel[25] = {1.0 / 16.0, 1.0 / 16.0, 1.0 / 16.0, 1.0 / 16.0, 1.0 / 16.0, 1.0 / 16.0, 1.0 / 4.0,
                         1.0 / 4.0,  1.0 / 4.0,  1.0 / 16.0, 1.0 / 16.0, 1.0 / 4.0,  3.0 / 8.0,  1.0 / 4.0,
                         1.0 / 16.0, 1.0 / 16.0, 1.0 / 4.0,  1.0 / 4.0,  1.0 / 4.0,  1.0 / 16.0, 1.0 / 16.0,
                         1.0 / 16.0, 1.0 / 16.0, 1.0 / 16.0, 1.0 / 16.0};

    double sigma_rt;
    double sigma_x;
    double sigma_n;

    int hole_width = 0;
};

void get_image_coordinates(int& x, int& y, const size_t idx);
size_t idx_from_coordinates(const size_t x, const size_t y, const size_t width);
bool is_out_of_bounds(const int x, const int y);
int clamp_x_coordinate(const int x);
int clamp_y_coordinate(const int y);

double compute_weight_component(const vec3& vec_p, const vec3& vec_q, const double standard_deviation);
double compute_weight(const size_t p, const size_t q, const KernelData& kernel_data, const PixelBuffers& buffers);
int expand_kernel_idx(const int idx, const int hole_width);

vec3 blur_pixel(const size_t p, const KernelData& kernel_data, PixelBuffers& buffers);
void one_denoising_iteration(const KernelData& kernel_data, PixelBuffers& buffers);
void median_filter(PixelBuffers& buffers);
void denoise(PixelBuffers& buffers, std::vector<DenoisingTask> pipeline);

#endif
