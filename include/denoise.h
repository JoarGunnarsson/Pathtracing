
#ifndef DENOISE_H
#define DENOISE_H

#include "vec3.h"
#include "utils.h"
#include "constants.h"

struct KernelData {
    double kernel[25] = {1.0 / 16.0, 1.0 / 16.0, 1.0 / 16.0, 1.0 / 16.0, 1.0 / 16.0, 1.0 / 16.0, 1.0 / 4.0,
                         1.0 / 4.0,  1.0 / 4.0,  1.0 / 16.0, 1.0 / 16.0, 1.0 / 4.0,  3.0 / 8.0,  1.0 / 4.0,
                         1.0 / 16.0, 1.0 / 16.0, 1.0 / 4.0,  1.0 / 4.0,  1.0 / 4.0,  1.0 / 16.0, 1.0 / 16.0,
                         1.0 / 16.0, 1.0 / 16.0, 1.0 / 16.0, 1.0 / 16.0};

    double sigma_rt = constants::sigma_rt;
    double sigma_x = constants::sigma_x;
    double sigma_n = constants::sigma_n;

    int hole_width = 0;
};

void get_image_coordinates(int& x, int& y, const int idx);
int idx_from_coordinates(const int x, const int y, const int width);
bool is_out_of_bounds(const int x, const int y);
int clamp_x_coordinate(const int x);
int clamp_y_coordinate(const int y);

double compute_weight(const int p, const int q, const KernelData& kernel_data, const PixelBuffers& buffers);
int expand_kernel_idx(const int idx, const int hole_width);

vec3 blur_pixel(const int p, const KernelData& kernel_data, PixelBuffers& buffers);
void one_denoising_iteration(const KernelData& kernel_data, PixelBuffers& buffers);
void median_filter(PixelBuffers& buffers);
void denoise(PixelBuffers& buffers);
#endif
