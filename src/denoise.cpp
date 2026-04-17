
#include "utils.h"
#include "denoise.h"
#include <cstring>
#include <algorithm>

void get_image_coordinates(int& x, int& y, const size_t idx) {
    x = static_cast<int>(idx % constants::WIDTH);
    y = static_cast<int>(idx / constants::WIDTH);
}

size_t idx_from_coordinates(const size_t x, const size_t y, const size_t width) {
    return width * y + x;
}

bool is_out_of_bounds(const int x, const int y) {
    return x < 0 || y < 0 || x > static_cast<int>(constants::WIDTH) - 1 || y > static_cast<int>(constants::HEIGHT) - 1;
}

int clamp_x_coordinate(int x) {
    if (x < 0) {
        return -x;
    }
    else if (x > static_cast<int>(constants::WIDTH) - 1) {
        return 2 * (static_cast<int>(constants::WIDTH) - 1) - x;
    }
    return x;
}

int clamp_y_coordinate(int y) {
    if (y < 0) {
        return -y;
    }
    else if (y > static_cast<int>(constants::HEIGHT) - 1) {
        return 2 * (static_cast<int>(constants::HEIGHT) - 1) - y;
    }
    return y;
}

double compute_weight(const size_t p, const size_t q, const KernelData& kernel_data, const PixelBuffers& buffers) {
    vec3 pixel_p, pixel_q;
    vec3 pos_p, pos_q;
    vec3 normal_p, normal_q;
    for (size_t j = 0; j < 3; j++) {
        pixel_p[j] = buffers.image[3 * p + j];
        pixel_q[j] = buffers.image[3 * q + j];

        pos_p[j] = buffers.position_buffer[3 * p + j];
        pos_q[j] = buffers.position_buffer[3 * q + j];

        normal_p[j] = buffers.normal_buffer[3 * p + j];
        normal_q[j] = buffers.normal_buffer[3 * q + j];
    }
    double w_rt = std::exp(-(pixel_p - pixel_q).length() / (kernel_data.sigma_rt * kernel_data.sigma_rt));
    double w_x = std::exp(-(pos_p - pos_q).length() / (kernel_data.sigma_x * kernel_data.sigma_x));
    double w_n = std::exp(-(normal_p - normal_q).length() / (kernel_data.sigma_n * kernel_data.sigma_n));

    return w_rt * w_x * w_n;
}

int expand_kernel_idx(const int idx, const int hole_width) {
    if (idx == 1) {
        return idx + hole_width;
    }
    else if (idx == -1) {
        return idx - hole_width;
    }
    else if (idx == 2) {
        return idx + 2 * hole_width;
    }
    else if (idx == -2) {
        return idx - 2 * hole_width;
    }
    return 0;
}

vec3 blur_pixel(const size_t p, const KernelData& kernel_data, PixelBuffers& buffers) {
    vec3 new_pixel_value = vec3(0, 0, 0);
    double normalization = 0;
    int global_x;
    int global_y;
    get_image_coordinates(global_x, global_y, p);
    for (int dx = -2; dx <= 2; dx++) {
        for (int dy = -2; dy <= 2; dy++) {
            int expanded_dx = expand_kernel_idx(dx, kernel_data.hole_width);
            int expanded_dy = expand_kernel_idx(dy, kernel_data.hole_width);
            size_t kernel_idx = idx_from_coordinates(static_cast<size_t>(dx + 2), static_cast<size_t>(dy + 2), 5);

            int x = global_x + expanded_dx;
            int y = global_y + expanded_dy;

            x = clamp_x_coordinate(x);
            y = clamp_y_coordinate(y);
            size_t q = idx_from_coordinates(static_cast<size_t>(x), static_cast<size_t>(y), constants::WIDTH);

            double weight = compute_weight(p, q, kernel_data, buffers);
            double kernel_value = kernel_data.kernel[kernel_idx];
            vec3 pixel_color = vec3(buffers.image[3 * q], buffers.image[3 * q + 1], buffers.image[3 * q + 2]);
            vec3 pixel_contribution = kernel_value * pixel_color * weight;
            if (std::isnan(pixel_contribution.length_squared())) {
                pixel_contribution = vec3(0.0);
                weight = 0.0;
            }
            new_pixel_value += pixel_contribution;
            normalization += kernel_value * weight;
        }
    }
    return new_pixel_value / normalization;
}

void one_denoising_iteration(const KernelData& kernel_data, PixelBuffers& buffers) {
    size_t buffer_size = constants::WIDTH * constants::HEIGHT * 3;
    double* tmp_image = new double[buffer_size];

    for (size_t j = 0; j < constants::WIDTH * constants::HEIGHT; j++) {
        vec3 blurred_pixel = blur_pixel(j, kernel_data, buffers);

        for (size_t i = 0; i < 3; i++) {
            tmp_image[3 * j + i] = blurred_pixel[i];
        }
    }

    std::memcpy(buffers.image, tmp_image, buffer_size * sizeof(double));

    delete[] tmp_image;
}

void atrous_filter(PixelBuffers& buffers) {
    KernelData kernel_data;
    for (int iteration = 0; iteration < constants::denoising_iterations; iteration++) {
        one_denoising_iteration(kernel_data, buffers);
        kernel_data.sigma_rt /= 2.0;
        kernel_data.sigma_x /= 2.0;
        kernel_data.sigma_n /= 2.0;
        kernel_data.hole_width += static_cast<int>(pow(2, iteration));
    }
}

void median_filter(PixelBuffers& buffers) {
    size_t buffer_size = constants::WIDTH * constants::HEIGHT * 3;
    double* tmp_image = new double[buffer_size];

    int offset = (constants::median_kernel_size - 1) / 2;
    size_t size = static_cast<size_t>(constants::median_kernel_size * constants::median_kernel_size);
    double* r = new double[size];
    double* g = new double[size];
    double* b = new double[size];
    for (size_t idx = 0; idx < constants::WIDTH * constants::HEIGHT; idx++) {
        int x, y;
        get_image_coordinates(x, y, idx);

        for (int dx = -offset; dx <= offset; dx++) {
            for (int dy = -offset; dy <= offset; dy++) {
                int local_x = clamp_x_coordinate(x + dx);
                int local_y = clamp_y_coordinate(y + dy);

                size_t neighbor_idx =
                    idx_from_coordinates(static_cast<size_t>(local_x), static_cast<size_t>(local_y), constants::WIDTH);
                size_t local_idx =
                    idx_from_coordinates(static_cast<size_t>(dx + offset), static_cast<size_t>(dy + offset),
                                         static_cast<size_t>(constants::median_kernel_size));

                r[local_idx] = buffers.image[3 * neighbor_idx];
                g[local_idx] = buffers.image[3 * neighbor_idx + 1];
                b[local_idx] = buffers.image[3 * neighbor_idx + 2];
            }
        }

        std::nth_element(r, r + size / 2, r + size);
        std::nth_element(g, g + size / 2, g + size);
        std::nth_element(b, b + size / 2, b + size);

        tmp_image[3 * idx] = r[size / 2];
        tmp_image[3 * idx + 1] = g[size / 2];
        tmp_image[3 * idx + 2] = b[size / 2];
    }
    std::memcpy(buffers.image, tmp_image, buffer_size * sizeof(double));
    delete[] r;
    delete[] g;
    delete[] b;
    delete[] tmp_image;
}

void denoise(PixelBuffers& buffers) {
    if (constants::enable_median_filtering) {
        median_filter(buffers);
    }

    if (constants::enable_atrous_filtering) {
        atrous_filter(buffers);
    }
}
