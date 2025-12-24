
#include "utils.h"
#include "denoise.h"

void get_image_coordinates(int& x, int& y, const int idx) {
    x = idx % constants::WIDTH;
    y = idx / constants::WIDTH;
}

int idx_from_coordinates(const int x, const int y, const int width) {
    return width * y + x;
}

bool is_out_of_bounds(const int x, const int y) {
    return x < 0 || y < 0 || x > constants::WIDTH - 1 || y > constants::HEIGHT - 1;
}

int clamp_x_coordinate(int x) {
    if (x < 0) {
        return -x;
    }
    else if (x > constants::WIDTH - 1) {
        return 2 * (constants::WIDTH - 1) - x;
    }
    return x;
}

int clamp_y_coordinate(int y) {
    if (y < 0) {
        return -y;
    }
    else if (y > constants::HEIGHT - 1) {
        return 2 * (constants::HEIGHT - 1) - y;
    }
    return y;
}

double compute_weight(const int p, const int q, const KernelData& kernel_data, const PixelBuffers& buffers) {
    vec3 pixel_p = vec3(buffers.image[3 * p], buffers.image[3 * p + 1], buffers.image[3 * p + 2]);
    vec3 pixel_q = vec3(buffers.image[3 * q], buffers.image[3 * q + 1], buffers.image[3 * q + 2]);
    double w_rt = std::exp(-(pixel_p - pixel_q).length() / (kernel_data.sigma_rt * kernel_data.sigma_rt));
    double w_x = std::exp(-(buffers.position_buffer[p] - buffers.position_buffer[q]).length() /
                          (kernel_data.sigma_x * kernel_data.sigma_x));
    double w_n = std::exp(-(buffers.normal_buffer[p] - buffers.normal_buffer[q]).length() /
                          (kernel_data.sigma_n * kernel_data.sigma_n));

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

vec3 blur_pixel(const int p, const KernelData& kernel_data, PixelBuffers& buffers) {
    vec3 new_pixel_value = vec3(0, 0, 0);
    double normalization = 0;
    int global_x;
    int global_y;
    get_image_coordinates(global_x, global_y, p);
    for (int dx = -2; dx <= 2; dx++) {
        for (int dy = -2; dy <= 2; dy++) {
            int expanded_dx = expand_kernel_idx(dx, kernel_data.hole_width);
            int expanded_dy = expand_kernel_idx(dy, kernel_data.hole_width);
            int kernel_idx = idx_from_coordinates(dx + 2, dy + 2, 5);

            int x = global_x + expanded_dx;
            int y = global_y + expanded_dy;

            x = clamp_x_coordinate(x);
            y = clamp_y_coordinate(y);
            int q = idx_from_coordinates(x, y, constants::WIDTH);

            double weight = compute_weight(p, q, kernel_data, buffers);
            double kernel_value = kernel_data.kernel[kernel_idx];
            vec3 pixel_color = vec3(buffers.image[3 * q], buffers.image[3 * q + 1], buffers.image[3 * q + 2]);
            vec3 pixel_contribution = kernel_value * pixel_color * weight;
            if (isnan(pixel_contribution.length_squared())) {
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
    size_t buffer_size = static_cast<size_t>(constants::WIDTH * constants::HEIGHT * 3);
    double* tmp_image = new double[buffer_size];

    for (int j = 0; j < constants::WIDTH * constants::HEIGHT; j++) {
        vec3 blurred_pixel = blur_pixel(j, kernel_data, buffers);

        for (int i = 0; i < 3; i++) {
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
    size_t buffer_size = static_cast<size_t>(constants::WIDTH * constants::HEIGHT * 3);
    double* tmp_image = new double[buffer_size];

    int offset = (constants::median_kernel_size - 1) / 2;
    size_t size = static_cast<size_t>(constants::median_kernel_size * constants::median_kernel_size);
    double* r = new double[size];
    double* g = new double[size];
    double* b = new double[size];
    for (int idx = 0; idx < constants::WIDTH * constants::HEIGHT; idx++) {
        int x, y;
        get_image_coordinates(x, y, idx);

        for (int dx = -offset; dx <= offset; dx++) {
            for (int dy = -offset; dy <= offset; dy++) {
                int local_x = clamp_x_coordinate(x + dx);
                int local_y = clamp_y_coordinate(y + dy);

                int neighbor_idx = idx_from_coordinates(local_x, local_y, constants::WIDTH);
                int local_idx = idx_from_coordinates(dx + offset, dy + offset, constants::median_kernel_size);

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
