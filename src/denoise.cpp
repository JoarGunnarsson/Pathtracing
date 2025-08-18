
#include "denoise.h"


void get_image_coordinates(int& x, int& y, const int idx){
    x = idx % constants::HEIGHT;
    y = (idx - x) / constants::WIDTH;
}


int idx_from_coordinates(const int x, const int y, const int width){
    return width * y + x;
}


void clamp_x_coordinate(int& x){
    if (x < 0){
        x = -x;
    }
    else if (x > constants::WIDTH-1){
        x = 2 * (constants::WIDTH-1) - x;
    }
}


void clamp_y_coordinate(int& y){
    if (y < 0){
        y = -y;
    }
    else if (y > constants::HEIGHT-1){
        y = 2 * (constants::HEIGHT-1) - y;
    }
}


double compute_weight(const int p, const int q, const KernelData& kernel_data, const double* pixel_buffer, const vec3* position_buffer, const vec3* normal_buffer){
    vec3 pixel_p = vec3(pixel_buffer[3*p], pixel_buffer[3*p+1], pixel_buffer[3*p+2]);
    vec3 pixel_q = vec3(pixel_buffer[3*q], pixel_buffer[3*q+1], pixel_buffer[3*q+2]);
    double w_rt = std::exp(-(pixel_p - pixel_q).length() / (kernel_data.sigma_rt * kernel_data.sigma_rt));
    double w_x = std::exp(-(position_buffer[p] - position_buffer[q]).length() / (kernel_data.sigma_x * kernel_data.sigma_x));
    double w_n = std::exp(-(normal_buffer[p] - normal_buffer[q]).length() / (kernel_data.sigma_n * kernel_data.sigma_n));

    return w_rt * w_x * w_rt;
}


int expand_kernel_idx(const int idx, const int hole_width){
    if(idx == 1){
        return idx + hole_width;
    }
    else if(idx == -1){
        return idx - hole_width;
    }
    else if(idx == 2){
        return idx + 2 * hole_width;
    }
    else if(idx == -2){
        return idx - 2 * hole_width;
    }
    return 0;
}


vec3 blur_pixel(const int p, const KernelData& kernel_data, const int iteration, const double* pixel_buffer, const vec3* position_buffer, const vec3* normal_buffer){
    vec3 new_pixel_value = vec3(0,0,0);
    double normalization = 0;
    int global_x;
    int global_y;
    get_image_coordinates(global_x, global_y, p);
    for (int dx = -2; dx <= 2; dx++){
        for (int dy = -2; dy <= 2; dy++){
            int expanded_dx = expand_kernel_idx(dx, kernel_data.hole_width);
            int expanded_dy = expand_kernel_idx(dy, kernel_data.hole_width);
            int kernel_idx = idx_from_coordinates(dx+2, dy+2, 5);

            int x = global_x + expanded_dx;
            int y = global_y + expanded_dy;

            clamp_x_coordinate(x);
            clamp_y_coordinate(y);
            int q = idx_from_coordinates(x, y, constants::WIDTH);

            double weight = compute_weight(p, q, kernel_data, pixel_buffer, position_buffer, normal_buffer);
            double kernel_value = kernel_data.kernel[kernel_idx];
            vec3 pixel_color = vec3(pixel_buffer[3*q], pixel_buffer[3*q+1], pixel_buffer[3*q+2]);
            vec3 pixel_contribution = kernel_value * pixel_color * weight;
            if (isnan(pixel_contribution.length_squared())){
                pixel_contribution = vec3(0.0); //TODO: Did this fix the issue?
                weight = 0.0;
            }
            new_pixel_value += pixel_contribution;
            normalization += kernel_value * weight;
        }
    }
    return new_pixel_value / normalization;

}


void one_denoising_iteration(const int iteration, const KernelData& kernel_data, double* pixel_buffer, const vec3* position_buffer, const vec3* normal_buffer){
    double* tmp_pixel_buffer = new double[constants::WIDTH * constants::HEIGHT * 3];
    for (int j = 0; j < constants::WIDTH * constants::HEIGHT; j++){
            vec3 blurred_pixel = blur_pixel(j, kernel_data, iteration, pixel_buffer, position_buffer, normal_buffer);

            for (int i = 0; i < 3; i++){
                tmp_pixel_buffer[3*j+i] = blurred_pixel[i];
            }
    }

    for (int j = 0; j < constants::WIDTH * constants::HEIGHT * 3; j++){
        pixel_buffer[j] = tmp_pixel_buffer[j];
    }

    delete[] tmp_pixel_buffer;

}

void denoise(double* pixel_buffer, const vec3* position_buffer, const vec3* normal_buffer){
    KernelData kernel_data;
    for (int iteration = 0; iteration < constants::denoising_iterations; iteration++){
        one_denoising_iteration(iteration, kernel_data, pixel_buffer, position_buffer, normal_buffer);
        kernel_data.sigma_rt /= 2.0;
        kernel_data.sigma_x /= 2.0;
        kernel_data.sigma_n /= 2.0;
        kernel_data.hole_width += pow(2, iteration);
    }
}
