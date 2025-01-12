
#ifndef DENOISE_H
#define DENOISE_H
#include "vec3.h"
#include "constants.h"


struct KernelData{
    double kernel[25] = {
        1.0/16.0, 1.0/16.0, 1.0/16.0, 1.0/16.0, 1.0/16.0, 
        1.0/16.0, 1.0/4.0,  1.0/4.0,  1.0/4.0,  1.0/16.0, 
        1.0/16.0, 1.0/4.0,  3.0/8.0,  1.0/4.0,  1.0/16.0, 
        1.0/16.0, 1.0/4.0,  1.0/4.0,  1.0/4.0,  1.0/16.0, 
        1.0/16.0, 1.0/16.0, 1.0/16.0, 1.0/16.0, 1.0/16.0};

    double sigma_rt = 10;
    double sigma_x = 10;
    double sigma_n = 10;
    int hole_width = 0;
};


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


double compute_weight(const int p, const int q, const KernelData& kernel_data, const vec3* pixel_buffer, const vec3* position_buffer, const vec3* normal_buffer){
    double w_rt = std::exp(-(pixel_buffer[p] - pixel_buffer[q]).length() / (kernel_data.sigma_rt * kernel_data.sigma_rt));
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


vec3 blur_pixel(const int p, const KernelData& kernel_data, const int iteration, const vec3* pixel_buffer, const vec3* position_buffer, const vec3* normal_buffer){
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
            new_pixel_value += kernel_value * pixel_buffer[q] * weight;
            normalization += kernel_value * weight;
        }
    }
    return new_pixel_value / normalization;

}


void one_denoising_iteration(const int iteration, const KernelData& kernel_data, vec3* pixel_buffer, const vec3* position_buffer, const vec3* normal_buffer){
    vec3* tmp_pixel_buffer = new vec3[constants::WIDTH*constants::HEIGHT];
    for (int j = 0; j < constants::WIDTH * constants::HEIGHT; j++){
            tmp_pixel_buffer[j] = blur_pixel(j, kernel_data, iteration, pixel_buffer, position_buffer, normal_buffer);
    }

    for (int j = 0; j < constants::WIDTH * constants::HEIGHT; j++){
            pixel_buffer[j] = tmp_pixel_buffer[j];
    }

    delete[] tmp_pixel_buffer;
    
}

void atrous_denoise(vec3* pixel_buffer, const vec3* position_buffer, const vec3* normal_buffer){
    KernelData kernel_data;
    for (int iteration = 0; iteration < constants::denoising_iterations; iteration++){
        one_denoising_iteration(iteration, kernel_data, pixel_buffer, position_buffer, normal_buffer);
        kernel_data.sigma_rt /= 2.0;
        kernel_data.sigma_x /= 2.0;
        kernel_data.sigma_n /= 2.0;
        kernel_data.hole_width += pow(2, iteration);
    }
}


void denoise(vec3* pixel_buffer, const vec3* position_buffer, const vec3* normal_buffer){
    if (!constants::use_denoising){
        return;
    }

    atrous_denoise(pixel_buffer, position_buffer, normal_buffer);
}

#endif