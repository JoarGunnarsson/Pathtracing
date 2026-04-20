#include <iostream>
#include <fstream>
#include <cstring>
#include <stdexcept>
#include "denoise.h"
#include "utils.h"
#include "constants.h"
#include "scene.h"

int main(int argc, char* argv[]) {
    if (argc != 2) {
        throw std::runtime_error(
            "Invalid arguments provided.\n"
            "Usage: main scene_directory\n\n"

            "positional arguments:\n"
            "   scene_directory                  path to the scene directory (relative to main project directory)");
    }

    load_settings(std::string(argv[1]) + "/settings.json");

    if (!constants::enable_atrous_filtering && !constants::enable_median_filtering) {
        std::clog << "No denoising mode set in scene settings, exiting..." << std::endl;
        return 0;
    }

    size_t pixel_count = constants::WIDTH * constants::HEIGHT;
    size_t FILESIZE = pixel_count * 3 * sizeof(double);

    PixelBuffers pixel_buffers;
    int image_fd, position_fd, normal_fd;
    pixel_buffers.image = create_mmap(constants::raw_pixeldata_file, FILESIZE, false, image_fd);
    pixel_buffers.position_buffer = create_mmap(constants::raw_positiondata_file, FILESIZE, false, position_fd);
    pixel_buffers.normal_buffer = create_mmap(constants::raw_normaldata_file, FILESIZE, false, normal_fd);

    std::clog << "Denoising..." << std::endl;
    PixelBuffers denoising_buffers;
    int denoised_image_fd;
    denoising_buffers.image = create_mmap(constants::raw_denoised_file_name, FILESIZE, true, denoised_image_fd);

    std::memcpy(denoising_buffers.image, pixel_buffers.image, FILESIZE);
    denoising_buffers.position_buffer = pixel_buffers.position_buffer;
    denoising_buffers.normal_buffer = pixel_buffers.normal_buffer;

    denoise(denoising_buffers);

    close_mmap(denoising_buffers.image, FILESIZE, denoised_image_fd);
    close_mmap(pixel_buffers.image, FILESIZE, image_fd);
    close_mmap(pixel_buffers.position_buffer, FILESIZE, position_fd);
    close_mmap(pixel_buffers.normal_buffer, FILESIZE, normal_fd);

    return 0;
}
