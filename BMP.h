#ifndef BMP
#define BMP

#include <iostream>
#include <fstream>
#include <cstdint>
#include <vector>

// Ensure that the structs are packed with no padding
#pragma pack(push, 1)

// BMP file header
struct BMPFileHeader {
    uint16_t file_type{0x4D42}; // "BM"
    uint32_t file_size{0};      // Size of the file (in bytes)
    uint16_t reserved1{0};      // Reserved
    uint16_t reserved2{0};      // Reserved
    uint32_t offset_data{54};   // Start position of pixel data (bytes from the beginning of the file)
};

// BMP information header
struct BMPInfoHeader {
    uint32_t size{40};          // Size of this header (in bytes)
    int32_t width{0};           // Width of the bitmap (in pixels)
    int32_t height{0};          // Height of the bitmap (in pixels)
    uint16_t planes{1};         // Number of color planes
    uint16_t bit_count{24};     // Number of bits per pixel
    uint32_t compression{0};    // Compression type
    uint32_t size_image{0};     // Size of the image data (in bytes)
    int32_t x_pixels_per_meter{0};
    int32_t y_pixels_per_meter{0};
    uint32_t colors_used{0};    // Number of colors used
    uint32_t colors_important{0}; // Number of important colors
};

#pragma pack(pop)

// Function to save BMP
void saveBMP(const char* filename, int width, int height, std::vector<float> colorArray) {
    BMPFileHeader file_header;
    BMPInfoHeader info_header;

    info_header.width = width;
    info_header.height = height;
    file_header.file_size = sizeof(BMPFileHeader) + sizeof(BMPInfoHeader) + width * height * 3;

    std::ofstream file(filename, std::ios::binary);
    if (!file) {
        std::cerr << "Error: Could not open file for writing: " << filename << std::endl;
        return;
    }

    // Write headers
    file.write(reinterpret_cast<const char*>(&file_header), sizeof(file_header));
    file.write(reinterpret_cast<const char*>(&info_header), sizeof(info_header));

    // Write pixel data (BGR format)
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int startIndex = 3*(y*height+x);
            int r = (int) (colorArray[startIndex] * 255);
            int g = (int) (colorArray[startIndex+1] * 255);
            int b =  (int) (colorArray[startIndex+2] * 255);
            file.put(b).put(g).put(r); // BMP uses BGR format
        }
    }

    file.close();
    std::cout << "Image saved successfully to " << filename << std::endl;
}
#endif
