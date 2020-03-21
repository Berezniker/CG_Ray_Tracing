#ifndef RAY_TRACING_IMAGE_H
#define RAY_TRACING_IMAGE_H

#include <iostream>  // std::cerr, std::endl, std::max(), std::min()
#include <fstream>   // std::ofstream
#include <vector>    // std::vector
#include <string>    // std::string
#include <cstdlib>   // exit()
#include "vector.h"  // float3

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"  // stbi_write_jpg()
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"        // stbi_load(), stbi_image_free()

enum CUBEMAP_POSITION {X_POS, X_NEG, Y_POS, Y_NEG, Z_POS, Z_NEG};
std::vector<std::vector<float3>> CUBEMAP(6);
int cubemapSize;

std::vector<float3> backgroundImage;
std::vector<int> backgroundImageSize(2);


float3 hex_to_rgb(uint32_t hexValue) {
    return float3((float)((hexValue >> 16) & 0xFF) / 255.0f,
                  (float)((hexValue >>  8) & 0xFF) / 255.0f,
                  (float)( hexValue        & 0xFF) / 255.0f);
}


void load_image(const char *imagePath, std::vector<float3> &image, int &width, int &height) {
    int nChannels = -1;
    unsigned char *STBimage = stbi_load(imagePath, &width, &height, &nChannels, 0);
    if (!STBimage || nChannels != 3) {
        std::cerr << "[Error]: can not load the image: \"" << imagePath  << "\"" << std::endl;
        exit(-1);
    }
    image = std::vector<float3>(width * height);
    for (uint y = 0; y < height ; ++y) {
        for (uint x = 0; x < width; ++x) {
            image[y * width + x] = float3(STBimage[(y * width + x) * 3 + 0],
                                          STBimage[(y * width + x) * 3 + 1],
                                          STBimage[(y * width + x) * 3 + 2]) / 255.0f;
        }
    }
    stbi_image_free(STBimage);
}


void load_background_image(const char *imagePath) {
    load_image(imagePath, backgroundImage, backgroundImageSize[1], backgroundImageSize[0]);
}


void load_cubemap(const std::vector<std::string> &imagePaths) {
    if (imagePaths.size() != 6) {
        std::cerr << "[!] Error: cubemap should consist of 6 images" << std::endl;
        exit(-1);
    }
    for (uint i = 0; i < imagePaths.size(); ++i) {
        load_image(imagePaths[i].c_str(), CUBEMAP[i], cubemapSize, cubemapSize);
    }
}


float3 get_background_color(float u, float v) {
    int x = (int)(u * (float)backgroundImageSize[1]);
    int y = (int)(v * (float)backgroundImageSize[0]);
    return backgroundImage[y * backgroundImageSize[1] + x];
}


float3 get_cubemap_color(CUBEMAP_POSITION pos, float u, float v) {
    int x = (int)(u * (float)cubemapSize), y = (int)(v * (float)cubemapSize);
    return CUBEMAP[pos][x * cubemapSize + y];
}


inline float clip(float x) { return std::max(0.0f, std::min(x, 1.0f)); }


void save_image(const std::vector<float3> &image, uint width, uint height, const std::string &outFile) {
    if (outFile.find("jpg") != std::string::npos) {
        std::vector<unsigned char> STBimage(width * height * 3);
        for (uint i = 0; i < height * width; ++i) {
            for (uint j = 0; j < 3; ++j) {
                STBimage[i * 3 + j] = (unsigned char)(255.0f * clip(image[i][j]));
            }
        }
        stbi_write_jpg(outFile.c_str(), width, height, 3, STBimage.data(), 100);
    }
    else if (outFile.find("ppm") != std::string::npos) {
        std::ofstream ofs;
        ofs.open(outFile.c_str());
        ofs << "P6\n" << width << " " << height << "\n255\n";
        for (uint i = 0; i < height * width; ++i) {
            for (uint j = 0; j < 3; ++j) {
                ofs << (unsigned char)(255.0f * clip(image[i][j]));
            }
        }
        ofs.close();
    } else {
        std::cerr << "[Error]: Can't save the image with the specified extension:"
                  << "\"" << outFile << "\"" << std::endl;
    }
}


#endif // RAY_TRACING_IMAGE_H
