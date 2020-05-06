#ifndef RAY_TRACING_TRACER_H
#define RAY_TRACING_TRACER_H

#include <vector>     // std::vector
#include <string>     // std::string
#include <limits>     // std::numeric_limits<float>::max()
#include <cmath>      // M_PI, std::fmax(), std::log(), std::exp(), tanf
#include <omp.h>      // #pragma
#include "objects.h"  // Material, Object
#include "vector.h"   // float3
#include "image.h"    // get_cubemap_color(), get_background_color(), save_image()


struct Config {
    uint height            = 512;
    uint width             = 512;
    uint maxDepth          = 3;
    const float fov        = M_PI / 3.0f; // fieldOfView
    bool antiAliasing      = false;
    bool toneMapping       = false;
    bool isBackgroundMap   = false;
    bool isCubemap         = false;
    float3 cameraOrigin;
    float3 backgroundColor;
    std::string outFile;
} config;


inline float3 reflect(const float3 &D, const float3 &N) { return D - N * 2.0f * (D * N); }


float3 refract(const float3 &D, float3 N, float eta_t, float eta_d=1.0f) {
    float eta = eta_d / eta_t;
    float cosI = - D * N;
    if (cosI < 0) {
        cosI = - cosI, N = - N, eta = 1.0f / eta;
    }
    float sinT2 = eta * eta * (1.0f - cosI * cosI);
    float cosT = sqrtf(1.0f - sinT2);
    return eta * D + (eta * cosI - cosT) * N;
}


float3 get_background_color(const float3 &dir) {
    float3 absDir = dir.abs();
    float u, v, max = absDir.max();
    if (config.isCubemap) {
        CUBEMAP_POSITION pos;
        if (absDir.x == max) {
            if (dir.x > 0) { pos = X_POS, u = -dir.y, v = dir.z; }
            else { pos = X_NEG, u = -dir.y, v = dir.z; }
        } else if (absDir.y == max) {
            if (dir.y > 0) { pos = Y_POS, u = dir.z, v = dir.x; }
            else { pos = Y_NEG, v = dir.x, u = -dir.z; };
        } else if (absDir.z == max) {
            if (dir.z > 0) { pos = Z_POS, u = -dir.y, v = dir.x; }
            else { pos = Z_NEG, u = -dir.y, v = dir.x; }
        }
        // [-1, 1] to [0, 1]
        u = 0.5f * (u / max + 1.0f), v = 0.5f * (v / max + 1.0f);
        return get_cubemap_color(pos, u, v);
    } else if (config.isBackgroundMap && absDir.z == max) {
        u = -dir.y / absDir.z, v = dir.x / absDir.z;
        u = 0.5f * (u + 1.0f), v = 0.5f * (v + 1.0f);
        return get_background_color(u, v);
    }
    return config.backgroundColor;
}


bool trace(const Ray &ray, const std::vector<Object *> &objects,
           float3 &hitPoint, float3 &hitNormal, float3 &hitColor,
           const Material **hitMaterial) {
    float curDist, closestDist = std::numeric_limits<float>::max();
    float3 P, N, C;
    for (const auto &objectPtr : objects) {
        if (objectPtr->is_ray_intersect(ray, curDist, P, N, C) && curDist < closestDist) {
            hitPoint  = P, hitNormal = N, hitColor  = C;
            *hitMaterial = &objectPtr->material;
            closestDist = curDist;
        }
    }
    return closestDist < 1024;
}


/* backward Ray Tracing */
float3 cast_ray(
        const Ray &ray, const std::vector<Object *> &objects,
        const std::vector<Light> &lights, uint depth=0) {
    float3 hitPoint, N, color;
    const Material *hitMaterial = nullptr;
    const float eps = 1e-3;

    if (depth >= config.maxDepth || !trace(ray, objects, hitPoint, N, color, &hitMaterial)) {
        return get_background_color(ray.direction);  // ray miss
    }
    if (hitMaterial->materialType == NORMAL) {
        return N;
    }
    float3 reflectRefract(0.0f, 0.0f, 0.0f);

    if (hitMaterial->materialType == MIRROR) {
        float3 reflectDir  = reflect(ray.direction, N).normalize();
        float3 reflectOrig = (reflectDir * N < 0.0f) ? hitPoint - N * eps : hitPoint + N * eps;
        Ray reflectRay(reflectOrig, reflectDir);
        float3 reflectColor = cast_ray(reflectRay, objects, lights, depth + 1);
        reflectRefract += reflectColor * hitMaterial->albedo[2];
    }
    if (hitMaterial->materialType == GLASS) {
        float3 refractDir  = refract(ray.direction, N, hitMaterial->eta).normalize();
        float3 refractOrig = (refractDir * N < 0.0f) ? hitPoint - N * eps : hitPoint + N * eps;
        Ray refractRay(refractOrig, refractDir);
        float3 refractColor = cast_ray(refractRay, objects, lights, depth + 1);
        reflectRefract += refractColor * hitMaterial->albedo[3];
    }

    float lightIntensity = 0.0f, specularItensity = 0.0f;
    for (const auto &light : lights) {
        float3 lightDir   = (light.position - hitPoint).normalize();
        float3 shadowOrig = (lightDir * N < 0.0f) ? hitPoint - N * eps : hitPoint + N * eps;
        Ray shadowRay(shadowOrig, lightDir);
        float lightDistance = (light.position - hitPoint).norm();
        float3 shadowHitPoint, tmpN, tmpC;
        const Material *tmpPtr = nullptr;
        if (trace(shadowRay, objects, shadowHitPoint, tmpN, tmpC, &tmpPtr) &&
            (shadowHitPoint - shadowOrig).norm() < lightDistance) { continue; }

        lightIntensity += light.intensity * std::max(0.0f, lightDir * N);
        if (hitMaterial->materialType != LAMBERT) {
            specularItensity += light.intensity *
                                powf(std::max(0.0f, reflect(lightDir, N) * ray.direction), hitMaterial->shine);
        }
    }

    return color * lightIntensity * hitMaterial->albedo[0] +
           specularItensity * hitMaterial->albedo[1] + reflectRefract;
}


void tone_mapping(std::vector<float3> &image) {
    std::vector<float> img(image.size());
    float average = 0.0f, eps = 0.001f, white = 0.0f;
    for (uint i = 0; i < img.size(); ++i) {
        img[i] = 0.212671f * image[i].x +
                 0.715600f * image[i].y +
                 0.072169f * image[i].z;
        average += std::log(eps + img[i]);
        white = std::fmax(white, img[i]);
    }
    average = 0.8f / std::exp(average / img.size());
    white = white * white;
    for (uint i = 0; i < image.size(); ++i) {
        image[i] *= average * ((1.0f + average * img[i] / white) / (1.0f + average * img[i]));
    }
}


float3 transform(float x, float y) {
    static float angle       = tanf(config.fov / 2.0f);
    static float aspectRatio = (float)config.width / (float)config.height;
    float u =  (2.0f * (x + 0.5f) / config.width  - 1.0f) * angle * aspectRatio;
    float v = -(2.0f * (y + 0.5f) / config.height - 1.0f) * angle;
    return float3(u, v, -1).normalize();
//    --- or ---
//    float u =   x - config.width  / 2.0f;
//    float v = -(y - config.height / 2.0f);
//    float w = -(config.height / 2.0f) / tanf(config.fov / 2.0f);
//    return float3(u, v, w).normalize();
}


void render(const std::vector<Object *> &objects, const std::vector<Light> &lights) {
    std::vector<float3> image(config.width * config.height);

    #pragma omp parallel for
    for (uint y = 0; y < config.height; ++y) {
        for (uint x = 0; x < config.width; ++x) {
            image[y * config.width + x] = cast_ray(Ray(config.cameraOrigin, transform(x, y)), objects, lights);
            if (config.antiAliasing) {
                const float bias = 0.5f;
                image[y * config.width + x] +=
                        cast_ray(Ray(config.cameraOrigin, transform((float)x + bias, y)), objects, lights) +
                        cast_ray(Ray(config.cameraOrigin, transform((float)x - bias, y)), objects, lights) +
                        cast_ray(Ray(config.cameraOrigin, transform(x, (float)y + bias)), objects, lights) +
                        cast_ray(Ray(config.cameraOrigin, transform(x, (float)y - bias)), objects, lights);
                image[y * config.width + x] /= 5.0f;
            }
        }
    }
    if (config.toneMapping) {
        tone_mapping(image);
    }
    save_image(image, config.width, config.height, config.outFile);
}


#endif // RAY_TRACING_TRACER_H
