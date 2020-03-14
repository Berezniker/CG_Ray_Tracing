#include <unordered_map>  // std::unordered_map
#include <iostream>       // std::max(), std::min(), std::fabs()
#include <cstdlib>        // srand(), rand()
#include <vector>         // std::vector
#include <string>         // std::string
#include <limits>         // std::numeric_limits<float>::max()
#include <cmath>          // sqrtf(), atanf2(), M_PI
#include <ctime>          // clock(), CLOCKS_PER_SEC
#include <omp.h>          // omp_set_num_threads()
#include "objects.h"
#include "vector.h"
#include "image.h"
#include "mesh.h"


struct Config {
    uint height            = 512;
    uint width             = 512;
    const float fov        = M_PI / 3.0f; // fieldOfView
    const uint maxDepth    = 3;
    bool antiAliasing      = false;
    bool isBackgroundMap   = false;
    bool isCubemap         = false;
    float3 cameraOrigin    = float3(0.0f);
    float3 backgroundColor = float3(0.0f);
    std::string outFile;
} config;


inline float3 reflect(const float3 &D, const float3 &N) { return D - N * 2.0f * (D * N); }


float3 refract(const float3 &I, const float3 &N, float eta_t, float eta_i=1.0f) {
    float eta = eta_i / eta_t;
    float cosI = - I * N;
    float sinT2 = eta * eta * (1.0f - cosI * cosI);
    if (sinT2 > 1.0f) return float3(1.0f);
    float cosT = sqrtf(1.0f - sinT2);
    return eta * I + (eta * cosI - cosT) * N;
}


float3 get_backgound_color(const float3 &dir) {
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
        u = 0.5f * (u / max + 1.0f), v = 0.5f * (v / max + 1.0f);  // [-1, 1] to [0, 1]
        return get_cubemap_color(pos, u, v);
    } else if (config.isBackgroundMap) {
        if (absDir.z == max) {
            u = -dir.y / absDir.z, v = dir.x / absDir.z;
            u = 0.5f * (u + 1.0f), v = 0.5f * (v + 1.0f);
            return get_background_color(u, v);
        }
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


/* backward ray-tracing */
float3 castRay(
        const Ray &ray, const std::vector<Object *> &objects,
        const std::vector<Light> &lights, uint depth=0) {
    float3 hitPoint, N, color;
    const Material *hitMaterial = nullptr;
    static float eps = 1e-3;

    if (depth >= config.maxDepth || !trace(ray, objects, hitPoint, N, color, &hitMaterial)) {
        return get_backgound_color(ray.direction);  // ray miss
    }
    if (hitMaterial->materialType == NORMAL) {
        return N;
    }
    float3 reflectRefract(0.0f, 0.0f, 0.0f);

    if (hitMaterial->materialType == MIRROR) {
        float3 reflectDir = reflect(ray.direction, N).normalize();
        float3 reflectOrig = (reflectDir * N < 0) ? hitPoint - N * eps : hitPoint + N * eps;
        Ray reflectRay(reflectOrig, reflectDir);
        float3 reflectColor = castRay(reflectRay, objects, lights, depth + 1);
        reflectRefract += reflectColor * hitMaterial->albedo[2];
    }
    if (hitMaterial->materialType == GLASS) {
        float3 refractDir = refract(ray.direction, N, hitMaterial->eta).normalize();
        float3 refractOrig = (refractDir * N < 0) ? hitPoint - N * eps : hitPoint + N * eps;
        Ray refractRay(refractOrig, refractDir);
        float3 refractColor = castRay(refractRay, objects, lights, depth + 1);
        reflectRefract += refractColor * hitMaterial->albedo[3];
    }

    float diffuse_light_intensity = 0.0f, specular_light_intensity = 0.0f;
    for (const auto &light : lights) {
        float3 lightDir = (light.position - hitPoint).normalize();
        float3 shadowOrig = (lightDir * N < 0) ? hitPoint - N * eps : hitPoint + N * eps;
        Ray shadowRay(shadowOrig, lightDir);
        float lightDistance = (light.position - hitPoint).norm();
        float3 shadowHitPoint, tmpN, tmpC;
        const Material *tmpPtr = nullptr;
        if (trace(shadowRay, objects, shadowHitPoint, tmpN, tmpC, &tmpPtr) &&
           (shadowHitPoint - shadowOrig).norm() < lightDistance) { continue; }
        diffuse_light_intensity += light.intensity * std::max(0.0f, lightDir * N);
        if (hitMaterial->materialType != LAMBERT) {
            specular_light_intensity += light.intensity *
                    powf(std::max(0.0f, reflect(lightDir, N) * ray.direction), hitMaterial->specularPower);
        }
    }

//    return N;
//    return float3(N.norm());

    return color  * diffuse_light_intensity * hitMaterial->albedo[0] +
           specular_light_intensity * hitMaterial->albedo[1] + reflectRefract;
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

//    #pragma omp parallel for shared(image)
    for (uint y = 0; y < config.height; ++y) {
        for (uint x = 0; x < config.width; ++x) {
//            std::cout << "( " << y << ", " << x << " )" << std::endl;
            image[y * config.width + x] = castRay(Ray(config.cameraOrigin, transform(x, y)), objects, lights);
            if (config.antiAliasing) {
                static float bias = 0.5f;
                image[y * config.width + x] +=
                        castRay(Ray(config.cameraOrigin, transform((float)x + bias, y)), objects, lights) +
                        castRay(Ray(config.cameraOrigin, transform((float)x - bias, y)), objects, lights) +
                        castRay(Ray(config.cameraOrigin, transform(x, (float)y + bias)), objects, lights) +
                        castRay(Ray(config.cameraOrigin, transform(x, (float)y - bias)), objects, lights);
                image[y * config.width + x] /= 5.0f;
            }
        }
    }
    save_image(image, config.width, config.height, config.outFile.c_str());
}


void scene_1() {
    std::cout << "scene 1" << std::endl;
    // environment setting
    config.width  = 720;
    config.height = 512;
    load_background_image("../env/stars_milky_way.jpg");
    config.isBackgroundMap = true;
    config.antiAliasing    = true;
    if (config.outFile.empty()) config.outFile = "../out_1_scene.jpg";

    std::vector<Object *> objects;
    std::vector<Light>    lights;

    lights.emplace_back(float3(-30, 30 ,  30), 0.9f);
    lights.emplace_back(float3(30 , -25,  40), 0.1f);

    Material mercury(FONG, "../env/planet/mercury.jpg",  5.0f, {1.0f, 0.1f});
    Material   venus(FONG, "../env/planet/venus.jpg",    5.0f, {1.0f, 0.1f});
    Material   earth(FONG, "../env/planet/earth.jpg",    5.0f, {1.0f, 0.1f});
    Material    moon(FONG, "../env/planet/moon.jpg",     5.0f, {1.0f, 0.1f});
    Material    mars(FONG, "../env/planet/mars.jpg",     5.0f, {1.0f, 0.1f});
    Material jupiter(FONG, "../env/planet/jupiter.jpg",  5.0f, {1.0f, 0.1f});
    Material  saturn(FONG, "../env/planet/saturn.jpg",   5.0f, {1.0f, 0.1f});
    Material  uranus(FONG, "../env/planet/uranus.jpg",   5.0f, {1.0f, 0.1f});
    Material neptune(FONG, "../env/planet/neptune.jpg",  5.0f, {1.0f, 0.1f});
    Material    make(FONG, "../env/planet/makemake.jpg", 5.0f, {1.0f, 0.1f});
    Material   ceres(FONG, "../env/planet/ceres.jpg",    5.0f, {1.0f, 0.1f});
    Material    eris(FONG, "../env/planet/eris.jpg",     5.0f, {1.0f, 0.1f});
    Material     sun(FONG, "../env/planet/sun.jpg",      5.0f, {1.0f, 0.1f});

    objects.push_back(new Sphere(mercury, float3(-6.0f , -10.0f, -25.0f), 2.0f ));
    objects.push_back(new Sphere(  venus, float3(-10.0f, -10.0f, -28.0f), 3.5f ));
    objects.push_back(new Sphere(  earth, float3(-15.0f, -9.0f , -31.0f), 4.0f ));
    objects.push_back(new Sphere(   moon, float3(-15.0f, -9.0f , -25.0f), 1.0f ));
    objects.push_back(new Sphere(   mars, float3(-18.0f, -5.0f , -34.0f), 4.0f ));
    objects.push_back(new Sphere(jupiter, float3(-15.0f, 0.0f  , -51.0f), 16.0f));
    objects.push_back(new Sphere( saturn, float3(-5.0f , 8.0f  , -62.0f), 17.0f));
    objects.push_back(new Sphere( uranus, float3(13.0f , 15.0f , -75.0f), 17.0f));
    objects.push_back(new Sphere(neptune, float3(33.0f , 13.0f , -85.0f), 17.0f));
    objects.push_back(new Sphere(   make, float3(55.0f , 2.0f  , -90.0f), 8.0f ));
    objects.push_back(new Sphere(  ceres, float3(60.0f , -9.0f , -95.0f), 5.0f ));
    objects.push_back(new Sphere(   eris, float3(5.0f  , -13.0f, -31.0f), 1.5f ));

    Material saturnRing(GLASS, hex_to_rgb(0xAF9F88), 10.0f , {1.0f, 0.0f, 0.0f, 0.5f}, 1.0f);
    float3 ringNormal(-2.0f, -2.75f, 0.8f);
    objects.push_back(new Ring(saturnRing, ringNormal, float3(-5.0f, 8.0f, -62.0f), 21.0f, 20.25f));
    objects.push_back(new Ring(saturnRing, ringNormal, float3(-5.0f, 8.0f, -62.0f), 25.0f, 21.25f));
    objects.push_back(new Ring(saturnRing, ringNormal, float3(-5.0f, 8.0f, -62.0f), 28.0f, 25.25f));

    render(objects, lights);
}


void scene_2() {
    std::cout << "scene 2" << std::endl;
    // environment setting
    config.width = 1024;
    config.height = 720;
    load_cubemap({"../env/cubmap/cloudy/px.jpg",
                  "../env/cubmap/cloudy/nx.jpg",
                  "../env/cubmap/cloudy/py.jpg",
                  "../env/cubmap/cloudy/ny.jpg",
                  "../env/cubmap/cloudy/pz.jpg",
                  "../env/cubmap/cloudy/nz.jpg"});
    config.isCubemap    = true;
    config.antiAliasing = false;
    if (config.outFile.empty()) config.outFile = "../out_2_scene.jpg";

    std::vector<Object *> objects;
    std::vector<Light>    lights;

    lights.emplace_back(float3(0.0f  , 40.0f, -5.0f ), 0.4f);
    lights.emplace_back(float3(10.0f , 30.0f, 10.0f ), 0.1f);
    lights.emplace_back(float3(10.0f , 30.0f, -20.0f), 0.1f);
    lights.emplace_back(float3(-10.0f, 30.0f, 10.0f ), 0.1f);
    lights.emplace_back(float3(-10.0f, 30.0f, -20.0f), 0.1f);

    Material mirror(MIRROR , float3(1.0f, 1.0f, 1.0f), 100.0f, {0.5f, 1.0f, 1.0f});
    Material  glass(GLASS  , float3(1.0f, 1.0f, 1.0f), 10.0f , {0.3f, 0.5f, 0.3f, 0.7f}, 1.05f);
    Material  floor(MIRROR , float3(0.3f, 0.3f, 1.0f), 0.0f  , {0.2f, 0.0f, 1.0f});

    objects.push_back(new Plane(floor, float3(0.0f, 1.0f, 0.0f), float3(0.0f, -1.0f, 0.0f)));

    objects.push_back(new Sphere(mirror, float3(1.0f , 0.0f , -3.25f),  1.0f ));
    objects.push_back(new Sphere(glass , float3(0.3f , 0.55f, -7.0f ),  1.55f));
    objects.push_back(new Sphere(mirror, float3(-2.0f, 1.5f , -15.0f),  2.5f ));

    render(objects, lights);
}


void scene_3() {
    std::cout << "scene 3" << std::endl;
    // environment setting
    config.width  = 512 * 2;
    config.height = 512 * 2;
//    load_cubemap({"../env/cubmap/SantaMariaDeiMiracoli/px.jpg",
//                  "../env/cubmap/SantaMariaDeiMiracoli/nx.jpg",
//                  "../env/cubmap/SantaMariaDeiMiracoli/py.jpg",
//                  "../env/cubmap/SantaMariaDeiMiracoli/ny.jpg",
//                  "../env/cubmap/SantaMariaDeiMiracoli/pz.jpg",
//                  "../env/cubmap/SantaMariaDeiMiracoli/nz.jpg"});
//    config.isCubemap    = true;
    config.antiAliasing = false;
    config.backgroundColor = hex_to_rgb(0xEECF8F);
    if (config.outFile.empty()) config.outFile = "../out_3_scene.jpg";

    std::vector<Object *> objects;
    std::vector<Light>    lights;

    lights.emplace_back(float3(-10.0f, 0.0f,  10.0f), 0.5f);
    lights.emplace_back(float3(10.0f , 0.0f,  10.0f), 0.5f);

    Material m1(MIRROR, hex_to_rgb(0x01FFFD), 10.0f, {0.8f, 0.5f, 0.2f});
    Material m2(MIRROR, hex_to_rgb(0xDA651F), 10.0f, {0.8f, 0.5f, 0.2f});
    Material s1(FONG, hex_to_rgb(0xFFCD69), 10.0f, {1.0f, 0.5f});
    Material s2(MIRROR, float3(1.0f, 1.0f, 1.0f), 0.0f  , {0.2f, 0.0f, 1.0f});
    Material floor(MIRROR, float3(1.0f, 1.0f, 1.0f), 0.0f  , {0.1f, 0.0f, 0.9f});

//    objects.push_back(new Plane(floor, float3(0.0f, 1.0f, 0.0f), float3(0.0f, -1.0f, 0.0f)));

    float3 p1(-1.25f, -0.92f, -1.75f);
    float3 p2(1.25f, -0.92f, -1.75f);
    float3 p3(1.25f, -0.92f, -3.25f);
    float3 p4(-1.25f, -0.92f, -3.25f);
    float3 sc(0.0f, 1.6f + 0.15f, -2.25f);

//    objects.push_back(new Sphere(s1, float3(0.0f, 1.6f, -2.5f), 0.5));
    objects.push_back(new Sphere(s2, float3(0.0f, 0.0f, -2.5f), 0.75));

//    objects.push_back(new Cylinder(m2, p1, p2, 0.08f));
//    objects.push_back(new Cylinder(m2, p2, p3, 0.08f));
//    objects.push_back(new Cylinder(m2, p3, p4, 0.08f));
//    objects.push_back(new Cylinder(m2, p4, p1, 0.08f));
//
//    objects.push_back(new Cylinder(m1, p1, sc, 0.08f));
//    objects.push_back(new Cylinder(m1, p2, sc, 0.08f));
//    objects.push_back(new Cylinder(m1, p3, sc, 0.08f));
//    objects.push_back(new Cylinder(m1, p4, sc, 0.08f));

    render(objects, lights);
}


void scene_4() {
    std::cout << "scene 4" << std::endl;
    // environment setting
    config.width  = 512;
    config.height = 512;
    config.antiAliasing = false;
    config.backgroundColor = hex_to_rgb(0xFFA6DC);
    if (config.outFile.empty()) config.outFile = "../out_4_scene.jpg";

    std::vector<Object *> objects;
    std::vector<Light>    lights;

    Material    red(FONG, float3(1.0f, 0.0f, 0.0f), 10.0f  , {1.0f, 0.5f});
    Material  green(FONG, hex_to_rgb(0x105D00)    , 10.0f  , {1.0f, 0.5f});
    Material normal(NORMAL, float3(), 0.0f, {});

    lights.emplace_back(float3(0.0f, 0.0f, 10.0f), 1.0f);
    objects.push_back(new Mesh(normal, "../obj/bun_small.obj", 3, float3(0, 0, -5)));
//    objects.push_back(new Sphere(red, float3(0, 0, -5), 3./2.));

    render(objects, lights);
}


void parsing_parameters(int argc, char **argv, int &nScene,
        std::unordered_map<std::string, std::string> &cmdLineParams) {
    for(uint i = 0; i < argc; ++i) {
        std::string key(argv[i]);
        if(!key.empty() && (key[0] == '-')) {
            if (i != argc - 1) {
                cmdLineParams[key] = argv[i + 1];
                i++;
            } else {
                cmdLineParams[key] = "";
            }
        }
    }

    if (cmdLineParams.find("-out") != cmdLineParams.end()) {
        config.outFile = cmdLineParams["-out"];
    }
    if (cmdLineParams.find("-threads") != cmdLineParams.end()) {
        int nThreads = atoi(cmdLineParams["-threads"].c_str());
        if (nThreads <= 0) {
            std::cerr << "Warning: Invalid <threads> value: " << nThreads << std::endl;
            std::cerr << "         Parallelization is disabled" << std::endl;
        } else {
            omp_set_num_threads(nThreads);
        }
    }
    if (cmdLineParams.find("-scene") != cmdLineParams.end()) {
        nScene = atoi(cmdLineParams["-scene"].c_str());
    }
}


int main(int argc, char **argv) {
    std::unordered_map<std::string, std::string> cmdLineParams;
    int nScene = 4;

    parsing_parameters(argc, argv, nScene, cmdLineParams);

    std::cout << "Run" << std::endl;
    clock_t start_time = clock();
    switch (nScene) {
        case 1: scene_1(); break;
        case 2: scene_2(); break;
        case 3: scene_3(); break;
        case 4: scene_4(); break;
        default: std::cerr << "Warning: Invalid <scene> value: " << nScene << std::endl;
    }
    std::cout << "End, time: " << (float)(clock() - start_time) / CLOCKS_PER_SEC << " sec" << std::endl;
    return 0;
}
