#include <unordered_map>  // std::unordered_map
#include <iostream>       // std::cout, std::cerr, std::endl
#include <cstdlib>        // atoi()
#include <vector>         // std::vector
#include <string>         // std::string
#include <omp.h>          // omp_set_num_threads(), omp_get_wtime()
#include "objects.h"
#include "vector.h"
#include "tracer.h"
#include "image.h"
#include "mesh.h"


/*    Planet    */
void scene_1() {
    std::cout << "scene 1" << std::endl;
    // environment setting
    config.width = 1280;
    config.height = 640;
    load_background_image("../env/planet/stars_milky_way.jpg");
    config.isBackgroundMap = true;
    config.antiAliasing    = true;
    if (config.outFile.empty()) config.outFile = "../out/out_1_scene.jpg";

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


/*    Mirror    */
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
    config.toneMapping  = true;
    config.maxDepth     = 5;
    if (config.outFile.empty()) config.outFile = "../out/out_2_scene.jpg";

    std::vector<Object *> objects;
    std::vector<Light>    lights;

    lights.emplace_back(float3(0.0f  , 40.0f, -5.0f ), 0.4f);
    lights.emplace_back(float3(10.0f , 30.0f, 10.0f ), 0.1f);
    lights.emplace_back(float3(10.0f , 30.0f, -20.0f), 0.1f);
    lights.emplace_back(float3(-10.0f, 30.0f, 10.0f ), 0.1f);
    lights.emplace_back(float3(-10.0f, 30.0f, -20.0f), 0.1f);

    Material mirror(MIRROR, float3(1.0f, 1.0f, 1.0f), 100.0f, {0.5f, 1.0f, 1.0f});
    Material  glass(GLASS , float3(1.0f, 1.0f, 1.0f), 0.0f  , {1.0f, 0.0f, 1.0f, 0.5f}, 1.0f);
    Material  floor(MIRROR, float3(0.3f, 0.3f, 1.0f), 0.0f  , {0.2f, 0.0f, 1.0f});

    objects.push_back(new Plane(floor, float3(0.0f, 1.0f, 0.0f), float3(0.0f, -1.0f, 0.0f)));

    objects.push_back(new Hyperboloid(glass, 1.0f, 1.5f, 1.0f,1.0f,9.5f, float3(0.0f, -1.0f, -10.0f)));

    float extRadius = 3.5f;
    for (float phi = 0.45; phi <= 2.0f * M_PI; phi += M_PI / 4.0f) {
        float3 center(extRadius * std::sin(phi), 0.5f, -10 + extRadius * std::cos(phi));
        objects.push_back(new Sphere(mirror, center, 1.0f));
    }
    for (float y = -1.0f; y <= 5.0f; y += 1.0f) {
        objects.push_back(new Sphere(mirror, float3(0.0f, y, -10.0f), 0.25f));
    }

    render(objects, lights);
}


/*     Mesh     */
void scene_3() {
    std::cout << "scene 3" << std::endl;
    // environment setting
    config.width  = 512;
    config.height = 512;
    config.backgroundColor = hex_to_rgb(0x3C9AFF);
    if (config.outFile.empty()) config.outFile = "../out/out_3_scene.jpg";

    std::vector<Object *> objects;
    std::vector<Light>    lights;

    Material normal(NORMAL, float3(), 0.0f, {});

    objects.push_back(new Mesh(normal, "../obj/deer.obj", 2.0f, float3(0.8f, 0.1f, -5.0f),
                               float3(0.0f, -25.0f, 0.0f)));
    objects.push_back(new Mesh(normal, "../obj/deer2.obj", 2.0f, float3(-0.2f, 0.1f, -5.0f),
                               float3(0.0f, -110.0f, 0.0f)));

    render(objects, lights);
}


/* Ray Marching */
void scene_4() {
    std::cout << "scene 4" << std::endl;
    // environment setting
    config.width  = 512;
    config.height = 512;
    config.antiAliasing = true;
    if (config.outFile.empty()) config.outFile = "../out/out_4_scene.jpg";

    std::vector<Object *> objects;
    std::vector<Light>    lights;

    lights.emplace_back(float3(0.0f, 10.0f, 10.0f), 1.0f);

    Material shine(FONG, float3(1.0f, 1.0f, 1.0f), 100.0f, {1.0f, 1.0f});

    objects.push_back(new Fractal(shine, float3(0, 0, -5), 1.5f));

    render(objects, lights);
}


void parsing_parameters(int argc, char **argv, int &nScene,
        std::unordered_map<std::string, std::string> &cmdLineParams) {
    for (uint i = 0; i < argc; ++i) {
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
            std::cerr << "[Warning]: Invalid <threads> value: " << nThreads << std::endl
                      << "           Parallelization is disabled" << std::endl;
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
    int nScene = 1;

    parsing_parameters(argc, argv, nScene, cmdLineParams);

    std::cout << "__Run__" << std::endl;
    double start_time = omp_get_wtime();
    switch (nScene) {
        case 1: scene_1(); break;
        case 2: scene_2(); break;
        case 3: scene_3(); break;
        case 4: scene_4(); break;
        default: std::cerr << "[Warning]: Invalid <scene> value: " << nScene << std::endl
                           << "           Please specify from 1 to 4" << std::endl;
    }
    std::cout << "End, time: " << omp_get_wtime() - start_time << " sec" << std::endl;
    return 0;
}
