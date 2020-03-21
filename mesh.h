#ifndef RAY_TRACING_MESH_H
#define RAY_TRACING_MESH_H

#include <algorithm>  // std::count()
#include <fstream>    // std::ifstream, is_open(), close()
#include <sstream>    // std::stringstream
#include <vector>     // std::vector
#include <string>     // std::string
#include <cstdlib>    // exit()
#include "objects.h"  // Object, Sphere, Triangle, MISS
#include "vector.h"   // float3


struct Mesh : Object {
    std::vector<Triangle> triangles;
    Sphere boundingSphere;
    
    Mesh(Material m, const std::string &path, float size, const float3 &position) :
            Object(std::move(m)), boundingSphere(m, position, size) {
        std::vector<float3> vertices, faces;
        obj_parser(path, vertices, faces);

        float3 center(0.0f, 0.0f, 0.0f);
        float3 minDist = vertices.front(), maxDist = vertices.front();
        for (auto &vertex : vertices) {
            center += vertex;
            for (uint k = 0; k < 3; ++k) {
                minDist[k] = std::fmin(minDist[k], vertex[k]);
                maxDist[k] = std::fmax(maxDist[k], vertex[k]);
            }
        }
        center /= vertices.size();

        float scale = size / (maxDist - minDist).max();
        for (auto &vertex : vertices) {
            vertex = (vertex - center) * scale + position;
        }

        for (auto &face : faces) {
            triangles.emplace_back(Triangle(m, vertices[face[0]], vertices[face[1]], vertices[face[2]]));
        }
    }

    static void obj_parser(const std::string &path,
                           std::vector<float3> &vertices,
                           std::vector<float3> &faces) {
        std::ifstream ifs(path);
        if (!ifs.is_open()) {
            std::cerr << "[Error]: can not open file: \"" << path  << "\"" << std::endl;
            exit(-1);
        }
        std::string buffer;
        for (ifs >> buffer; !ifs.eof(); ifs >> buffer) {
            if (buffer == "v") {
                // template: v x y z
                float x, y, z;
                ifs >> x >> y >> z;
                vertices.emplace_back(x, y, z);
            } else if (buffer == "f") {
                int x, y, z, t;
                std::string line;
                getline(ifs, line);
                if (line.find('/') != std::string::npos &&
                    std::count(line.begin(), line.end(), '/') == 6) {
                    // template: f v1/vt1/vn1 v2/vt2/vn2 v3/vt3/vn3
                    sscanf(line.c_str(), " %d/%d/%d %d/%d/%d %d/%d/%d\n", &x, &t, &t, &y, &t, &t, &z, &t, &t);
                } else if (line.find('/') != std::string::npos &&
                           std::count(line.begin(), line.end(), '/') == 3) {
                    // template: f v1/vn1 v2/vn2 v3/vn3
                    sscanf(line.c_str(), " %d/%d %d/%d %d/%d\n", &x, &t, &y, &t, &z, &t);
                } else if (line.find("//") != std::string::npos &&
                           std::count(line.begin(), line.end(), '/') == 6) {
                    // template: f v1//vn1 v2//vn2 v3//vn3
                    sscanf(line.c_str(), " %d//%d %d//%d %d//%d\n", &x, &t, &y, &t, &z, &t);
                } else {
                    // template: f v1 v2 v3
                    sscanf(line.c_str(), " %d %d %d\n", &x, &y, &z);
                }
                faces.emplace_back(x - 1, y - 1, z - 1);
            }
        }
        ifs.close();
    }

    bool is_ray_intersect(const Ray &ray, float &t, float3 &P, float3 &N, float3 &color) const override {
        float3 tmpP, tmpN, tmpC;
        if (boundingSphere.is_ray_intersect(ray, t, tmpP, tmpN, tmpC)) {
            float curDist, closestDist = std::numeric_limits<float>::max();
            for (const auto &triangle : triangles) {
                if (triangle.is_ray_intersect(ray, curDist, tmpP, tmpN, tmpC) && curDist < closestDist) {
                    closestDist = curDist;
                    P = tmpP, N = tmpN, color = tmpC;
                }
            }
            t = closestDist;
            return t >= 0;
        }
        return MISS;
    }
};


#endif // RAY_TRACING_MESH_H
