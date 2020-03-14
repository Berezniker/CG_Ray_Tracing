#ifndef RAY_TRACING_MESH_H
#define RAY_TRACING_MESH_H


#include <fstream>    // std::ifstream
#include <sstream>    // std::stringstream
#include "objects.h"  // Object, Triangle, MISS


struct Mesh : Object {
    std::vector<Triangle> triangles;
    Sphere boundingSphere;
    
    Mesh(Material m, const std::string &path, float size, const float3 &position) :
        Object(std::move(m)), boundingSphere(m, position, size) {
        // read .obj
        std::ifstream ifs(path, std::ios::in);
        std::vector<float3> vertices;
        std::vector<float3> faces;
        for (std::string buffer; ifs >> buffer; ) {
            float x, y, z;
            if (buffer == "v") {
                ifs >> x >> y >> z;
                vertices.emplace_back(x, y, z);
            } else if (buffer == "f") {
                ifs >> x >> y >> z;
                faces.emplace_back(x - 1, y - 1, z - 1);
            }
        }
        ifs.close();

        float3 center(0.0f);
        float3 minDist = vertices.front(), maxDist = vertices.front();
        for (auto &vertex : vertices) {
            center += vertex;
            for (uint i = 0; i < 3; ++i) {
                if (vertex[i] < minDist[i]) {
                    minDist[i] = vertex[i];
                }
                if (vertex[i] > maxDist[i]) {
                    maxDist[i] = vertex[i];
                }
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

    bool is_ray_intersect(const Ray &ray, float &t, float3 &P, float3 &N, float3 &color) const override {
        float3 tmpP, tmpN, tmpC;
        if (boundingSphere.is_ray_intersect(ray, t, tmpP, tmpN, tmpC)) {
            float closestDist = std::numeric_limits<float>::max();
            for (const auto &triangle : triangles) {
                float curDist;
                if (triangle.is_ray_intersect(ray, curDist, tmpP, tmpN, tmpC) && curDist < closestDist) {
                    closestDist = curDist;
                    P = tmpP;
                    N = tmpN;
                    color = tmpC;
                }
            }
            t = closestDist;
            return t >= 0;
        }
        return MISS;
    }
};


#endif // RAY_TRACING_MESH_H
