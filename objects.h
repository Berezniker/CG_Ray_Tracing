#ifndef RAY_TRACING_OBJECTS_H
#define RAY_TRACING_OBJECTS_H

#include <iostream>  // std::cerr, std::endl, std::move()
#include <limits>    // std::numeric_limits<float>::max()
#include <vector>    // std::vector
#include <cmath>     // M_PI, std::fmin(), std::sin(), sqrtf(), atan2f()
#include "vector.h"  // float3, cross(), norm(), normalize()
#include "image.h"   // load_image()


struct Ray {
    float3 origin;
    float3 direction;

    Ray(const float3 &origin, const float3 &direction) :
            origin(origin), direction(direction.normalize()) {}

    float3 substitute(float t) const { return origin + t * direction; }
};


struct Light {
    float3 position;
    float intensity;

    Light(const float3 &position, float intensity) :
           position(position), intensity(intensity) {}
};


enum MATERIAL_TYPE {NORMAL, LAMBERT, FONG, MIRROR, GLASS};

struct Material {
    MATERIAL_TYPE materialType;
    float3 color;
    float shine;
    float eta;
    std::vector<float> albedo;

    std::vector<float3> texture;
    int textureWidth  = 0;
    int textureHeight = 0;


    explicit Material(MATERIAL_TYPE materialType, const float3 &color,
                      float shine, const std::vector<float> &albedo, float eta=1.0f) :
            materialType(materialType), color(color), shine(shine), eta(eta) {
        set_albedo(albedo);
    }

    explicit Material(MATERIAL_TYPE materialType, const char *texturePath,
                      float shine, const std::vector<float> &albedo, float eta=1.0f) :
            materialType(materialType), color(),  shine(shine), eta(eta) {
        set_albedo(albedo);
        load_image(texturePath, texture, textureWidth, textureHeight);
    }

    void set_albedo(const std::vector<float> &initAlbedo) {
        albedo = std::vector<float>(4, 0.0f);
        switch (materialType) {
            case GLASS  : albedo[3] = initAlbedo[3];
            case MIRROR : albedo[2] = initAlbedo[2];
            case FONG   : albedo[1] = initAlbedo[1];
            case LAMBERT: albedo[0] = initAlbedo[0];
            case NORMAL : break;
        }
    }

    float3 get_color(float x=0.0f, float y=0.0f) const {
        if (texture.empty()) {
            return color;
        } else {
            int u = (int)(x * (float)textureWidth), v = (int)(y * (float)textureHeight);
            return texture[v * textureWidth + u];
        }
    }
};


#define HIT  true
#define MISS false

struct Object {
    Material material;

    explicit Object(Material m) : material(std::move(m)) {};
    virtual ~Object() = default;
    virtual bool is_ray_intersect(const Ray &, float &, float3 &, float3 &, float3&) const = 0;
};

struct Triangle : Object {
    float3 A, B, C;
    float3 normal;

    Triangle(Material m, const float3 &A, const float3 &B, const float3 &C) :
            Object(std::move(m)), A(A), B(B), C(C) {
        normal = (B - A).cross(C - A).normalize();
        if (normal == 0.0f) {
            std::cerr << "[Warning]: the points of the triangle are lying on one line" << std::endl;
        }
    }

    bool is_ray_intersect(const Ray &ray, float &t, float3 &P, float3 &N, float3 &color) const override {
        t = (normal * ray.origin + normal * A) / (normal * ray.direction);
        if (t < 0) return MISS;
        P = ray.substitute(t);
        N = get_normal();
        color = Object::material.get_color();
        // inside-outside test:
        return normal * (B - A).cross(P - A) > 0 &&
               normal * (C - B).cross(P - B) > 0 &&
               normal * (A - C).cross(P - C) > 0;
    }

    inline float3 get_normal() const { return normal; }
};

struct Sphere : Object{
    const float3 center;
    const float radius;

    Sphere(Material m, const float3 &center, float radius) :
            Object(std::move(m)), center(center), radius(radius) {}

    bool is_ray_intersect(const Ray &ray, float &t, float3 &P, float3 &N, float3 &color) const override {
        float3 oc = ray.origin - center;
        float a = ray.direction * ray.direction;
        float b = oc * ray.direction;
        float c = oc * oc - radius * radius;
        float D = b * b - a * c;
        if (D < 0) return MISS;
        t = std::fmin((-b - sqrtf(D)) / a,(-b + sqrtf(D)) / a);
        if (t < 0) return MISS;
        P = ray.substitute(t);
        N = get_normal(P);
        color = get_color(N);
        return HIT;
    }

    inline float3 get_normal(const float3 &hit) const { return (hit - center).normalize(); }

    float3 get_color(const float3 &N) const {
        float phi = ((atan2f(N.x, N.z)) + (float)M_PI) / 2.0f / (float)M_PI;
        float theta = ((atan2f(N.y, sqrtf(N.x * N.x + N.z * N.z))) + (float)M_PI / 2.0f) / (float)M_PI;
        theta = 1.0f - theta;
        return Object::material.get_color(phi, theta);
    }
};

struct Plane : Object {
    float3 normal;
    float3 p0;

    Plane(Material m, const float3 &normal, const float3 &p) :
            Object(std::move(m)), normal(normal.normalize()), p0(p) {}

    bool is_ray_intersect(const Ray &ray, float &t, float3 &P, float3 &N, float3 &color) const override {
        t = ((p0 - ray.origin) * normal) / (normal * ray.direction + (float)1e-5);
        if (t < 0) return MISS;
        P = ray.substitute(t);
        N = get_normal();
        color = Object::material.get_color();
        return HIT;
    }

    inline virtual float3 get_normal() const { return normal; }
};

struct Disk : Plane {
    float radius;

    Disk(Material m, const float3 &normal, const float3 &p, float radius) :
        Plane(std::move(m), normal.normalize(), p), radius(radius) {}

    bool is_ray_intersect(const Ray &ray, float &t, float3 &P, float3 &N, float3 &color) const override {
        if (Plane::is_ray_intersect(ray, t, P, N, color)) {
            float3 v = P - Plane::p0;
            return sqrtf(v * v) <= radius;
        }
        return MISS;
    }
};

struct Ring : Plane {
    float extRadius;
    float intRadius;

    Ring(Material m, const float3 &normal, const float3 &p, float extRadius, float intRadius) :
            Plane(std::move(m), normal.normalize(), p), extRadius(extRadius), intRadius(intRadius) {}

    bool is_ray_intersect(const Ray &ray, float &t, float3 &P, float3 &N, float3 &color) const override {
        if (Plane::is_ray_intersect(ray, t, P, N, color)) {
            float3 v = P - Plane::p0;
            float svv = sqrtf(v * v);
            return (intRadius <= svv) && (svv <= extRadius);
        }
        return MISS;
    }
};

struct Cylinder : Object {
    float3 pTop, pBottom;
    Disk    top,  bottom;
    float radius;
    Ray axis;

    Cylinder(Material m, const float3 &pB, const float3 &pT, float radius) :
            Object(std::move(m)), pTop(pT), pBottom(pB), radius(radius),
            axis(pBottom, (pT - pB).normalize()),
            bottom(std::move(m), -(pT - pB).normalize(), pB, radius),
            top(   std::move(m),  (pT - pB).normalize(), pT, radius) {}

    bool is_ray_intersect(const Ray &ray, float &t, float3 &P, float3 &N, float3 &color) const override {
        float3 dp = ray.origin - axis.origin;
        float3 tmp_a = ray.direction - (ray.direction * axis.direction) * axis.direction;
        float3 tmp_b = dp - (dp * axis.direction) * axis.direction;
        float a = tmp_a * tmp_a;
        float b = tmp_a * tmp_b;
        float c = tmp_b * tmp_b - radius * radius;
        float D = b * b - a * c;
        if (D < 0) return MISS;

        float t0, tT, tB, *interPoint = nullptr;
        t0 = std::fmin((-b - sqrtf(D)) / a, (-b + sqrtf(D)) / a);

        float3 q = ray.substitute(t0);
        if (t0 >= 0 && axis.direction * (q - pTop) < 0 &&
            axis.direction * (q - pBottom) > 0) {
            interPoint = &t0;
        }

        float3 tmpP, tmpN, tmpC;

        if (top.is_ray_intersect(ray, tT, tmpP, tmpN, tmpC) && tT >= 0) {
            if (interPoint) {
                *interPoint = std::fmin(tT, *interPoint);
            } else {
                interPoint = &tT;
            }
        }
        if (bottom.is_ray_intersect(ray, tB, tmpP, tmpN, tmpC) && tB >= 0) {
            if (interPoint) {
                *interPoint = std::fmin(tB, *interPoint);
            } else {
                interPoint = &tB;
            }
        }

        if (interPoint) {
            t = *interPoint;
            P = ray.substitute(t);
            N = get_normal(P);
            color = Object::material.get_color();
            return t >= 0;
        } else {
            return MISS;
        }
    }

    float3 get_normal(const float3 &hit) const {
        float lenBottomToHit = (hit - pBottom).norm();
        float lenTopToHit    = (hit - pTop   ).norm();
        if (lenBottomToHit <= radius && lenTopToHit <= radius) {
            return (lenBottomToHit < lenTopToHit) ? bottom.get_normal() : top.get_normal();
        }
        if (lenBottomToHit <= radius) return bottom.get_normal();
        if (lenTopToHit    <= radius) return top.get_normal();

        float3 BotToHit = hit - pBottom;
        return (BotToHit - axis.direction * (BotToHit * axis.direction)).normalize();
    }
};

struct Hyperboloid : Object {
    float3 abc2;
    float3 position;
    float  radius;
    float  halfSize;

    Hyperboloid(Material m, float a, float b, float c, float radius, float size, float3 pos) :
        Object(std::move(m)), radius(radius), position(pos), halfSize(size / 2.0f),
        abc2(1.0f / (a * a), -1.0f / (b * b), 1.0f / (c * c)) {}

    bool is_ray_intersect(const Ray &ray, float &t, float3 &P, float3 &N, float3 &color) const override {
        float3 center = ray.origin - position;
        float a = float3(ray.direction.x * ray.direction.x,
                         ray.direction.y * ray.direction.y,
                         ray.direction.z * ray.direction.z) * abc2;
        float b = float3(center.x * ray.direction.x,
                         center.y * ray.direction.y,
                         center.z * ray.direction.z) * abc2;
        float c = float3(center.x * center.x,
                         center.y * center.y,
                         center.z * center.z) * abc2 - radius;
        float D = b * b - a * c;
        if (D < 0) return MISS;
        t = std::fmin((-b - sqrtf(D)) / a, (-b + sqrtf(D)) / a);
        if (t < 0) return MISS;
        P = ray.substitute(t);
        if (std::fabs(P.y - position.y) >= halfSize) return MISS;
        N = get_normal(P);
        color = Object::material.get_color();
        return HIT;
    }

    float3 get_normal(const float3 &hit) const {
        return float3((hit.x - position.x) * abc2.x,
                      (hit.y - position.y) * abc2.y,
                      (hit.z - position.z) * abc2.z).normalize();
    }
};

struct Fractal: Object {
    const  float3 center;
    const  float  radius;
    static float  minHitDist;
    static float  maxTraceDist;
    static int    numSteps;

    Fractal(Material m, const float3 &center, float radius) :
        Object(std::move(m)), center(center), radius(radius) {}

    float SDF(const float3 &P, float k=10.0f) const {
        float distortion = std::sin(k * P.x) * std::sin(k * P.y) * std::sin(k * P.z) * 0.05f;
        float sphereDist = (P - center).norm() - radius;
        return sphereDist + distortion;
    }

    bool is_ray_intersect(const Ray &ray, float &t, float3 &P, float3 &N, float3 &color) const override {
        float totalDist = 0.0f, currentDist;
        for (uint i = 0; (i < numSteps) && (totalDist < maxTraceDist); ++i) {
            P = ray.substitute(totalDist);
            currentDist = SDF(P);
            if (currentDist < minHitDist) {
                N = get_normal(P);
                color = 0.5f * N + 0.5f;  // loophole
            }
            totalDist += currentDist;
        }
//        color = Object::material.get_color();
        return totalDist < maxTraceDist;
    }

    float3 get_normal(const float3 &hit) const {
        const float3 dx = float3(0.001f, 0.0f, 0.0f);
        const float3 dy = float3(0.0f, 0.001f, 0.0f);
        const float3 dz = float3(0.0f, 0.0f, 0.001f);
        float gradX = SDF(hit + dx) - SDF(hit - dx);
        float gradY = SDF(hit + dy) - SDF(hit - dy);
        float gradZ = SDF(hit + dz) - SDF(hit - dz);
        return float3(gradX, gradY, gradZ).normalize();
    }
};

int   Fractal::numSteps     = 64;
float Fractal::minHitDist   = 1e-3;
float Fractal::maxTraceDist = 100.0f;


#endif // RAY_TRACING_OBJECTS_H
