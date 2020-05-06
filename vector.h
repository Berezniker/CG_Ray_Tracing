#ifndef RAY_TRACING_VECTOR_H
#define RAY_TRACING_VECTOR_H

#include <iostream>  // std::ostream
#include <cmath>     // M_PI, std::sin(), std::cos(), std::sqrt()
                     // std::fmin(), std::fmax(), std::abs()

struct float3 {
    /* Coordinates */
    float x, y, z;

    /* Constructor */
    explicit float3() : x(0.0f), y(0.0f), z(0.0f) {}
    explicit float3(float xyz) : x(xyz), y(xyz), z(xyz) {}
    explicit float3(float x, float y, float z) : x(x), y(y), z(z) {}

    /* Operator Overloading */
    bool    operator == (float s) const { return (x == s) && (y == s) && (z == s); }
    float3  operator +  (float s) const { return float3(x + s, y + s, z + s); }
    float3  operator -  (float s) const { return float3(x - s, y - s, z - s); }
    float3  operator *  (float s) const { return float3(x * s, y * s, z * s); }
    float3  operator /  (float s) const { return float3(x / s, y / s, z / s); }
    float3& operator += (float s) { x += s, y += s, z += s; return *this; }
    float3& operator -= (float s) { x -= s, y -= s, z -= s; return *this; }
    float3& operator *= (float s) { x *= s, y *= s, z *= s; return *this; }
    float3& operator /= (float s) { x /= s, y /= s, z /= s; return *this; }

    bool    operator == (const float3 &v) const { return (x == v.x) && (y == v.y) && (z == v.z); }
    float3  operator +  (const float3 &v) const { return float3(x + v.x, y + v.y, z + v.z); }
    float3  operator -  (const float3 &v) const { return float3(x - v.x, y - v.y, z - v.z); }
    float3& operator += (const float3 &v) { x += v.x, y += v.y, z += v.z; return *this; }
    float3& operator -= (const float3 &v) { x -= v.x, y -= v.y, z -= v.z; return *this; }
    float   operator *  (const float3 &v) const { return x * v.x + y * v.y + z * v.z; }

    friend float3 operator + (float s, const float3 &v) { return v + s; }
    friend float3 operator - (float s, const float3 &v) { return v - s; }
    friend float3 operator * (float s, const float3 &v) { return v * s; }
    friend float3 operator / (float s, const float3 &v) { return v / s; }

    float3 operator -  ()      const { return float3(-x, -y, -z); }
    float  operator [] (int i) const { return i ? ((i == 2) ? z : y) : x; }
    float& operator [] (int i)       { return i ? ((i == 2) ? z : y) : x; }
    friend std::ostream& operator << (std::ostream &os, const float3 &v)
    { return os << "( " << v.x << ", " << v.y << ", " << v.z << " )"; }

    /* Methods */
    bool   is_zero()   const { return (x == 0.0f) && (y == 0.0f) && (z == 0.0f); }
    float  norm()      const { return std::sqrt(x * x + y * y + z * z); }
    float3 normalize() const { return this->is_zero() ? *this : *this / this->norm(); }
    float3 abs()       const { return float3(std::abs(x), std::abs(y), std::abs(z)); }
    float  max()       const { return std::fmax(x, std::fmax(y, z)); }
    float  min()       const { return std::fmin(x, std::fmin(y, z)); }
    float3 cross(const float3 &v) const
    { return float3(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x); }

    float3 rotX(float angle) const { // Y -> Z
        if (angle == 0) return *this;
        angle *= M_PI / 180.0f;
        float c = std::cos(angle), s = std::sin(angle);
        return float3(x, y * c - z * s, y * s + z * c);
    }
    float3 rotY(float angle) const { // Z -> X
        if (angle == 0) return *this;
        angle *= M_PI / 180.0f;
        float c = std::cos(angle), s = std::sin(angle);
        return float3(x * c + z * s, y, - x * s + z * c);
    }
    float3 rotZ(float angle) const { // X -> Y
        if (angle == 0) return *this;
        angle *= M_PI / 180.0f;
        float c = std::cos(angle), s = std::sin(angle);
        return float3(x * c - y * s, x * s + y * c, z);
    }
};


#endif // RAY_TRACING_VECTOR_H
