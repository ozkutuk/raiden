#include "tinymath.h"

#include <cmath>

namespace tmath {

vec3f::vec3f(float x, float y, float z) : x(x), y(y), z(z) {
}

vec3f &vec3f::operator+=(const vec3f &rhs) {
    x += rhs.x;
    y += rhs.y;
    z += rhs.z;

    return *this;
}

bool operator==(const vec3f &lhs, const vec3f &rhs) {
    return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
}

bool operator!=(const vec3f &lhs, const vec3f &rhs) {
    return !(lhs == rhs);
}

vec3f operator+(vec3f lhs, const vec3f &rhs) {
    lhs += rhs;
    return lhs;
}

vec3f &vec3f::operator-=(const vec3f &rhs) {
    x -= rhs.x;
    y -= rhs.y;
    z -= rhs.z;

    return *this;
}

vec3f operator-(vec3f lhs, const vec3f &rhs) {
    lhs -= rhs;
    return lhs;
}

vec3f &vec3f::operator*=(float scalar) {
    x *= scalar;
    y *= scalar;
    z *= scalar;

    return *this;
}

vec3f operator*(vec3f lhs, float scalar) {
    lhs *= scalar;
    return lhs;
}

vec3f operator*(float scalar, vec3f rhs) {
    rhs *= scalar;
    return rhs;
}
vec3f &vec3f::operator/=(float scalar) {
    x /= scalar;
    y /= scalar;
    z /= scalar;

    return *this;
}

vec3f operator/(vec3f lhs, float scalar) {
    lhs /= scalar;
    return lhs;
}

vec3f &vec3f::operator*=(const vec3f &rhs) {
    x *= rhs.x;
    y *= rhs.y;
    z *= rhs.z;

    return *this;
}

vec3f operator*(vec3f lhs, const vec3f &rhs) {
    lhs *= rhs;
    return lhs;
}

float length(const vec3f &vec) {
    float result = std::sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
    return result;
}

vec3f normalize(vec3f vec) {
    float vecLength = length(vec);
    if (vecLength == 0.0f)
        return vec3f(); // return zero vector

    vec /= vecLength;
    return vec;
}

float dot(const vec3f &lhs, const vec3f &rhs) {
    return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

vec3f cross(const vec3f &lhs, const vec3f &rhs) {
    vec3f result;
    result.x = lhs.y * rhs.z - lhs.z * rhs.y;
    result.y = lhs.z * rhs.x - lhs.x * rhs.z;
    result.z = lhs.x * rhs.y - lhs.y * rhs.x;
    return result;
}

vec3f reflect(const vec3f &incident, const vec3f &normal) {
    return incident - 2.0f * normal * dot(incident, normal);
}

vec3f clamp(vec3f v, float min_value, float max_value) {
    v.x = std::min(std::max(v.x, min_value), max_value);
    v.y = std::min(std::max(v.y, min_value), max_value);
    v.z = std::min(std::max(v.z, min_value), max_value);
    return v;
}

} // namespace tmath
