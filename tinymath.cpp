#include "tinymath.h"

#include <algorithm>
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

#if 0
vec3f refract(const vec3f &incident, const vec3f &normal, float n, float nt) {
    float d_dot_n = tmath::dot(incident, normal);
    float n_nt = n / nt;
    return n_nt * (incident - (normal * d_dot_n)) -
           normal * std::sqrt(1.0f - (n_nt * n_nt) * (1.0f - d_dot_n * d_dot_n));
}
#endif
#if 0
Vec3f refract(const vec3f &I, const vec3f &N, float ior) {
    float cosi = (dot(I, N);
    float etai = 1, etat = ior;
    vec3f n = N;
    if (cosi < 0) {
        cosi = -cosi;
    } else {
        std::swap(etai, etat);
        n = -1.0f * N;
    }
    float eta = etai / etat;
    float k = 1 - eta * eta * (1 - cosi * cosi);
    return k < 0 ? 0 : eta * I + (eta * cosi - sqrtf(k)) * n;
}
#endif

vec3f refract2(const vec3f &I, const vec3f &N, const float &refractive_index) { // Snell's law
    float cosi = -std::max(-1.f, std::min(1.f, dot(I, N)));
    float etai = 1, etat = refractive_index;
    vec3f n = N;
    if (cosi < 0) { // if the ray is inside the object, swap the indices and invert the normal to
                    // get the correct result
        cosi = -cosi;
        std::swap(etai, etat);
        n = -1.0f * N;
    }
    float eta = etai / etat;
    float k = 1 - eta * eta * (1 - cosi * cosi);
    return k < 0 ? vec3f(0, 0, 0) : I * eta + n * (eta * cosi - sqrtf(k));
}

} // namespace tmath
