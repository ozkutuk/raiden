#pragma once

namespace tmath {

class vec3f {

  public:
    float x;
    float y;
    float z;

    explicit vec3f(float x = 0, float y = 0, float z = 0);
    vec3f(const vec3f &) = default;
    vec3f &operator=(const vec3f &) = default;

    // equality
    friend bool operator==(const vec3f &lhs, const vec3f &rhs);
    friend bool operator!=(const vec3f &lhs, const vec3f &rhs);

    // addition
    vec3f &operator+=(const vec3f &rhs);
    friend vec3f operator+(vec3f lhs, const vec3f &rhs);

    // substraction
    vec3f &operator-=(const vec3f &rhs);
    friend vec3f operator-(vec3f lhs, const vec3f &rhs);

    // scalar multiplication
    vec3f &operator*=(float scalar);
    friend vec3f operator*(vec3f lhs, float scalar);
    friend vec3f operator*(float scalar, vec3f rhs);

    vec3f &operator/=(float scalar);
    friend vec3f operator/(vec3f lhs, float scalar);

    vec3f &operator*=(const vec3f &rhs);
    friend vec3f operator*(vec3f lhs, const vec3f &rhs);
};

float length(const vec3f &vec);
vec3f normalize(vec3f vec);
float dot(const vec3f &lhs, const vec3f &rhs);
vec3f cross(const vec3f &lhs, const vec3f &rhs);

} // namespace tmath
