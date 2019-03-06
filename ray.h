#pragma once

#include "tinymath.h"

class Ray {
  public:
    explicit Ray(tmath::vec3f origin, tmath::vec3f direction);
    Ray(const Ray &) = default;
    Ray &operator=(const Ray &) = default;
    tmath::vec3f at(float t) const;

    tmath::vec3f origin;
    tmath::vec3f direction;
};
