#pragma once

#include "material.h"
#include "tinymath.h"

class Ray;

class Sphere {
  public:
    explicit Sphere(tmath::vec3f center, float radius, Material material);
    float hit(const Ray &ray) const;

    tmath::vec3f center;
    float radius;
    Material material;
};

