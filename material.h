#pragma once

#include "tinymath.h"

// TODO: add specular component
class Material {
  public:
    explicit Material(tmath::vec3f diffuse, tmath::vec3f specular, tmath::vec3f ambient,
                      tmath::vec3f mirror_reflectance, float specular_exponent);

    tmath::vec3f diffuse;
    tmath::vec3f specular;
    tmath::vec3f ambient;
    tmath::vec3f mirror_reflectance;
    float specular_exponent;
};

