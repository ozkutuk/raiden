#pragma once

#include "tinymath.h"

class Material {
  public:
    explicit Material(tmath::vec3f diffuse, tmath::vec3f specular, tmath::vec3f ambient,
                      tmath::vec3f mirror_reflectance, float specular_exponent, float refractive_index, tmath::vec3f transparency);

    tmath::vec3f diffuse;
    tmath::vec3f specular;
    tmath::vec3f ambient;
    tmath::vec3f mirror_reflectance;
    float specular_exponent;
    float refractive_index;
    tmath::vec3f transparency;
};

