#pragma once

#include "tinymath.h"

// TODO: add specular component
class Material {
  public:
    explicit Material(tmath::vec3f diffuse, tmath::vec3f specular, float specular_exponent);

    tmath::vec3f diffuse;
    tmath::vec3f specular;
    float specular_exponent;
};

