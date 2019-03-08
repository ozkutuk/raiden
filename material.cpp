#include "material.h"

#include <utility>

#include "tinymath.h"

Material::Material(tmath::vec3f diffuse, tmath::vec3f specular, tmath::vec3f ambient,
                   float specular_exponent)
    : diffuse(std::move(diffuse)), specular(std::move(specular)), ambient(std::move(ambient)),
      specular_exponent(specular_exponent) {
}
