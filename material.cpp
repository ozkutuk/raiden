#include "material.h"

#include <utility>

#include "tinymath.h"

Material::Material(tmath::vec3f diffuse, tmath::vec3f specular, float specular_exponent)
    : diffuse(std::move(diffuse)), specular(std::move(specular)),
      specular_exponent(specular_exponent) {
}
