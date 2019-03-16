#include "material.h"

#include <utility>

#include "tinymath.h"

Material::Material(tmath::vec3f diffuse, tmath::vec3f specular, tmath::vec3f ambient,
                   tmath::vec3f mirror_reflectance, float specular_exponent, float refractive_index,
                   tmath::vec3f transparency)
    : diffuse(std::move(diffuse)), specular(std::move(specular)), ambient(std::move(ambient)),
      mirror_reflectance(std::move(mirror_reflectance)), specular_exponent(specular_exponent),
      refractive_index(refractive_index), transparency(std::move(transparency)) {
}
