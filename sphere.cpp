#include "sphere.h"
#include "ray.h"
#include "material.h"
#include "tinymath.h"

#include <cmath>
#include <utility>

Sphere::Sphere(tmath::vec3f center, float radius, Material material)
    : center(std::move(center)), radius(radius), material(std::move(material)) {
}

float Sphere::hit(const Ray &ray) const {
    float a = tmath::dot(ray.direction, ray.direction);
    tmath::vec3f o_minus_c = ray.origin - this->center;
    float b = 2 * tmath::dot(o_minus_c, ray.direction);
    float c = tmath::dot(o_minus_c, o_minus_c) - this->radius * this->radius;

    float discriminant = b * b - 4 * a * c;
    if (discriminant > 0) {
        float b2a = -b / (2 * a);
        float sqr2a = std::sqrt(discriminant) / (2 * a);
        float t = b2a - sqr2a;
        if (t < 0)
            t = b2a + sqr2a;
        // we do not handle both of them being negative, negative
        // return values are handled as not-intersected anyway
        return t;
    } else
        return -1.0f; // NOTE: this could be handled better
}
