#include "surface.h"
#include "material.h"
#include "ray.h"
#include "tinymath.h"

#include <cmath>
#include <vector>
#include <utility>

SurfaceList::SurfaceList(std::vector<std::unique_ptr<Surface>> surfaces) : surfaces(std::move(surfaces)) {
}

std::optional<HitRecord> SurfaceList::hit(const Ray &ray) const {
    // float distance = std::numeric_limits<float>::max();
    // const Surface *closest = nullptr;
    std::optional<HitRecord> hit_record;
    for (const auto &surface : surfaces) {
        auto intersected = surface->hit(ray);
        // if (t > T_MIN && t < distance) {
        if (intersected) {
            hit_record = intersected;
            // distance = t;
            // closest = &surface;
        }
    }

    return hit_record;

#if 0
    HitRecord hit;
    if (distance < T_MAX)
        hit.sphere = closest;
    else
        hit.sphere = nullptr;
    hit.t = distance;
    return hit;
#endif

}

Sphere::Sphere(tmath::vec3f center, float radius, Material material)
    : center(std::move(center)), radius(radius), material(std::move(material)) {
}

std::optional<HitRecord> Sphere::hit(const Ray &ray) const {
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
        if (t > T_MIN && t < T_MAX) {
            HitRecord hit_record;
            hit_record.t = t;
            hit_record.normal = tmath::normalize(ray.at(t) - this->center);
            hit_record.material = &(this->material);
            return hit_record;
        }
    }

    return {};
}
