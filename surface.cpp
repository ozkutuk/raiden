#include "surface.h"
#include "material.h"
#include "ray.h"
#include "tinymath.h"

#include <cmath>
#include <utility>
#include <vector>

uint64_t Triangle::test_count = 0;
uint64_t Triangle::hit_count = 0;

Box::Box(tmath::vec3f min_point, tmath::vec3f max_point)
    : min_point(std::move(min_point)), max_point(std::move(max_point)) {
}

bool Box::hit(const Ray &ray) const {
    float tmin, tmax, tymin, tymax, tzmin, tzmax;
    if (ray.direction.x >= 0) {
        tmin = (min_point.x - ray.origin.x) / ray.direction.x;
        tmax = (max_point.x - ray.origin.x) / ray.direction.x;
    } else {
        tmin = (max_point.x - ray.origin.x) / ray.direction.x;
        tmax = (min_point.x - ray.origin.x) / ray.direction.x;
    }
    if (ray.direction.y >= 0) {
        tymin = (min_point.y - ray.origin.y) / ray.direction.y;
        tymax = (max_point.y - ray.origin.y) / ray.direction.y;
    } else {
        tymin = (max_point.y - ray.origin.y) / ray.direction.y;
        tymax = (min_point.y - ray.origin.y) / ray.direction.y;
    }
    if ((tmin > tymax) || (tymin > tmax))
        return false;
    if (tymin > tmin)
        tmin = tymin;
    if (tymax < tmax)
        tmax = tymax;
    if (ray.direction.z >= 0) {
        tzmin = (min_point.z - ray.origin.z) / ray.direction.z;
        tzmax = (max_point.z - ray.origin.z) / ray.direction.z;
    } else {
        tzmin = (max_point.z - ray.origin.z) / ray.direction.z;
        tzmax = (min_point.z - ray.origin.z) / ray.direction.z;
    }
    if ((tmin > tzmax) || (tzmin > tmax))
        return false;
    if (tzmin > tmin)
        tmin = tzmin;
    if (tzmax < tmax)
        tmax = tzmax;
    return ((tmin < T_MAX) && (tmax > T_MIN));
}

void Box::update(const tmath::vec3f &point) {
    if (point.x < min_point.x)
        min_point.x = point.x;
    if (point.y < min_point.y)
        min_point.y = point.y;
    if (point.z < min_point.z)
        min_point.z = point.z;

    if (point.x > max_point.x)
        max_point.x = point.x;
    if (point.y > max_point.y)
        max_point.y = point.y;
    if (point.z > max_point.z)
        max_point.z = point.z;
}

SurfaceList::SurfaceList(std::vector<std::unique_ptr<Surface>> surfaces)
    : surfaces(std::move(surfaces)) {
}

std::optional<HitRecord> SurfaceList::hit(const Ray &ray) const {
    float t = std::numeric_limits<float>::max();
    std::optional<HitRecord> hit_record;
    for (const auto &surface : surfaces) {
        auto intersected = surface->hit(ray);
        if (intersected && intersected->t < t) {
            hit_record = intersected;
            t = intersected->t;
        }
    }

    return hit_record;
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

Triangle::Triangle(tmath::vec3f v1, tmath::vec3f v2, tmath::vec3f v3, Material material)
    : v1(std::move(v1)), v2(std::move(v2)), v3(std::move(v3)), material(std::move(material)) {
}

std::optional<HitRecord> Triangle::hit(const Ray &ray) const {
    // [ a d g ][ beta  ]   [ j ]
    // [ b e h ][ gamma ] = [ k ]
    // [ c f i ][   t   ]   [ l ]

    test_count += 1;

    float a = v1.x - v2.x;
    float b = v1.y - v2.y;
    float c = v1.z - v2.z;

    float d = v1.x - v3.x;
    float e = v1.y - v3.y;
    float f = v1.z - v3.z;

    float g = ray.direction.x;
    float h = ray.direction.y;
    float i = ray.direction.z;

    float j = v1.x - ray.origin.x;
    float k = v1.y - ray.origin.y;
    float l = v1.z - ray.origin.z;

    float ei_hf = e * i - h * f;
    float gf_di = g * f - d * i;
    float dh_eg = d * h - e * g;
    float denominator = a * ei_hf + b * gf_di + c * dh_eg;

    float beta = (j * ei_hf + k * gf_di + l * dh_eg) / denominator;
    if (beta < 0.0f || beta > 1.0f)
        return {};

    float ak_jb = a * k - j * b;
    float jc_al = j * c - a * l;
    float bl_kc = b * l - k * c;
    float gamma = (i * ak_jb + h * jc_al + g * bl_kc) / denominator;
    if (gamma < 0.0f || gamma > 1.0f)
        return {};

    if (beta + gamma > 1.0f)
        return {};

    float t = -(f * ak_jb + e * jc_al + d * bl_kc) / denominator;
    if (t > T_MIN && t < T_MAX) {
        HitRecord hit_record;
        hit_record.t = t;
        hit_record.normal = tmath::normalize(tmath::cross(v2 - v1, v3 - v1));
        hit_record.material = &(this->material);
        hit_count += 1;
        return hit_record;
    }

    return {};
}

Mesh::Mesh(std::vector<Triangle> triangles, Material material)
    : triangles(std::move(triangles)), material(std::move(material)) {

    for (const auto &triangle : this->triangles) {
        bounding_box.update(triangle.v1);
        bounding_box.update(triangle.v2);
        bounding_box.update(triangle.v3);
    }
}

std::optional<HitRecord> Mesh::hit(const Ray &ray) const {
    float t = std::numeric_limits<float>::max();
    std::optional<HitRecord> hit_record;
    if (bounding_box.hit(ray)) {
        for (const auto &triangle : triangles) {
            auto intersected = triangle.hit(ray);
            if (intersected && intersected->t < t) {
                hit_record = intersected;
                t = intersected->t;
            }
        }
    }

    return hit_record;
}

