#include "surface.h"
#include "material.h"
#include "ray.h"
#include "tinymath.h"

#include <cmath>
#include <utility>
#include <vector>
#include <algorithm>

uint64_t Face::test_count = 0;
uint64_t Face::hit_count = 0;
uint64_t Box::test_count = 0;
uint64_t Box::hit_count = 0;
uint64_t Sphere::test_count = 0;
uint64_t Sphere::hit_count = 0;

Box::Box()
    : min_point(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(),
                std::numeric_limits<float>::max()),
      max_point(std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest(),
                std::numeric_limits<float>::lowest()) {
}

Box::Box(tmath::vec3f min_point, tmath::vec3f max_point)
    : min_point(std::move(min_point)), max_point(std::move(max_point)) {
}

Box::Box(const Box &left, const Box &right) : Box() {
    update(left.min_point);
    update(left.max_point);
    update(right.min_point);
    update(right.max_point);
}

bool Box::hit(const Ray &ray) const {
    float tmin, tmax, tymin, tymax, tzmin, tzmax, divx, divy, divz;

    test_count += 1;
    divx = 1 / ray.direction.x;
    if (divx >= 0) {
        tmin = (min_point.x - ray.origin.x) * divx;
        tmax = (max_point.x - ray.origin.x) * divx;
    } else {
        tmin = (max_point.x - ray.origin.x) * divx;
        tmax = (min_point.x - ray.origin.x) * divx;
    }

    divy = 1 / ray.direction.y;
    if (divy >= 0) {
        tymin = (min_point.y - ray.origin.y) * divy;
        tymax = (max_point.y - ray.origin.y) * divy;
    } else {
        tymin = (max_point.y - ray.origin.y) * divy;
        tymax = (min_point.y - ray.origin.y) * divy;
    }
    if ((tmin > tymax) || (tymin > tmax))
        return false;
    if (tymin > tmin)
        tmin = tymin;
    if (tymax < tmax)
        tmax = tymax;

    divz = 1 / ray.direction.z;
    if (divz >= 0) {
        tzmin = (min_point.z - ray.origin.z) * divz;
        tzmax = (max_point.z - ray.origin.z) * divz;
    } else {
        tzmin = (max_point.z - ray.origin.z) * divz;
        tzmax = (min_point.z - ray.origin.z) * divz;
    }
    if ((tmin > tzmax) || (tzmin > tmax))
        return false;
    if (tzmin > tmin)
        tmin = tzmin;
    if (tzmax < tmax)
        tmax = tzmax;

    if ((tmin < T_MAX) && (tmax > T_MIN)) {
        hit_count += 1;
        return true;
    }
    return false;
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

Surface::Surface(Box bounding_box) : bounding_box(std::move(bounding_box)) {
}

SurfaceList::SurfaceList(std::vector<std::shared_ptr<Surface>> surfaces)
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
    : Surface(Box(tmath::vec3f(center.x - radius, center.y - radius, center.z - radius),
                  tmath::vec3f(center.x + radius, center.y + radius, center.z + radius))),
      center(std::move(center)), radius(radius), material(std::move(material)) {
}

std::optional<HitRecord> Sphere::hit(const Ray &ray) const {
    if (bounding_box.hit(ray)) {
        test_count += 1;
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
                hit_count += 1;
                return hit_record;
            }
        }
    }
    return {};
}

Face::Face(tmath::vec3f v1, tmath::vec3f v2, tmath::vec3f v3)
    : v1(std::move(v1)), v2(std::move(v2)), v3(std::move(v3)) {

    bounding_box.update(v1);
    bounding_box.update(v2);
    bounding_box.update(v3);
}

Triangle::Triangle(Face face, Material material)
    : face(std::move(face)), material(std::move(material)) {

    // TODO: Bounding box is duplicate between face and triangle
    bounding_box.update(face.v1);
    bounding_box.update(face.v2);
    bounding_box.update(face.v3);
}

std::optional<HitRecord> Face::hit(const Ray &ray) const {
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
        // material to be set on parent object
        hit_record.material = nullptr;
        hit_count += 1;
        return hit_record;
    }

    return {};
}

std::optional<HitRecord> Triangle::hit(const Ray &ray) const {
    auto hit_record = face.hit(ray);
    hit_record->material = &(this->material);
    return hit_record;
}

// TODO: possibly add axis as argument
Mesh::Mesh(std::vector<std::shared_ptr<Face>> triangles, Material material)
    : triangles(std::move(triangles)), material(std::move(material)) {

    for (const auto &triangle : this->triangles) {
        bounding_box.update(triangle->v1);
        bounding_box.update(triangle->v2);
        bounding_box.update(triangle->v3);
    }

    std::vector<std::shared_ptr<Surface>> faces(this->triangles.begin(), this->triangles.end());
    bvh = BVH(faces, Axis::X);
}

std::optional<HitRecord> Mesh::hit(const Ray &ray) const {
    float t = std::numeric_limits<float>::max();
    std::optional<HitRecord> hit_record;

    auto intersected = bvh.hit(ray);
    if (intersected && intersected->t < t) {
        hit_record = intersected;
        t = intersected->t;
        hit_record->material = &(this->material);
    }

    return hit_record;
}

inline Axis next_axis(Axis axis) {
    return Axis((static_cast<int>(axis) + 1) % 3);
}

std::optional<HitRecord> BVH::hit(const Ray &ray) const {
    if (bounding_box.hit(ray)) {
        std::optional<HitRecord> right_hit;
        std::optional<HitRecord> left_hit;

        if (right)
            right_hit = right->hit(ray);
        if (left)
            left_hit = left->hit(ray);

        if (left_hit && right_hit) {
            if (right_hit->t < left_hit->t)
                return right_hit;
            else
                return left_hit;
        }

        if (left_hit)
            return left_hit;

        if (right_hit)
            return right_hit;
    }

    return {};
}

BVH::BVH(const std::vector<std::shared_ptr<Surface>> &surfaces, Axis axis) {
    std::size_t len = surfaces.size();
    if (len == 0) {
        left = nullptr;
        right = nullptr;
    } else if (len == 1) {
        left = surfaces[0];
        right = nullptr;
        bounding_box = left->bounding_box;
    } else if (len == 2) {
        left = surfaces[0];
        right = surfaces[1];
        bounding_box = Box(left->bounding_box, right->bounding_box);
    } else {
        for (const auto &surface : surfaces) {
            bounding_box.update(surface->bounding_box.min_point);
            bounding_box.update(surface->bounding_box.max_point);
        }

        std::size_t half_len = len / 2;

        auto surfaces_copy = surfaces;
        std::sort(surfaces_copy.begin(), surfaces_copy.end(),
                  [axis](std::shared_ptr<Surface> s1, std::shared_ptr<Surface> s2) {

                      auto s1_midpoint = (s1->bounding_box.min_point + s1->bounding_box.max_point) / 2;
                      auto s2_midpoint = (s2->bounding_box.min_point + s2->bounding_box.max_point) / 2;

                      float s1_mid, s2_mid;
                      if (axis == Axis::X) {
                          s1_mid = s1_midpoint.x;
                          s2_mid = s2_midpoint.x;
                      } else if (axis == Axis::Y) {
                          s1_mid = s1_midpoint.y;
                          s2_mid = s2_midpoint.y;
                      } else {
                          s1_mid = s1_midpoint.z;
                          s2_mid = s2_midpoint.z;
                      }

                      return s1_mid < s2_mid;
                  });
        std::vector<std::shared_ptr<Surface>> left_surfaces(surfaces_copy.begin(),
                                                            surfaces_copy.begin() + half_len);
        std::vector<std::shared_ptr<Surface>> right_surfaces(surfaces_copy.begin() + half_len,
                                                             surfaces_copy.end());

        left = std::make_shared<BVH>(left_surfaces, next_axis(axis));
        right = std::make_shared<BVH>(right_surfaces, next_axis(axis));
    }
}

