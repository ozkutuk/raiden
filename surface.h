#pragma once

#include <memory>
#include <optional>
#include <vector>

#include "material.h"
#include "tinymath.h"

class Ray;

constexpr float T_MIN = 0.0f;
constexpr float T_MAX = 1000.0f;

class HitRecord {
  public:
    float t;
    tmath::vec3f normal;
    const Material *material;
};

class Box {
  public:
    explicit Box(tmath::vec3f min_point, tmath::vec3f max_point);
    explicit Box(const Box &left, const Box &right);
    Box();
    bool hit(const Ray &ray) const;
    void update(const tmath::vec3f &point);

    tmath::vec3f min_point;
    tmath::vec3f max_point;

    static uint64_t test_count;
    static uint64_t hit_count;
};

class Surface {
  public:
    explicit Surface(Box bounding_box);
    Surface() = default;
    virtual std::optional<HitRecord> hit(const Ray &ray) const = 0;

    Box bounding_box;
};

class SurfaceList : public Surface {
  public:
    explicit SurfaceList(std::vector<std::shared_ptr<Surface>> surfaces);
    std::optional<HitRecord> hit(const Ray &ray) const override;

    std::vector<std::shared_ptr<Surface>> surfaces;
};

enum class Axis {
    X,
    Y,
    Z
};

class BVH : public Surface {
  public:
    explicit BVH(const std::vector<std::shared_ptr<Surface>> &surfaces, Axis axis);
    std::optional<HitRecord> hit(const Ray &ray) const override;
    
    std::shared_ptr<Surface> left;
    std::shared_ptr<Surface> right;
};

class Sphere : public Surface {
  public:
    explicit Sphere(tmath::vec3f center, float radius, Material material);
    std::optional<HitRecord> hit(const Ray &ray) const override;

    tmath::vec3f center;
    float radius;
    Material material;

    static uint64_t test_count;
    static uint64_t hit_count;
};

class Face : public Surface {
  public:
    explicit Face(tmath::vec3f v1, tmath::vec3f v2, tmath::vec3f v3);
    std::optional<HitRecord> hit(const Ray &ray) const override;
    Face() = delete;

    tmath::vec3f v1;
    tmath::vec3f v2;
    tmath::vec3f v3;

    // benchmark stuff
    static uint64_t test_count;
    static uint64_t hit_count;
};

class Triangle : public Surface {
  public:
    explicit Triangle(Face face, Material material);
    std::optional<HitRecord> hit(const Ray &ray) const override;

    Face face;
    Material material;
};

class Mesh : public Surface {
  public:
    explicit Mesh(std::vector<Face> triangles, Material material);
    std::optional<HitRecord> hit(const Ray &ray) const override;

    std::vector<Face> triangles;
    Material material;
};
