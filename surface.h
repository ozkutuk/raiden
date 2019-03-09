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

class Surface {
  public:
    virtual std::optional<HitRecord> hit(const Ray &ray) const = 0;
};

class SurfaceList : public Surface {
  public:
    explicit SurfaceList(std::vector<std::unique_ptr<Surface>> surfaces);
    std::optional<HitRecord> hit(const Ray &ray) const override;

    std::vector<std::unique_ptr<Surface>> surfaces;
};

class Sphere : public Surface {
  public:
    explicit Sphere(tmath::vec3f center, float radius, Material material);
    std::optional<HitRecord> hit(const Ray &ray) const override;

    tmath::vec3f center;
    float radius;
    Material material;
};

class Triangle : public Surface {
  public:
    explicit Triangle(tmath::vec3f v1, tmath::vec3f v2, tmath::vec3f v3, Material material);
    std::optional<HitRecord> hit(const Ray &ray) const override;

    tmath::vec3f v1;
    tmath::vec3f v2;
    tmath::vec3f v3;
    Material material;
};

