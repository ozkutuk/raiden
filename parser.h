#ifndef __HW1__PARSER__
#define __HW1__PARSER__

#include "tinymath.h"
#include <string>
#include <vector>

namespace parser {
// Notice that all the structures are as simple as possible
// so that you are not enforced to adopt any style or design.

using Vec3f = tmath::vec3f;

struct Vec3i {
    int x, y, z;
};

struct Vec4f {
    float x, y, z, w;
};

struct Camera {
    Vec3f position;
    Vec3f gaze;
    Vec3f up;
    Vec4f near_plane;
    float near_distance;
    int image_width, image_height;
    std::string image_name;
};

struct Light {
    Vec3f position;
    Vec3f intensity; // Actually radiance for non-point lights,
                     // but we need not be too pedantic.
    enum class LightType { POINT, AREA }; // TODO: just use inheritance for this

    LightType type;
};

struct PointLight : public Light { 
};


struct AreaLight : public Light {
    Vec3f normal;
    float size;
};

struct Material {
    Vec3f ambient;
    Vec3f diffuse;
    Vec3f specular;
    Vec3f mirror;
    float phong_exponent;
    float refractive_index;
    Vec3f transparency;
    float roughness;
};

struct Face {
    int v0_id;
    int v1_id;
    int v2_id;
};

struct Mesh {
    int material_id;
    std::vector<Face> faces;
};

struct Triangle {
    int material_id;
    Face indices;
};

struct Sphere {
    int material_id;
    int center_vertex_id;
    float radius;
};

struct Scene {
    // Data
    Vec3f background_color;
    float shadow_ray_epsilon;
    int max_recursion_depth;
    std::vector<Camera> cameras;
    Vec3f ambient_light;
    std::vector<PointLight> point_lights;
    std::vector<AreaLight> area_lights;
    std::vector<Material> materials;
    std::vector<Vec3f> vertex_data;
    std::vector<Mesh> meshes;
    std::vector<Triangle> triangles;
    std::vector<Sphere> spheres;

    // Functions
    void loadFromXml(const std::string &filepath);
};
} // namespace parser

#endif
