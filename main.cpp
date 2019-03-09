#define STB_IMAGE_WRITE_IMPLEMENTATION

// stb_image requires this to be defined
// with Microsoft compilers
#ifdef _MSC_VER
#define STBI_MSC_SECURE_CRT
#endif

#include "stb_image_write.h"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <memory>
#include <vector>

#include "parser.h"
#include "ray.h"
#include "surface.h"
#include "tinymath.h"

struct Color {
    uint8_t r;
    uint8_t g;
    uint8_t b;
};

struct Light {
    tmath::vec3f intensity;
    tmath::vec3f position;
};

tmath::vec3f cast_ray(const Ray &ray, const SurfaceList &surfaces, const std::vector<Light> &lights,
                      const tmath::vec3f &ambient_light, int recursion_depth) {

    auto intersected = surfaces.hit(ray);
    if (!intersected) {
        // TODO: get bg color from scene
        return tmath::vec3f(0.0f, 0.0f, 0.0f);
    }

    tmath::vec3f total_light = intersected->material->ambient * ambient_light;

    tmath::vec3f diffuse_light_intensity(0.0f, 0.0f, 0.0f);
    tmath::vec3f specular_light_intensity(0.0f, 0.0f, 0.0f);
    tmath::vec3f hit_point = ray.at(intersected->t);
    // tmath::vec3f normal = tmath::normalize(hit_point - hit.sphere->center);
    tmath::vec3f normal = intersected->normal;
    tmath::vec3f to_eye = tmath::normalize(ray.origin - hit_point);
    for (const auto &light : lights) {
        tmath::vec3f light_vec = light.position - hit_point;
        float light_distance = tmath::length(light_vec);
        // calculate this way, because we can reuse the distance
        tmath::vec3f light_dir = light_vec / light_distance;

        // TODO: get epsilon from xml
        Ray shadow_ray(hit_point + 1e-3 * normal, light_dir);
        // HitRecord shadow_hit = scene_hit(shadow_ray, spheres);
        auto shadow_hit = surfaces.hit(shadow_ray);

        if (shadow_hit && shadow_hit->t < light_distance)
            continue;

        diffuse_light_intensity += light.intensity * std::max(0.0f, tmath::dot(light_dir, normal)) /
                                   (light_distance * light_distance);

        // NOTE: this may be calculated without half vector
        tmath::vec3f half = tmath::normalize(light_dir + to_eye);
        float nh_angle = std::max(0.0f, tmath::dot(normal, half));
        specular_light_intensity += light.intensity *
                                    std::pow(nh_angle, intersected->material->specular_exponent) /
                                    (light_distance * light_distance);
    }

    total_light = intersected->material->diffuse * diffuse_light_intensity +
                  intersected->material->specular * specular_light_intensity;

    if (tmath::length(intersected->material->mirror_reflectance) != 0 && recursion_depth > 0) {
        Ray reflected(hit_point + 1e-3 * normal, tmath::reflect(-1.0f * to_eye, normal));
        // TODO: get rid of this awful re-calculation
        if (surfaces.hit(reflected))
            total_light +=
                intersected->material->mirror_reflectance *
                cast_ray(reflected, surfaces, lights, ambient_light, recursion_depth - 1);
    }

    return total_light;
}

void render(const SurfaceList &surfaces, const std::vector<Light> &lights,
            const tmath::vec3f &ambient_light, const parser::Camera &camera) {
    int image_width = camera.image_width;
    int image_height = camera.image_height;

    tmath::vec3f origin = camera.position;
    tmath::vec3f right = tmath::normalize(tmath::cross(camera.gaze, camera.up));

    float l, r, b, t;
    l = camera.near_plane.x;
    r = camera.near_plane.y;
    b = camera.near_plane.z;
    t = camera.near_plane.w;

    tmath::vec3f top_left = origin + camera.near_distance * camera.gaze + l * right + t * camera.up;

    float rl = r - l;
    float tb = t - b;

    std::vector<Color> image(image_width * image_height);

    for (int y = 0; y < image_height; y++) {
        for (int x = 0; x < image_width; x++) {

            float su = (rl * (x + 0.5f)) / image_width;
            float sv = (tb * (y + 0.5f)) / image_height;

            Ray r(origin, tmath::normalize((top_left + su * right - sv * camera.up) - origin));
            // TODO: get recursion depth from xml
            tmath::vec3f value =
                tmath::clamp(cast_ray(r, surfaces, lights, ambient_light, 6), 0.0f, 255.0f);
            Color c;
            c.r = static_cast<uint8_t>(value.x);
            c.g = static_cast<uint8_t>(value.y);
            c.b = static_cast<uint8_t>(value.z);

            image[y * image_width + x] = c;
        }
    }

    stbi_write_png("out.png", image_width, image_height, 3, image.data(),
                   image_width * sizeof(Color));
}

int main(int argc, char **argv) {

    if (argc != 2) {
        fprintf(stderr, "usage: ./raiden [scene_file]\n");
        return EXIT_FAILURE;
    }

    parser::Scene scene;
    scene.loadFromXml(argv[1]);

    std::vector<std::unique_ptr<Surface>> surface_vector;
    for (const auto &sphere : scene.spheres) {
        Material m(scene.materials[sphere.material_id - 1].diffuse,
                   scene.materials[sphere.material_id - 1].specular,
                   scene.materials[sphere.material_id - 1].ambient,
                   scene.materials[sphere.material_id - 1].mirror,
                   scene.materials[sphere.material_id - 1].phong_exponent);
        std::unique_ptr<Surface> s = std::make_unique<Sphere>(
            scene.vertex_data[sphere.center_vertex_id - 1], sphere.radius, m);
        surface_vector.emplace_back(std::move(s));
    }

    for (const auto &triangle : scene.triangles) {
        Material m(scene.materials[triangle.material_id - 1].diffuse,
                   scene.materials[triangle.material_id - 1].specular,
                   scene.materials[triangle.material_id - 1].ambient,
                   scene.materials[triangle.material_id - 1].mirror,
                   scene.materials[triangle.material_id - 1].phong_exponent);
        std::unique_ptr<Surface> s =
            std::make_unique<Triangle>(scene.vertex_data[triangle.indices.v0_id - 1],
                                       scene.vertex_data[triangle.indices.v1_id - 1],
                                       scene.vertex_data[triangle.indices.v2_id - 1], m);
        surface_vector.emplace_back(std::move(s));
    }

    std::vector<Triangle> triangles;
    for (const auto &mesh : scene.meshes) {
        Material m(scene.materials[mesh.material_id - 1].diffuse,
                   scene.materials[mesh.material_id - 1].specular,
                   scene.materials[mesh.material_id - 1].ambient,
                   scene.materials[mesh.material_id - 1].mirror,
                   scene.materials[mesh.material_id - 1].phong_exponent);

        for (const auto &face : mesh.faces) {
            triangles.emplace_back(Triangle(scene.vertex_data[face.v0_id - 1],
                                            scene.vertex_data[face.v1_id - 1],
                                            scene.vertex_data[face.v2_id - 1], m));
        }

        std::unique_ptr<Surface> s = std::make_unique<Mesh>(triangles, m);
        surface_vector.emplace_back(std::move(s));
    }

    SurfaceList surfaces(std::move(surface_vector));

    std::vector<Light> lights;
    for (const auto &light : scene.point_lights) {
        Light l;
        l.position = light.position;
        l.intensity = light.intensity;
        lights.emplace_back(l);
    }

    // TODO: render for each camera in xml
    render(surfaces, lights, scene.ambient_light, scene.cameras[0]);
    return EXIT_SUCCESS;
}
