#define STB_IMAGE_WRITE_IMPLEMENTATION

// stb_image requires this to be defined
// with Microsoft compilers
#ifdef _MSC_VER
#define STBI_MSC_SECURE_CRT
#endif

#include "stb_image_write.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
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

float fresnel_reflection(const tmath::vec3f &incident, const tmath::vec3f &normal, float eta) {
    float cos_theta = std::abs(tmath::dot(incident, normal));
    float r_0 = ((eta - 1.0f) * (eta - 1.0f)) / ((eta + 1.0f) * (eta + 1.0f));
    float r = r_0 + (1.0f - r_0) * std::pow(1.0f - cos_theta, 5);
    return r;
}

tmath::vec3f cast_ray(const Ray &ray, const BVH &surfaces, const std::vector<Light> &lights,
                      const tmath::vec3f &ambient_light, const tmath::vec3f &background,
                      float epsilon, int recursion_depth) {

    auto intersected = surfaces.hit(ray);
    if (!intersected) {
        return background;
    }

    tmath::vec3f total_light = intersected->material->ambient * ambient_light;

    tmath::vec3f diffuse_light_intensity(0.0f, 0.0f, 0.0f);
    tmath::vec3f specular_light_intensity(0.0f, 0.0f, 0.0f);
    tmath::vec3f hit_point = ray.at(intersected->t);
    tmath::vec3f normal = intersected->normal;
    tmath::vec3f to_eye = tmath::normalize(ray.origin - hit_point);
    for (const auto &light : lights) {
        tmath::vec3f light_vec = light.position - hit_point;
        float light_distance = tmath::length(light_vec);
        // calculate this way, because we can reuse the distance
        tmath::vec3f light_dir = light_vec / light_distance;

        Ray shadow_ray(hit_point + epsilon * normal, light_dir);
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

    total_light += intersected->material->diffuse * diffuse_light_intensity +
                   intersected->material->specular * specular_light_intensity;

    if (intersected->material->refractive_index != 0 && recursion_depth > 0) {
        tmath::vec3f dir = -1.0f * to_eye;
        bool outside = tmath::dot(normal, dir) < 0;
        tmath::vec3f origin = hit_point;
        if (outside)
            origin -= epsilon * normal;
        else
            origin += epsilon * normal;

        Ray refracted(origin,
                      tmath::refract2(dir, normal, intersected->material->refractive_index));
        Ray reflected(hit_point + epsilon * normal, tmath::reflect(-1.0f * to_eye, normal));

        float reflection = fresnel_reflection(dir, normal, intersected->material->refractive_index);
        float refraction = 1.0f - reflection;

        // TODO: get rid of this awful re-calculation
        if (auto hit = surfaces.hit(refracted)) {
            tmath::vec3f attenuation =
                tmath::vec3f(std::pow(intersected->material->transparency.x, hit->t),
                             std::pow(intersected->material->transparency.y, hit->t),
                             std::pow(intersected->material->transparency.z, hit->t));

            total_light += refraction * attenuation *
                           cast_ray(refracted, surfaces, lights, ambient_light, background, epsilon,
                                    recursion_depth);
        }
        if (surfaces.hit(reflected))
            total_light += reflection * cast_ray(reflected, surfaces, lights, ambient_light,
                                                 background, epsilon, recursion_depth - 1);

    }

    else if (tmath::length(intersected->material->mirror_reflectance) != 0 && recursion_depth > 0) {
        Ray reflected(hit_point + epsilon * normal, tmath::reflect(-1.0f * to_eye, normal));
        // TODO: get rid of this awful re-calculation
        if (surfaces.hit(reflected))
            total_light += intersected->material->mirror_reflectance *
                           cast_ray(reflected, surfaces, lights, ambient_light, background, epsilon,
                                    recursion_depth - 1);
    }

    return total_light;
}

void render(const BVH &surfaces, const std::vector<Light> &lights,
            const tmath::vec3f &ambient_light, const tmath::vec3f &background, float epsilon,
            int max_recursion, const parser::Camera &camera) {
    int image_width = camera.image_width;
    int image_height = camera.image_height;

    tmath::vec3f origin = camera.position;
    tmath::vec3f gaze = tmath::normalize(camera.gaze);
    tmath::vec3f right = tmath::normalize(tmath::cross(gaze, camera.up));
    tmath::vec3f up = tmath::cross(right, gaze);

    float l, r, b, t;
    l = camera.near_plane.x;
    r = camera.near_plane.y;
    b = camera.near_plane.z;
    t = camera.near_plane.w;

    tmath::vec3f top_left = origin + camera.near_distance * gaze + l * right + t * up;

    float rl = r - l;
    float tb = t - b;

    std::vector<Color> image(image_width * image_height);

    for (int y = 0; y < image_height; y++) {
        for (int x = 0; x < image_width; x++) {

            float su = (rl * (x + 0.5f)) / image_width;
            float sv = (tb * (y + 0.5f)) / image_height;

            Ray r(origin, tmath::normalize((top_left + su * right - sv * up) - origin));
            tmath::vec3f value = tmath::clamp(
                cast_ray(r, surfaces, lights, ambient_light, background, epsilon, max_recursion),
                0.0f, 255.0f);
            Color c;
            c.r = static_cast<uint8_t>(value.x);
            c.g = static_cast<uint8_t>(value.y);
            c.b = static_cast<uint8_t>(value.z);

            image[y * image_width + x] = c;
        }
    }

    stbi_write_png(camera.image_name.c_str(), image_width, image_height, 3, image.data(),
                   image_width * sizeof(Color));
}

int main(int argc, char **argv) {

    if (argc != 2) {
        fprintf(stderr, "usage: ./raiden [scene_file]\n");
        return EXIT_FAILURE;
    }

    auto start = std::chrono::high_resolution_clock::now();

    parser::Scene scene;
    scene.loadFromXml(argv[1]);

    std::vector<std::shared_ptr<Surface>> surface_vector;
    for (const auto &sphere : scene.spheres) {
        Material m(scene.materials[sphere.material_id - 1].diffuse,
                   scene.materials[sphere.material_id - 1].specular,
                   scene.materials[sphere.material_id - 1].ambient,
                   scene.materials[sphere.material_id - 1].mirror,
                   scene.materials[sphere.material_id - 1].phong_exponent,
                   scene.materials[sphere.material_id - 1].refractive_index,
                   scene.materials[sphere.material_id - 1].transparency);
        std::shared_ptr<Surface> s = std::make_shared<Sphere>(
            scene.vertex_data[sphere.center_vertex_id - 1], sphere.radius, m);
        surface_vector.emplace_back(std::move(s));
    }

    for (const auto &triangle : scene.triangles) {
        Material m(scene.materials[triangle.material_id - 1].diffuse,
                   scene.materials[triangle.material_id - 1].specular,
                   scene.materials[triangle.material_id - 1].ambient,
                   scene.materials[triangle.material_id - 1].mirror,
                   scene.materials[triangle.material_id - 1].phong_exponent,
                   scene.materials[triangle.material_id - 1].refractive_index,
                   scene.materials[triangle.material_id - 1].transparency);
        std::shared_ptr<Surface> s =
            std::make_shared<Triangle>(Face(scene.vertex_data[triangle.indices.v0_id - 1],
                                            scene.vertex_data[triangle.indices.v1_id - 1],
                                            scene.vertex_data[triangle.indices.v2_id - 1]),
                                       m);
        surface_vector.emplace_back(std::move(s));
    }

    for (const auto &mesh : scene.meshes) {
        Material m(scene.materials[mesh.material_id - 1].diffuse,
                   scene.materials[mesh.material_id - 1].specular,
                   scene.materials[mesh.material_id - 1].ambient,
                   scene.materials[mesh.material_id - 1].mirror,
                   scene.materials[mesh.material_id - 1].phong_exponent,
                   scene.materials[mesh.material_id - 1].refractive_index,
                   scene.materials[mesh.material_id - 1].transparency);

        std::vector<Face> faces;
        for (const auto &face : mesh.faces) {
            faces.emplace_back(scene.vertex_data[face.v0_id - 1], scene.vertex_data[face.v1_id - 1],
                               scene.vertex_data[face.v2_id - 1]);
        }

        std::shared_ptr<Surface> s = std::make_shared<Mesh>(faces, m);
        surface_vector.emplace_back(std::move(s));
    }

    SurfaceList surfaces(std::move(surface_vector));
    BVH bvh(surfaces.surfaces, Axis::X);

    std::vector<Light> lights;
    for (const auto &light : scene.point_lights) {
        Light l;
        l.position = light.position;
        l.intensity = light.intensity;
        lights.emplace_back(l);
    }

    for (const auto &camera : scene.cameras)
        render(bvh, lights, scene.ambient_light, scene.background_color,
               scene.shadow_ray_epsilon, scene.max_recursion_depth, camera);

    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - start);

    std::cout << "Render time                    : " << elapsed.count() / 1000.0 << "s\n"
              << "# of ray-triangle tests        : " << Face::test_count << "\n"
              << "# of ray-triangle intersections: " << Face::hit_count << "\n"
              << "# of ray-sphere tests          : " << Sphere::test_count << "\n"
              << "# of ray-sphere intersections  : " << Sphere::hit_count << "\n"
              << "# of ray-box tests             : " << Box::test_count << "\n"
              << "# of ray-box intersections     : " << Box::hit_count << std::endl;
    return EXIT_SUCCESS;
}
