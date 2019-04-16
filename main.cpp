#define STB_IMAGE_WRITE_IMPLEMENTATION

// stb_image requires this to be defined
// with Microsoft compilers
#ifdef _MSC_VER
#define STBI_MSC_SECURE_CRT
#endif

#include "stb_image_write.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <limits>
#include <memory>
#include <tuple>
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


using parser::Light, parser::PointLight, parser::AreaLight; // TODO: reorganize this

inline double mrandom() {
    return ((double)rand() / (RAND_MAX));
}

float fresnel_reflection(const tmath::vec3f &incident, const tmath::vec3f &normal, float eta) {
    float cos_theta = std::abs(tmath::dot(incident, normal));
    float r_0 = ((eta - 1.0f) * (eta - 1.0f)) / ((eta + 1.0f) * (eta + 1.0f));
    float r = r_0 + (1.0f - r_0) * std::pow(1.0f - cos_theta, 5);
    return r;
}

std::tuple<tmath::vec3f, tmath::vec3f, tmath::vec3f> orthonormal_basis(const tmath::vec3f &vec) {
    tmath::vec3f w = tmath::normalize(vec);
    tmath::vec3f not_colinear = w;

    std::array<float, 3> magnitudes = {std::abs(not_colinear.x), std::abs(not_colinear.y),
                                       std::abs(not_colinear.z)};

    auto min_element = std::min_element(std::begin(magnitudes), std::end(magnitudes));
    std::size_t min_index = std::distance(std::begin(magnitudes), min_element);

    if (min_index == 0)
        not_colinear.x = 1.0f;
    else if (min_index == 1)
        not_colinear.y = 1.0f;
    else
        not_colinear.z = 1.0f;

    tmath::vec3f u = tmath::normalize(tmath::cross(not_colinear, w));
    tmath::vec3f v = tmath::cross(w, u);

    return {u, v, w};
}

tmath::vec3f to_canonical(const std::tuple<tmath::vec3f, tmath::vec3f, tmath::vec3f> &basis,
                          tmath::vec3f point) {
    auto [u, v, w] = basis;
    return tmath::vec3f(tmath::dot(u, point), tmath::dot(v, point), tmath::dot(w, point));
}

tmath::vec3f cast_ray(const Ray &ray, const BVH &surfaces, const std::vector<std::unique_ptr<Light>> &lights,
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
        tmath::vec3f light_intensity = light->intensity;
        tmath::vec3f light_vec = light->position - hit_point;
        float light_distance;
        tmath::vec3f light_dir;

        if (light->type == Light::LightType::POINT) {
            // calculate this way, because we can reuse the distance
            light_distance = tmath::length(light_vec);
            light_dir = light_vec / light_distance;
        } else if (light->type == Light::LightType::AREA) {
            auto area_ptr = static_cast<AreaLight*>(light.get()); // wow this is ugly
            const float light_size = area_ptr->size;
            const tmath::vec3f light_normal = area_ptr->normal;
            auto [u, v, w] = orthonormal_basis(light_normal);

            float u_offset = -(light_size / 2.0f) + light_size * mrandom();
            float v_offset = -(light_size / 2.0f) + light_size * mrandom();

            light_dir = light_vec + u_offset * u + v_offset * v;
            light_distance = tmath::length(light_dir);
            light_dir = light_dir / light_distance;
            light_intensity *=
                std::abs(tmath::dot(light_dir, light_normal)) * light_size * light_size;
        } else {
            std::cerr << "Unexpected light type" << std::endl;
            std::exit(EXIT_FAILURE);
        }

        Ray shadow_ray(hit_point + epsilon * normal, light_dir);
        auto shadow_hit = surfaces.hit(shadow_ray);

        if (shadow_hit && shadow_hit->t < light_distance)
            continue;

        diffuse_light_intensity += light_intensity * std::max(0.0f, tmath::dot(light_dir, normal)) /
                                   (light_distance * light_distance);

        // NOTE: this may be calculated without half vector
        tmath::vec3f half = tmath::normalize(light_dir + to_eye);
        float nh_angle = std::max(0.0f, tmath::dot(normal, half));
        specular_light_intensity += light_intensity *
                                    std::pow(nh_angle, intersected->material->specular_exponent) /
                                    (light_distance * light_distance);
    }

    total_light += intersected->material->diffuse * diffuse_light_intensity +
                   intersected->material->specular * specular_light_intensity;

    const bool glossy = false;
    const float blur_window = 0.5f;

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

        auto reflected_direction = tmath::reflect(-1.0f * to_eye, normal);
        if (glossy) {
            auto [u, v, w] = orthonormal_basis(reflected_direction);

            float u_offset = -(blur_window / 2.0f) + blur_window * mrandom();
            float v_offset = -(blur_window / 2.0f) + blur_window * mrandom();

            reflected_direction =
                tmath::normalize(reflected_direction + u_offset * u + v_offset * v);
        }

        Ray reflected(hit_point + epsilon * normal, reflected_direction);

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

        auto reflected_direction = tmath::reflect(-1.0f * to_eye, normal);
        if (glossy) {
            auto [u, v, w] = orthonormal_basis(reflected_direction);

            float u_offset = -(blur_window / 2.0f) + blur_window * mrandom();
            float v_offset = -(blur_window / 2.0f) + blur_window * mrandom();

            reflected_direction =
                tmath::normalize(reflected_direction + u_offset * u + v_offset * v);
        }

        Ray reflected(hit_point + epsilon * normal, reflected_direction);

        // TODO: get rid of this awful re-calculation
        if (surfaces.hit(reflected))
            total_light += intersected->material->mirror_reflectance *
                           cast_ray(reflected, surfaces, lights, ambient_light, background, epsilon,
                                    recursion_depth - 1);
    }

    return total_light;
}

void render(const BVH &surfaces, const std::vector<std::unique_ptr<Light>> &lights,
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

    srand(time(0));
    for (int y = 0; y < image_height; y++) {
        for (int x = 0; x < image_width; x++) {

            float su = rl / image_width;
            float sv = tb / image_height;

            const int n_samples = 1;
            const float bin_dimension = 1.0f / std::sqrt(n_samples);
            tmath::vec3f value;

            const bool dof = false;

            if (dof) {
                float aperture_size = 1.5f;
                float focal_distance = 21.0f;

                tmath::vec3f pixel_coordinate =
                    top_left + (su * (x + 0.5f) * right) - (sv * (y + 0.5f) * up);

                Ray ray_to_pixel(origin, tmath::normalize(pixel_coordinate - origin));
                tmath::vec3f focal_point = ray_to_pixel.at(focal_distance);

                for (int sub_y = 0; sub_y < std::sqrt(n_samples); sub_y++) {
                    for (int sub_x = 0; sub_x < std::sqrt(n_samples); sub_x++) {
                        float lens_origin_x_offset =
                            aperture_size * (sub_x + mrandom()) * bin_dimension;
                        float lens_origin_y_offset =
                            aperture_size * (sub_y + mrandom()) * bin_dimension;

                        tmath::vec3f pixel_top_left =
                            pixel_coordinate - tmath::vec3f(0.5f, 0.5f, 0.0f);
                        tmath::vec3f lens_origin = pixel_top_left + lens_origin_x_offset * right -
                                                   lens_origin_y_offset * up;

                        Ray r(lens_origin, tmath::normalize(focal_point - lens_origin));
                        value += cast_ray(r, surfaces, lights, ambient_light, background, epsilon,
                                          max_recursion);
                    }
                }
            } else {

                for (int sub_y = 0; sub_y < std::sqrt(n_samples); sub_y++) {
                    for (int sub_x = 0; sub_x < std::sqrt(n_samples); sub_x++) {
                        float ray_x = su * (x + (sub_x + mrandom()) * bin_dimension);
                        float ray_y = sv * (y + (sub_y + mrandom()) * bin_dimension);

                        Ray r(origin,
                              tmath::normalize((top_left + ray_x * right - ray_y * up) - origin));
                        value += cast_ray(r, surfaces, lights, ambient_light, background, epsilon,
                                          max_recursion);
                    }
                }
            }

            value = tmath::clamp(value / n_samples, 0.0f, 255.0f);

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

        std::vector<std::shared_ptr<Face>> faces;
        for (const auto &face : mesh.faces) {
            faces.emplace_back(std::make_shared<Face>(scene.vertex_data[face.v0_id - 1],
                                                      scene.vertex_data[face.v1_id - 1],
                                                      scene.vertex_data[face.v2_id - 1]));
        }

        std::shared_ptr<Surface> s = std::make_shared<Mesh>(faces, m);
        surface_vector.emplace_back(std::move(s));
    }

    SurfaceList surfaces(std::move(surface_vector));
    BVH bvh(surfaces.surfaces, Axis::X);

    std::vector<std::unique_ptr<Light>> lights;
    for (const auto &light : scene.point_lights) {
        PointLight l;
        l.position = light.position;
        l.intensity = light.intensity;
        l.type = Light::LightType::POINT;
        lights.emplace_back(std::make_unique<PointLight>(l));
    }
    for (const auto &light : scene.area_lights) {
        AreaLight l;
        l.position = light.position;
        l.intensity = light.intensity;
        l.normal = light.normal;
        l.size = light.size;
        l.type = Light::LightType::AREA;
        lights.emplace_back(std::make_unique<AreaLight>(l));
    }

    for (const auto &camera : scene.cameras)
        render(bvh, lights, scene.ambient_light, scene.background_color, scene.shadow_ray_epsilon,
               scene.max_recursion_depth, camera);

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
