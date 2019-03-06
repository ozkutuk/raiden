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
#include <vector>

#include "ray.h"
#include "sphere.h"
#include "tinymath.h"

struct Color {
    uint8_t r;
    uint8_t g;
    uint8_t b;
};

struct HitRecord {
    const Sphere *sphere;
    float t;
};

struct Light {
    tmath::vec3f intensity;
    tmath::vec3f position;
};

HitRecord scene_hit(const Ray &ray, const std::vector<Sphere> &spheres) {
    float distance = std::numeric_limits<float>::max();
    const Sphere *closest = nullptr;
    for (const auto &sphere : spheres) {
        float t = sphere.hit(ray);
        if (t > 0 && t < distance) {
            distance = t;
            closest = &sphere;
        }
    }

    HitRecord hit;
    if (distance < 1000)
        hit.sphere = closest;
    else
        hit.sphere = nullptr;
    hit.t = distance;
    return hit;
}

tmath::vec3f cast_ray(const Ray &ray, const std::vector<Sphere> &spheres,
                      const std::vector<Light> &lights) {

    HitRecord hit = scene_hit(ray, spheres);
    if (!hit.sphere) {
        tmath::vec3f unit = tmath::normalize(ray.direction);
        float k = 0.5f * (unit.y + 1.0f);
        tmath::vec3f white = tmath::vec3f(1.0f, 1.0f, 1.0f);
        tmath::vec3f blue = tmath::vec3f(0.5f, 0.7f, 1.0f);

        return (1.0f - k) * white + k * blue;
    }

    tmath::vec3f diffuse_light_intensity(0.0f, 0.0f, 0.0f);
    tmath::vec3f hit_point = ray.at(hit.t);
    tmath::vec3f normal = tmath::normalize(hit_point - hit.sphere->center);
    for (const auto &light : lights) {
        tmath::vec3f light_dir = tmath::normalize(light.position - hit_point);
        diffuse_light_intensity += light.intensity * std::max(0.0f, tmath::dot(light_dir, normal));
    }

    return hit.sphere->material.diffuse * diffuse_light_intensity;
}

void render(const std::vector<Sphere> &spheres, const std::vector<Light> &lights) {
    constexpr int image_width = 1920;
    constexpr int image_height = 1080;

    // NOTE: assume camera is at 0,0,0 for now
    // TODO: create a proper camera
    tmath::vec3f origin(0.0f, 0.0f, 0.0f);
    tmath::vec3f top_left(-8.0f, 4.5f, -4.0f);
    tmath::vec3f horizontal(16.0f, 0.0f, 0.0f);
    tmath::vec3f vertical(0.0f, -9.0f, 0.0f);
    tmath::vec3f u = horizontal / image_width;
    tmath::vec3f v = vertical / image_height;

    std::vector<Color> image(image_width * image_height);

    for (int y = 0; y < image_height; y++) {
        for (int x = 0; x < image_width; x++) {
            Ray r(origin, top_left + x * u + y * v);
            tmath::vec3f value = cast_ray(r, spheres, lights);
            Color c;
            c.r = static_cast<uint8_t>(value.x * 255.99f);
            c.g = static_cast<uint8_t>(value.y * 255.99f);
            c.b = static_cast<uint8_t>(value.z * 255.99f);

            image[y * image_width + x] = c;
        }
    }

    stbi_write_png("out.png", image_width, image_height, 3, image.data(),
                   image_width * sizeof(Color));
}

int main(void) {
    Material darkred(tmath::vec3f(0.7f, 0.0f, 0.3f));
    Material beige(tmath::vec3f(0.9f, 0.9f, 0.7f));
    Sphere sphere(tmath::vec3f(-4.0f, 0.0f, -7.0f), 1.75f, darkred);
    Sphere sphere2(tmath::vec3f(4.0f, 2.0f, -9.0f), 2.0f, beige);
    std::vector<Sphere> spheres = {sphere, sphere2};

    Light light;
    light.position = tmath::vec3f(-6.0f, 4.0f, -5.0f);
    tmath::vec3f white = tmath::vec3f(0.9f, 0.9f, 0.9f);
    light.intensity = white;
    std::vector<Light> lights = {light};

    render(spheres, lights);
    return EXIT_SUCCESS;
}
