#define STB_IMAGE_WRITE_IMPLEMENTATION

// stb_image requires this to be defined
// with Microsoft compilers
#ifdef _MSC_VER
#define STBI_MSC_SECURE_CRT
#endif

#include "stb_image_write.h"

#include <cstdint>
#include <vector>

#include "ray.h"
#include "tinymath.h"

struct Color {
    uint8_t r;
    uint8_t g;
    uint8_t b;
};

struct Sphere {
    tmath::vec3f center;
    float radius;
};

bool hit_sphere(const Sphere &sphere, const Ray &ray) {
    float a = tmath::dot(ray.direction, ray.direction);
    tmath::vec3f o_minus_c = ray.origin - sphere.center;
    float b = 2 * tmath::dot(o_minus_c, ray.direction);
    float c = tmath::dot(o_minus_c, o_minus_c) - sphere.radius * sphere.radius;

    float discriminant = b * b - 4 * a * c;
    return discriminant > 0;
}

tmath::vec3f color(const Ray &ray) {
    Sphere sphere;
    sphere.center = tmath::vec3f(0.0f, 0.0f, -1.0f);
    sphere.radius = 0.5f;

    if (hit_sphere(sphere, ray))
        return tmath::vec3f(1.0f, 0.0f, 0.0f);

    tmath::vec3f unit = tmath::normalize(ray.direction);
    float t = 0.5f * (unit.y + 1.0f);
    tmath::vec3f white = tmath::vec3f(1.0f, 1.0f, 1.0f);
    tmath::vec3f blue = tmath::vec3f(0.5f, 0.7f, 1.0f);

    return (1.0f - t) * white + t * blue;
}

int main(void) {
    constexpr int image_width = 1024;
    constexpr int image_height = 512;

    // assume camera is at 0,0,0
    tmath::vec3f origin(0.0f, 0.0f, 0.0f);
    tmath::vec3f top_left(-2.0f, 1.0f, -1.0f);
    tmath::vec3f horizontal(4.0f, 0.0f, 0.0f);
    tmath::vec3f vertical(0.0f, -2.0f, 0.0f);
    tmath::vec3f u = horizontal / image_width;
    tmath::vec3f v = vertical / image_height;

    std::vector<Color> image(image_width * image_height);

    for (int y = 0; y < image_height; y++) {
        for (int x = 0; x < image_width; x++) {
            Ray r(origin, top_left + x * u + y * v);
            tmath::vec3f value = color(r);
            Color c;
            c.r = static_cast<uint8_t>(value.x * 255.99f);
            c.g = static_cast<uint8_t>(value.y * 255.99f);
            c.b = static_cast<uint8_t>(value.z * 255.99f);

            image[y * image_width + x] = c;
        }
    }

    stbi_write_png("out.png", image_width, image_height, 3, image.data(),
                   image_width * sizeof(Color));
    return EXIT_SUCCESS;
}
