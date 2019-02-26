#define STB_IMAGE_WRITE_IMPLEMENTATION

// stb_image requires this to be defined
// with Microsoft compilers
#ifdef _MSC_VER
#define STBI_MSC_SECURE_CRT
#endif

#include "stb_image_write.h"

#include <cstdint>
#include <vector>

// TODO: enforce 4-byte alignment
struct Color {
  uint8_t r;
  uint8_t g;
  uint8_t b;
};

int main(void) {
  constexpr int width = 1280;
  constexpr int height = 720;

  std::vector<Color> image(width * height);
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      double s = static_cast<double>(x) / width;
      s = static_cast<uint8_t>(s * 255);
      double t = static_cast<double>(y) / height;
      t = static_cast<uint8_t>(t * 255);
      Color c;
      c.r = s;
      c.g = t;
      image[y * width + x] = c;
    }
  }
  stbi_write_png("out.png", width, height, 3, image.data(),
                 width * sizeof(Color));
  return EXIT_SUCCESS;
}