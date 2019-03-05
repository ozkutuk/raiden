#include "ray.h"

Ray::Ray(tmath::vec3f origin, tmath::vec3f direction)
    : origin(origin), direction(direction) {
}

tmath::vec3f Ray::at(float t) {
  return origin + t * direction;
}
