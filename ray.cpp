#include "ray.h"

#include <utility>

Ray::Ray(tmath::vec3f origin, tmath::vec3f direction)
    : origin(std::move(origin)), direction(std::move(direction)) {
}

tmath::vec3f Ray::at(float t) const {
    return origin + t * direction;
}
