#ifndef COLORS_H
#define COLORS_H
#include "vec3.h"

namespace colors {
const vec3 WHITE = vec3(1.0, 1.0, 1.0);
const vec3 BLACK = vec3(0.0, 0.0, 0.0);
const vec3 RED = vec3(1.0, 0.0, 0.0);
const vec3 GREEN = vec3(0.0, 1.0, 0.0);
const vec3 BLUE = vec3(0.0, 0.0, 1.0);
const vec3 YELLOW = vec3(1.0, 1.0, 0.0);
const vec3 CYAN = vec3(0.0, 1.0, 1.0);
const vec3 SKY_BLUE = vec3(0.251, 0.624, 0.769);
const vec3 GREY = vec3(0.5, 0.5, 0.5);
const vec3 BLUEISH_WHITE = vec3(0.7215, 0.8274, 0.8705);
const vec3 WARM_WHITE = vec3(0.9922, 0.9569, 0.8627);
const vec3 GOLD = vec3(1, 0.84, 0);
const vec3 CHOCOLATE_BROWN = vec3(0.25490196078, 0.0980392156, 0) * 0.5;
}

inline double apply_gamma_correction(const double x) {
    double threshold = 0.04045;
    return x <= threshold ? x / 12.92 : std::pow((x + 0.055) / 1.055, 2.4);
}

inline vec3 apply_gamma_correction(const vec3& rgb) {
    vec3 out = vec3(0);
    for (int i = 0; i < 3; i++) {
        out[i] = apply_gamma_correction(rgb[i]);
    }
    return out;
}

#endif
