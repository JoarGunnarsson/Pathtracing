#ifndef COLORS_H
#define COLORS_H
#include "vec3.h"


namespace colors{
    const vec3 WHITE = vec3(1,1,1);
    const vec3 BLACK = vec3(0,0,0);
    const vec3 RED = vec3(1.0, 0.0, 0.0);
    const vec3 GREEN = vec3(0.0, 1.0, 0.0);
    const vec3 BLUE =  vec3(0.0, 0.0, 1.0);
    const vec3 YELLOW = vec3(1.0, 1.0, 0.0);
    const vec3 CYAN = vec3(0.0, 1.0, 1.0);
    const vec3 SKY_BLUE = vec3(0.251, 0.624, 0.769);
    const vec3 GREY = vec3(0.5, 0.5, 0.5);
    const vec3 BLUEISH_WHITE = vec3(0.7215, 0.8274, 0.8705);
    const vec3 WARM_WHITE = vec3(0.9922, 0.9569, 0.8627);
    const vec3 GOLD = vec3(1, 0.84, 0);
    const vec3 CHOCOLATE_BROWN = vec3(0.25490196078, 0.0980392156, 0) * 0.5;
}


inline vec3 tone_map(vec3& rgb){
    return rgb / (rgb.max() + 1.0);
}

#endif
