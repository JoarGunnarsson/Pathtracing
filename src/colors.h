#ifndef colors
#define colors
#include "vec3.h"


namespace colors{
    vec3 const WHITE = vec3(1,1,1);
    vec3 const BLACK = vec3(0,0,0);
    vec3 const RED = vec3(1.0, 0.0, 0.0);
    vec3 const GREEN = vec3(0.0, 1.0, 0.0);
    vec3 const BLUE =  vec3(0.0, 0.0, 1.0);
    vec3 const YELLOW = vec3(1.0, 1.0, 0.0);
    vec3 const CYAN = vec3(0.0, 1.0, 1.0);
    vec3 const SKY_BLUE = vec3(0.251, 0.624, 0.769);
    vec3 const GREY = vec3(0.5, 0.5, 0.5);
    vec3 const BLUEISH_WHITE = vec3(0.7215, 0.8274, 0.8705);
    vec3 const WARM_WHITE = vec3(0.9922, 0.9569, 0.8627);
}

vec3 colorClip(vec3 rgb){
    return rgb / (rgb.max()+1);
}

#endif