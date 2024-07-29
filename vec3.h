#ifndef VEC3_H
#define VEC3_H
#include <cmath>
#include <iostream>


class vec3 {
    public:
        float e[3];

        vec3() : e{0,0,0} {}
        vec3(float e0, float e1, float e2) : e{e0, e1, e2} {}

        float x() const {
            return e[0]; 
        }

        float y() const { 
            return e[1]; 
        }
        
        float z() const { 
            return e[2];
        }

        vec3 operator-() const { 
            return vec3(-e[0], -e[1], -e[2]); 
        }

        float operator[](int i) const { 
            return e[i]; 
        }

        float& operator[](int i) { 
            return e[i]; 
        }

        float length() const {
            return std::sqrt(length_squared());
        }

        float length_squared() const {
            return e[0]*e[0] + e[1]*e[1] + e[2]*e[2];
        }
};


vec3 addVectors(vec3 v1, vec3 v2){
    return vec3(v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2]);
}


vec3 subtractVectors(vec3 v1, vec3 v2){
    return vec3(v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]);
}


vec3 multiplyVector(vec3 v, float value){
    return vec3(v[0] * value, v[1] * value, v[2] * value);
}


vec3 divideVector(vec3 v, float value){
    return vec3(v[0] / value, v[1] / value, v[2] / value);
}



float dotVectors(vec3 v1, vec3 v2){
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}


vec3 crossVectors(vec3 v1, vec3 v2){
    return vec3(v1[1] * v2[2] - v1[2] * v2[1], v1[2] * v2[0] - v1[0] * v2[2], v1[0] * v2[1] - v1[1] * v2[0]);
}


vec3 normalizeVector(vec3 v){
    float norm = v.length();
    return divideVector(v, norm);
}


void displayVector(vec3 v){
    std::clog << "[";
    for (int i = 0; i < 3; i++){
        std::clog << v[i];
        if (i != 2){
            std::clog << ", ";
        }
    }
    std::clog << "]\n";
}

#endif
