#ifndef VEC3_H
#define VEC3_H
#include <cmath>
#include <iostream>


class vec3 {
    public:
        double e[3];

        vec3() : e{0,0,0} {}
        vec3(double e0, double e1, double e2) : e{e0, e1, e2} {}

        inline double x() const {
            return e[0]; 
        }

        inline double y() const { 
            return e[1]; 
        }
        
        inline double z() const { 
            return e[2];
        }

        vec3 operator-() const { 
            return vec3(-e[0], -e[1], -e[2]); 
        }

        inline double operator[](int i) const { 
            return e[i]; 
        }

        inline double length() const {
            return std::sqrt(length_squared());
        }

        inline double length_squared() const {
            return e[0]*e[0] + e[1]*e[1] + e[2]*e[2];
        }

        inline double max() const{
            return std::max(std::max(e[0], e[1]), e[2]);
        }
};


inline vec3 addVectors(const vec3& v1, const vec3& v2){
    return vec3(v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2]);
}


inline vec3 subtractVectors(const vec3& v1, const vec3& v2){
    return vec3(v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]);
}


inline vec3 multiplyVector(const vec3& v, double value){
    return vec3(v[0] * value, v[1] * value, v[2] * value);
}

inline vec3 multiplyVectorElementwise(const vec3& v1, const vec3& v2){
    return vec3(v1[0] * v2[0], v1[1] * v2[1], v1[2] * v2[2]);
}


inline vec3 divideVector(const vec3& v, double value){
    return multiplyVector(v, 1/value);
}



inline double dotVectors(const vec3& v1, const vec3& v2){
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}


vec3 crossVectors(const vec3& v1, const vec3& v2){
    return vec3(v1[1] * v2[2] - v1[2] * v2[1], v1[2] * v2[0] - v1[0] * v2[2], v1[0] * v2[1] - v1[1] * v2[0]);
}


inline vec3 normalizeVector(const vec3& v){
    double norm = v.length();
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
