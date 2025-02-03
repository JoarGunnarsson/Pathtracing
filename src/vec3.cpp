#include "vec3.h"


void display_vector(const vec3& v){
    std::cout << "[";
    for (int i = 0; i < 3; i++){
        std::cout << v[i];
        if (i != 2){
            std::cout << ", ";
        }
    }
    std::cout << "]\n";
}
