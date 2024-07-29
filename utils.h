#ifndef utils
#define utils
#include "vec3.h"
#include "Sphere.h"


struct Ray{
    vec3 startingPosition;
    vec3 directionVector;
};


float solveQuadratic(float b, float c){
    float discriminant = pow(b, 2) - 4 * c;
    if (discriminant < 0){
        return -1.0;
    }
    float root_discriminant = sqrt(discriminant);
    float maximum_solution = - 1.0 / 2.0 * (b - root_discriminant);
    float minimum_solution = - 1.0 / 2.0 * (b + root_discriminant);


    if (minimum_solution > 0){
        return minimum_solution;
    }

    if (maximum_solution > 0){
        return maximum_solution;
    }
    return -1.0;

}


vec3 sampleSpherical(){
    float r1 = ((float) rand() / (RAND_MAX));
    float r2 = ((float) rand() / (RAND_MAX));
    float r3 = ((float) rand() / (RAND_MAX));
    vec3 sample = vec3(r1, r2, r3);
    sample = normalizeVector(sample);
    return sample;
}


#endif