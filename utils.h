#ifndef utils
#define utils
#include "vec3.h"
#include <random>


std::default_random_engine normal_generator;
std::normal_distribution<double> normal_distribution(0, 1);


struct Hit{
    int intersectedObjectIndex;
    int objectID;
    double distance;
    vec3 intersectionPoint;
    vec3 incomingVector;
    vec3 normalVector;
};


struct Ray{
    vec3 startingPosition;
    vec3 directionVector;
    bool specular = false;
};


double solveQuadratic(double b, double c){
    double discriminant = pow(b, 2) - 4 * c;
    if (discriminant < 0){
        return -1.0;
    }
    double root_discriminant = sqrt(discriminant);
    double maximum_solution = - 1.0 / 2.0 * (b - root_discriminant);
    double minimum_solution = - 1.0 / 2.0 * (b + root_discriminant);

    if (minimum_solution > 0){
        return minimum_solution;
    }

    if (maximum_solution > 0){
        return maximum_solution;
    }
    return -1.0;

}


vec3 sampleSpherical(){
    double r1 = normal_distribution(normal_generator);
    double r2 = normal_distribution(normal_generator);
    double r3 = normal_distribution(normal_generator);
    vec3 sample = vec3(r1, r2, r3);
    sample = normalizeVector(sample);
    return sample;
}


vec3 sampleHemisphere(vec3 normal){
    vec3 sample = sampleSpherical();
    if (dotVectors(normal, sample) < 0){
        return -sample;
    }
    return sample;
}


vec3 sampleCosineHemisphere(vec3& normalVector){
    vec3 nonParallelVector = vec3(1.0, 0.0, 0.0);
    if (std::abs(dotVectors(nonParallelVector, normalVector)) == 1.0){
        nonParallelVector = vec3(0.0, 1.0, 0.0);
    }
    
    vec3 xHat = crossVectors(normalVector, nonParallelVector);
    xHat = normalizeVector(xHat);
    vec3 yHat = crossVectors(normalVector, xHat);
    yHat = normalizeVector(yHat);

    double theta = ((double) rand() / (RAND_MAX)) * 2 * M_PI;
    double radius = sqrt((double) rand() / (RAND_MAX));
    double x = cos(theta) * radius;
    double y = sin(theta) * radius;
    double z = sqrt(1 - pow(x, 2) - pow(y,2));
    vec3 scaledX = multiplyVector(xHat, x);
    vec3 scaledY = multiplyVector(yHat, y);
    vec3 scaledZ = multiplyVector(normalVector, z);
    vec3 xPlusY = addVectors(scaledX, scaledY);
    vec3 sphere_points = addVectors(xPlusY, scaledZ);
    return sphere_points;
}


#endif