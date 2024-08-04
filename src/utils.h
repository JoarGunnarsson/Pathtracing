#ifndef utils
#define utils
#include "vec3.h"
#include <random>


class VirtualMethodNotAllowedException : public std::logic_error {
public:
    explicit VirtualMethodNotAllowedException(const std::string& message)
        : std::logic_error(message) {}
};


std::minstd_rand normal_generator;
std::normal_distribution<double> normal_distribution(0, 1);

std::minstd_rand uniform_generator;
std::uniform_real_distribution<double> uniform_dist(0, 1);


inline double randomUniform(const double low, const double high){
    return (high - low) * uniform_dist(uniform_generator) + low;
}


inline int randomInt(const int low, const int high){
    return (int) randomUniform(low, high);
}


inline double randomNormal(){
    return normal_distribution(normal_generator);
}


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
    double r1 = randomNormal();
    double r2 = randomNormal();
    double r3 = randomNormal();
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
    return xHat * x + yHat * y + normalVector * z;;
}


vec3 reflectVector(const vec3& directionVector, const vec3& normalVector){
    return directionVector - normalVector * 2.0 * dotVectors(normalVector, directionVector);;
}


vec3 refractVector(const vec3& normalVector, const vec3& incidentVector, const double n1, const double n2){
    double mu = n1 / n2;
    double cosIncident = dotVectors(normalVector, incidentVector);
    double lengthInNormalDirectionSquared = 1 - pow(mu, 2) * (1 - pow(cosIncident, 2));
    if (lengthInNormalDirectionSquared < 0){
        return vec3(0,0,0);
    }
    vec3 perpendicularVectors = incidentVector - normalVector * cosIncident;
    return normalVector * sqrt(lengthInNormalDirectionSquared) + perpendicularVectors * mu;
}


double schlickApproximation(const double R0, const double cosTheta){
    double R = R0 + (1 - R0) * pow((1 - cosTheta), 5);
    return R;
}


double fresnelMultiplier(const vec3& incidentVector, const vec3& transmittedVector, const vec3& normalVector, const double n1, const double n2){
    double R0 = pow(((n1 - n2) / (n1 + n2)), 2);
    if (n2 >= n1){
        double cosIncident = dotVectors(incidentVector, normalVector);
        return schlickApproximation(R0, cosIncident);
    }
    else{
        double cosTransmitted = dotVectors(transmittedVector, normalVector);
        return schlickApproximation(R0, cosTransmitted);
    }
}

#endif