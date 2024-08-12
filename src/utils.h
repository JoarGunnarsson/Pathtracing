#ifndef utils
#define utils
#include "vec3.h"
#include <random>


class VirtualMethodNotAllowedException : public std::logic_error {
public:
    explicit VirtualMethodNotAllowedException(const std::string& message)
        : std::logic_error(message) {}
};


std::random_device rand_dev;
std::minstd_rand normal_generator(rand_dev());
std::normal_distribution<double> normal_distribution(0, 1);

std::minstd_rand uniform_generator(rand_dev());
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


enum reflectionType{
    DIFFUSE = 0,
    REFLECTED = 1,
    TRANSMITTED = 2
};


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
    int type = DIFFUSE;
};


inline double posFmod(const double a, const double b){
    return fmod((fmod(a, b) + b), b);
}


inline double clamp(const double value, const double min, const double max){
    return std::max(std::min(value, max), min);
}


double solveQuadratic(const double b, const double c){
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


vec3 sampleHemisphere(const vec3& normal){
    vec3 sample = sampleSpherical();
    if (dotVectors(normal, sample) < 0){
        return -sample;
    }
    return sample;
}

void setPerpendicularVectors(const vec3& zHat, vec3& xHat, vec3& yHat){
    vec3 nonParallelVector = vec3(1.0, 0.0, 0.0);
    if (std::abs(dotVectors(nonParallelVector, zHat)) == 1.0){
        nonParallelVector = vec3(0.0, 1.0, 0.0);
    }
    
    xHat = crossVectors(zHat, nonParallelVector);
    xHat = normalizeVector(xHat);
    yHat = crossVectors(zHat, xHat);
    yHat = normalizeVector(yHat);
}

vec3 sampleAngledHemisphere(const vec3& normalVector, const double cosMax){
    vec3 xHat;
    vec3 yHat;
    setPerpendicularVectors(normalVector, xHat, yHat);
    double phi = randomUniform(0, 2.0 * M_PI);
    double cosTheta = randomUniform(cosMax, 1);
    double sinTheta = sqrt(1 - (cosTheta * cosTheta));
    double x = sinTheta * cos(phi);
    double y = sinTheta * sin(phi); 
    double z = cosTheta;
    return xHat * x + yHat * y + normalVector * z;
}


vec3 sampleCosineHemisphere(const vec3& normalVector){
    vec3 xHat;
    vec3 yHat;
    setPerpendicularVectors(normalVector, xHat, yHat);

    double theta = ((double) rand() / (RAND_MAX)) * 2 * M_PI;
    double radius = sqrt((double) rand() / (RAND_MAX));
    double x = cos(theta) * radius;
    double y = sin(theta) * radius;
    double z = sqrt(1 - x * x - y * y);
    return xHat * x + yHat * y + normalVector * z;
}


vec3 reflectVector(const vec3& directionVector, const vec3& normalVector){
    return directionVector - normalVector * 2.0 * dotVectors(normalVector, directionVector);
}


vec3 refractVector(const vec3& normalVector, const vec3& incidentVector, const double n1, const double n2){
    double mu = n1 / n2;
    double cosIncident = dotVectors(normalVector, incidentVector);
    double lengthInNormalDirectionSquared = 1 - mu * mu * (1 - cosIncident * cosIncident);
    if (lengthInNormalDirectionSquared < 0){
        return vec3(0,0,0);
    }
    vec3 perpendicularVectors = incidentVector - normalVector * cosIncident;
    return normalVector * sqrt(lengthInNormalDirectionSquared) + perpendicularVectors * mu;
}


double fresnelDielectric(const double cosIncident, const double n1, const double n2){
    double sinIncident = sqrt(1 - cosIncident * cosIncident);
    double cosTransmitted = sqrt(1 - pow(n1 / n2 * sinIncident, 2));
    double n1CosIncident = n1 * cosIncident;
    double n2CosTransmitted = n2 * cosTransmitted;
    double n1CosTransmitted = n1 * cosTransmitted;
    double n2CosIncident = n2 * cosIncident;
    double R_s = pow((n1CosIncident - n2CosTransmitted) / (n1CosIncident + n2CosTransmitted), 2);
    double R_p = pow((n1CosTransmitted - n2CosIncident) / (n1CosTransmitted + n2CosIncident), 2);
    return 0.5 * (R_s + R_p);
}


double fresnelConductor(const double cosTheta, const double n1, const double k1, const double n2, const double k2){
    double eta;
    double k;
    if (k1 == 0){
        eta = n2 / n1;
        k = k2 / n1;
    }
    else{
        double complexIndexNorm = (n1 * n1 * k1 * k1);
        eta = n2 * n1 / complexIndexNorm;
        k = - k1 * n2 / complexIndexNorm;
    }

    double cosTheta2 = cosTheta * cosTheta;
    double sinTheta2 = 1 - cosTheta2;
    double f0 = sqrt(pow(eta * eta - k * k - sinTheta2, 2) + 4 * eta * eta * k * k);
    double a2b2 = f0;
    double a = sqrt(0.5 * f0 + eta * eta - k * k - sinTheta2);
    double f1 = a2b2 + cosTheta2;
    double f2 = 2 * a * cosTheta;
    double f3 = cosTheta2 * a2b2 + sinTheta2 * sinTheta2;
    double f4 = f2 * sinTheta2;

    double R_p = (f1 - f2) / (f1 + f2);
    double R_s = R_p * (f3 - f4) / (f3 + f4); 
    return 0.5 * (R_p + R_s);
}


double fresnelMultiplier(const double cosIncident, const double n1, const double k1, const double n2, const double k2, const bool isDielectric){
    if (isDielectric || (k1==0 && k2==0)){
        return fresnelDielectric(cosIncident, n1, n2);
    }

    return fresnelConductor(cosIncident, n1, k1, n2, k2);
}
#endif