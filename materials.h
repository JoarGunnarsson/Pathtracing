#ifndef materials
#define materials
#include "vec3.h"
#include "colors.h"
#include "utils.h"

struct brdfData{
    vec3 outgoingVector;
    vec3 brdfMultiplier;
};

class Material{
    public:
        vec3 albedo;
        double refractiveIndex;
        double attenuationCoefficient;
        vec3 absorptionColor;
        vec3 emmissionColor;
        double lightIntensity;
        Material() :
            albedo(WHITE),
            refractiveIndex(1.0f),
            attenuationCoefficient(0.0f),
            absorptionColor(WHITE),
            emmissionColor(WHITE),
            lightIntensity(0.0f) {}
        Material(vec3 _diffuseColor, double _diffuseCoefficient=1, double _refractiveIndex=1, double _attenuationCoefficient=0, vec3 _absorptionColor=WHITE,
        vec3 _emmissionColor=WHITE, double _lightIntensity=0){
            albedo = multiplyVector(_diffuseColor, _diffuseCoefficient);
            refractiveIndex = _refractiveIndex;
            attenuationCoefficient = _attenuationCoefficient;
            absorptionColor = _absorptionColor;
            emmissionColor = _emmissionColor;
            lightIntensity = _lightIntensity;
        }

    vec3 eval(){
        return divideVector(albedo, M_PI);
    }

    double pdf(){
        return 1 / (2 * M_PI);
    }

    brdfData sample(Hit& hit){
        vec3 outgoingVector = sampleCosineHemisphere(hit.normalVector);
        vec3 brdfMultiplier = albedo;
        brdfData data;
        data.outgoingVector = outgoingVector;
        data.brdfMultiplier = brdfMultiplier;
        return data;
    }
};

#endif