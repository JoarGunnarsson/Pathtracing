#ifndef materials
#define materials
#include "vec3.h"
#include "colors.h"


class Material{
    public:
        vec3 diffuseColor;
        float diffuseCoefficient;
        float refractiveIndex;
        float attenuationCoefficient;
        vec3 absorptionColor;
        vec3 emmissionColor;
        float lightIntensity;
        Material() :
            diffuseColor(WHITE),
            diffuseCoefficient(0.0f),
            refractiveIndex(1.0f),
            attenuationCoefficient(0.0f),
            absorptionColor(WHITE),
            emmissionColor(WHITE),
            lightIntensity(0.0f) {}
        Material(vec3 _diffuseColor, float _diffuseCoefficient=1, float _refractiveIndex=1, float _attenuationCoefficient=0, vec3 _absorptionColor=WHITE,
        vec3 _emmissionColor=WHITE, float _lightIntensity=0){
            diffuseColor = _diffuseColor;
            diffuseCoefficient = _diffuseCoefficient;
            refractiveIndex = _refractiveIndex;
            attenuationCoefficient = _attenuationCoefficient;
            absorptionColor = _absorptionColor;
            emmissionColor = _emmissionColor;
            lightIntensity = _lightIntensity;
        }
};

#endif