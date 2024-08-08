#ifndef materials
#define materials
#include "vec3.h"
#include "colors.h"
#include "utils.h"
#include "constants.h"
#include "objects.h"

class Object;

Hit findClosestHit(const Ray& ray, Object** objects, const int size);

struct brdfData{
    vec3 outgoingVector;
    vec3 brdfMultiplier;
    bool specular = false;
};

class Material{
    public:
        vec3 albedo;
        double refractiveIndex;
        double attenuationCoefficient;
        vec3 absorptionAlbedo;
        vec3 emmissionColor;
        double lightIntensity;
        Material(){
            albedo = WHITE;
            refractiveIndex = 1;
            attenuationCoefficient = 0;
            absorptionAlbedo = WHITE;
            attenuationCoefficient = 0;
            emmissionColor = WHITE;
            lightIntensity = 0;
        }
        Material(vec3 _diffuseColor, double _diffuseCoefficient=0.8, double _refractiveIndex=1, double _attenuationCoefficient=0,
        vec3 _emmissionColor=WHITE, double _lightIntensity=0){
            albedo = _diffuseColor * _diffuseCoefficient;
            refractiveIndex = _refractiveIndex;
            attenuationCoefficient = _attenuationCoefficient;
            absorptionAlbedo = vec3(1,1,1) - albedo;
            emmissionColor = _emmissionColor;
            lightIntensity = _lightIntensity;
        }

    virtual vec3 eval(){
        throw VirtualMethodNotAllowedException("this is a pure virtual method and should not be called.");
        vec3 vec;
        return vec;
    }

    virtual brdfData sample(const Hit& hit, Object** objectPtrList, const int numberOfObjects){
        throw VirtualMethodNotAllowedException("this is a pure virtual method and should not be called.");
        brdfData data;
        return data;
    }

    vec3 getLightEmittance(){
        return emmissionColor * lightIntensity;
    }
};


class DiffuseMaterial : public Material{
    public:
        using Material::Material;

    vec3 eval() override{
        return albedo / M_PI;
    }

    brdfData sample(const Hit& hit, Object** objectPtrList, const int numberOfObjects) override{
        vec3 outgoingVector = sampleCosineHemisphere(hit.normalVector);
        vec3 brdfMultiplier = albedo;
        brdfData data;
        data.outgoingVector = outgoingVector;
        data.brdfMultiplier = brdfMultiplier;
        data.specular = false;
        return data;
    }
};


class ReflectiveMaterial : public Material{
    public:
        using Material::Material;

    vec3 eval() override{
        return BLACK;
    }

    brdfData sample(const Hit& hit, Object** objectPtrList, const int numberOfObjects) override{
        vec3 outgoingVector = reflectVector(hit.incomingVector, hit.normalVector);
        brdfData data;
        data.outgoingVector = outgoingVector;
        data.brdfMultiplier = albedo;
        data.specular = true;
        return data;
    }
};


class TransparentMaterial : public Material{
    public:
        using Material::Material;

    vec3 eval() override{
        return BLACK;
    }

    brdfData sample(const Hit& hit, Object** objectPtrList, const int numberOfObjects) override{
        double incomingDotNormal = dotVectors(hit.incomingVector, hit.normalVector);
        vec3 fresnelNormal;
        bool inside = incomingDotNormal > 0.0;
        double n1;
        double n2;
        if (!inside){
            fresnelNormal = -hit.normalVector;
            n1 = constants::airRefractiveIndex;
            n2 = refractiveIndex;
        }
        else{
            fresnelNormal = hit.normalVector;
            n1 = refractiveIndex;
            n2 = constants::airRefractiveIndex;
        }

        vec3 transmittedVector = refractVector(fresnelNormal, hit.incomingVector, n1, n2);

        double F_r = 1;
        if (transmittedVector.length_squared() != 0){
            F_r = fresnelMultiplier(hit.incomingVector, -fresnelNormal, n1, n2);
        }

        double randomNum = randomUniform(0, 1);
        bool isReflected = randomNum <= F_r;
        
        vec3 brdfMultiplier;
        vec3 outgoingVector;
        if (isReflected){
            brdfMultiplier = albedo;
            outgoingVector = reflectVector(hit.incomingVector, -fresnelNormal);
        }
        else{
            bool enableRandomTransmission = false;
            if (enableRandomTransmission){
                vec3 randomHemispherePoint = sampleHemisphere(fresnelNormal);
                double smoothness = 0.5;
                vec3 scaledTransmittedVector = transmittedVector * smoothness;
                vec3 scaledHemispherePoint = randomHemispherePoint * (1 - smoothness);
                transmittedVector = scaledTransmittedVector + scaledHemispherePoint;
                transmittedVector = normalizeVector(transmittedVector);
            }

            Ray transmissionRay;
            transmissionRay.directionVector = transmittedVector;
            transmissionRay.startingPosition = hit.intersectionPoint;
            Hit transmissionHit = findClosestHit(transmissionRay, objectPtrList, numberOfObjects);
            vec3 attenuationColor;
            double distance = transmissionHit.distance;
            if (distance > 0 && !inside){
                vec3 log_attenuation = absorptionAlbedo * attenuationCoefficient * (-distance);
                attenuationColor = expVector(log_attenuation);
            }
            else{
                attenuationColor = albedo;
            }

            double refractionIntensityFactor = pow(n2 / n1, 2);
            brdfMultiplier = attenuationColor * refractionIntensityFactor;
            outgoingVector = transmittedVector;
        }
        brdfData data;
        data.outgoingVector = outgoingVector;
        data.brdfMultiplier = brdfMultiplier;
        data.specular = true;
        return data;
    }
};


#endif