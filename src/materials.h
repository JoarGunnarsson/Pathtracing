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
    int type = DIFFUSE;
};


class Material{
    public:
        vec3 albedo;
        double refractiveIndex;
        double attenuationCoefficient;
        vec3 absorptionAlbedo;
        vec3 emmissionColor;
        double lightIntensity;
        bool isDielectric;
        double imaginaryRefractiveIndex;
        Material(){
            albedo = WHITE;
            refractiveIndex = 1;
            attenuationCoefficient = 0;
            absorptionAlbedo = WHITE;
            attenuationCoefficient = 0;
            emmissionColor = WHITE;
            lightIntensity = 0;
            isDielectric = true;
            imaginaryRefractiveIndex = 0;
        }
        Material(vec3 _diffuseColor, double _diffuseCoefficient=0.8, double _refractiveIndex=1, double _attenuationCoefficient=0,
        vec3 _emmissionColor=WHITE, double _lightIntensity=0, bool _isDielectric=true, double _imaginaryRefractiveIndex=0){
            albedo = _diffuseColor * _diffuseCoefficient;
            refractiveIndex = _refractiveIndex;
            attenuationCoefficient = _attenuationCoefficient;
            absorptionAlbedo = vec3(1,1,1) - albedo;
            emmissionColor = _emmissionColor;
            lightIntensity = _lightIntensity;
            isDielectric = _isDielectric;
            if (isDielectric){
                imaginaryRefractiveIndex = 0;
            }
            else{
                imaginaryRefractiveIndex = _imaginaryRefractiveIndex;
            }
        }

    virtual vec3 eval(const vec3& incidentVector, const vec3& outgoingVector, const vec3& normalVector){
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


class DiffuseMaterial : public virtual Material{
    public:
        using Material::Material;

    vec3 eval(const vec3& incidentVector, const vec3& outgoingVector, const vec3& normalVector) override{
        return albedo / M_PI;
    }

    brdfData sample(const Hit& hit, Object** objectPtrList, const int numberOfObjects) override{
        vec3 adjustedNormal = dotVectors(hit.incomingVector, hit.normalVector) < 0 ? hit.normalVector : -hit.normalVector;
        vec3 outgoingVector = sampleCosineHemisphere(adjustedNormal);
        vec3 brdfMultiplier = albedo;
        brdfData data;
        data.outgoingVector = outgoingVector;
        data.brdfMultiplier = brdfMultiplier;
        data.type = DIFFUSE;
        return data;
    }
};


class ReflectiveMaterial : public virtual Material{
    public:
        using Material::Material;

    vec3 eval(const vec3& incidentVector, const vec3& outgoingVector, const vec3& normalVector) override{
        return BLACK;
    }

    brdfData sample(const Hit& hit, Object** objectPtrList, const int numberOfObjects) override{
        vec3 adjustedNormal = dotVectors(hit.incomingVector, hit.normalVector) < 0 ? hit.normalVector : -hit.normalVector;
        vec3 outgoingVector = reflectVector(hit.incomingVector, adjustedNormal);
        brdfData data;
        data.outgoingVector = outgoingVector;
        data.brdfMultiplier = isDielectric ? WHITE : albedo;
        data.type = REFLECTED;
        return data;
    }
};


class TransparentMaterial : public virtual Material{
    public:
        using Material::Material;

    vec3 eval(const vec3& incidentVector, const vec3& outgoingVector, const vec3& normalVector) override{
        return BLACK;
    }

    brdfData sample(const Hit& hit, Object** objectPtrList, const int numberOfObjects) override{
        double incomingDotNormal = dotVectors(hit.incomingVector, hit.normalVector);
        vec3 fresnelNormal;
        bool inside = incomingDotNormal > 0.0;
        double n1;
        double k1;
        double n2;
        double k2;
        if (!inside){
            fresnelNormal = -hit.normalVector;
            n1 = constants::airRefractiveIndex;
            k1 = 0;
            n2 = refractiveIndex;
            k2 = imaginaryRefractiveIndex;
        }
        else{
            fresnelNormal = hit.normalVector;
            n1 = refractiveIndex;
            k1 = imaginaryRefractiveIndex;
            n2 = constants::airRefractiveIndex;
            k2 = 0;
        }

        vec3 transmittedVector = refractVector(fresnelNormal, hit.incomingVector, n1, n2);

        double F_r = 1;
        if (transmittedVector.length_squared() != 0){
            double cosIncident = dotVectors(hit.incomingVector, fresnelNormal);
            F_r = fresnelMultiplier(cosIncident, n1, k1, n2, k2, isDielectric);
        }

        double randomNum = randomUniform(0, 1);
        bool isReflected = randomNum <= F_r;
        
        vec3 brdfMultiplier;
        vec3 outgoingVector;
        if (isReflected && isDielectric){
            brdfMultiplier = WHITE;
            outgoingVector = reflectVector(hit.incomingVector, -fresnelNormal);
        }
        else if (isReflected && !isDielectric){
            brdfMultiplier = albedo;
            outgoingVector = reflectVector(hit.incomingVector, -fresnelNormal);
        }
        else{
            Ray transmissionRay;
            transmissionRay.directionVector = transmittedVector;
            transmissionRay.startingPosition = hit.intersectionPoint;
            Hit transmissionHit = findClosestHit(transmissionRay, objectPtrList, numberOfObjects);
            vec3 attenuationColor;
            double distance = transmissionHit.distance;
            if (distance > 0 && !inside){
                vec3 log_attenuation = absorptionAlbedo * attenuationCoefficient * (-distance);
                attenuationColor = albedo * expVector(log_attenuation);
            }
            else{
                attenuationColor = albedo;
            }

            double refractionIntensityFactor = n2 * n2 / (n1 * n1);
            brdfMultiplier = attenuationColor * refractionIntensityFactor;
            outgoingVector = transmittedVector;
        }

        brdfData data;
        data.outgoingVector = outgoingVector;
        data.brdfMultiplier = brdfMultiplier;
        data.type = TRANSMITTED;
        return data;
    }
};


class GlossyMaterial : public DiffuseMaterial, public ReflectiveMaterial{
    public:
        double roughness = 0.9;
        double percentageSpecular = 0.3;
        using DiffuseMaterial::DiffuseMaterial;

    vec3 eval(const vec3& incidentVector, const vec3& outgoingVector, const vec3& normalVector) override{
        double randomNum = randomUniform(0, 1);
        bool doSpecular = randomNum <= percentageSpecular;
        if (doSpecular){
            return ReflectiveMaterial::eval(incidentVector, outgoingVector, normalVector);
        }
        else{
            return DiffuseMaterial::eval(incidentVector, outgoingVector, normalVector);
        }
    }

    brdfData sample(const Hit& hit, Object** objectPtrList, const int numberOfObjects) override{
        double randomNum = randomUniform(0, 1);
        bool doSpecular = randomNum <= percentageSpecular;
        if (doSpecular){
            brdfData data = ReflectiveMaterial::sample(hit, objectPtrList, numberOfObjects);
            vec3 randomRay = sampleCosineHemisphere(hit.normalVector);
            data.outgoingVector = normalizeVector(data.outgoingVector * (1 - roughness) + randomRay * roughness);
            return data;
        }
        else{
            return DiffuseMaterial::sample(hit, objectPtrList, numberOfObjects);
        }
    }
};


class FrostyMaterial : public DiffuseMaterial, public TransparentMaterial{
    public:
        double roughness = 0.05;
        double percentageSpecular = 1;
        using DiffuseMaterial::DiffuseMaterial;

    vec3 eval(const vec3& incidentVector, const vec3& outgoingVector, const vec3& normalVector) override{
        double randomNum = randomUniform(0, 1);
        bool doSpecular = randomNum <= percentageSpecular;
        if (doSpecular){
            return TransparentMaterial::eval(incidentVector, outgoingVector, normalVector);
        }
        else{
            return DiffuseMaterial::eval(incidentVector, outgoingVector, normalVector);
        }
    }

    brdfData sample(const Hit& hit, Object** objectPtrList, const int numberOfObjects) override{
        double randomNum = randomUniform(0, 1);
        bool doSpecular = randomNum <= percentageSpecular;
        if (doSpecular){
            brdfData data = TransparentMaterial::sample(hit, objectPtrList, numberOfObjects);
            vec3 adjustedNormal;
            if (data.type == TRANSMITTED){
                adjustedNormal = -hit.normalVector;
            }   
            else{
                adjustedNormal = hit.normalVector;
            }        
            vec3 randomRay = sampleCosineHemisphere(adjustedNormal);
            data.outgoingVector = normalizeVector(data.outgoingVector * (1 - roughness) + randomRay * roughness);
            return data;
        }
        else{
            return DiffuseMaterial::sample(hit, objectPtrList, numberOfObjects);
        }
    }
};
#endif