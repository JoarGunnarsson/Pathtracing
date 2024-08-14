#ifndef materials
#define materials
#include "vec3.h"
#include "colors.h"
#include "utils.h"
#include "constants.h"
#include "objects.h" // TODO: is this needed?


class Object;


Hit findClosestHit(const Ray& ray, Object** objects, const int size);


struct brdfData{
    vec3 outgoingVector;
    vec3 brdfMultiplier;
    int type = DIFFUSE;
};


class ValueMap{
    public:
        double* data;
        int width;
        double uMax;
        int height;
        double vMax;
        ValueMap(){
            data = new double(0);
            width = 1;
            height = 1;
            uMax = 1;
            vMax = 1;
        }
        ValueMap(double* _data, const int _width=1, const int _height=1, const double _uMax=1, const double _vMax=1){
            data = _data;
            width = _width;
            height = _height;
            uMax = _uMax;
            vMax = _vMax;
        }
        ~ValueMap(){
            delete[] data;
        }
};


class ValueMap1D : public ValueMap{
    public:
        using ValueMap::ValueMap;

    double get(const double u, const double v) {
        int uIdx = int((double) width * posFmod(u / uMax, 1.0));
        int vIdx = int((double) height * posFmod(v / vMax, 1.0));
        int index = (vIdx * width + uIdx);
        return data[index];
    }
};


class ValueMap3D : public ValueMap{
    public:
        using ValueMap::ValueMap;

    vec3 get(const double u, const double v){
        int uIdx = int((double) width * posFmod(u / uMax, 1.0));
        int vIdx = int((double) height * posFmod(v / vMax, 1.0));
        int startIndex = 3 * (vIdx * width + uIdx);
        return vec3(data[startIndex], data[startIndex + 1], data[startIndex + 2]);
    }
};


struct MaterialData{
    ValueMap3D* whiteMap = new ValueMap3D(WHITE);
    double* zero = new double(0);
    ValueMap1D* zeroMap = new ValueMap1D(zero);
    double* ohFive = new double(0.5);
    ValueMap1D* ohFiveMap = new ValueMap1D(ohFive);
    double* one = new double(1);
    ValueMap1D* oneMap = new ValueMap1D(one);

    ValueMap3D* albedoMap = whiteMap;
    double refractiveIndex = 1;
    double attenuationCoefficient = 0;
    vec3 absorptionAlbedo = WHITE;
    ValueMap3D* emmissionColorMap = whiteMap;
    ValueMap1D* lightIntensityMap = zeroMap;
    bool isDielectric = true;
    double imaginaryRefractiveIndex = 0;
    ValueMap1D* roughnessMap = ohFiveMap;
    ValueMap1D* percentageSpecularMap = oneMap;
    bool isLightSource = false;
};

class Material{
    public:
        ValueMap3D* albedoMap;
        double refractiveIndex;
        double attenuationCoefficient;
        vec3 absorptionAlbedo;
        ValueMap3D* emmissionColorMap;
        ValueMap1D* lightIntensityMap;
        bool isDielectric;
        double imaginaryRefractiveIndex;
        bool isLightSource;
        Material(){}
        Material(MaterialData data){
            albedoMap = data.albedoMap;
            refractiveIndex = data.refractiveIndex;
            attenuationCoefficient = data.attenuationCoefficient;
            absorptionAlbedo = data.absorptionAlbedo;
            emmissionColorMap = data.emmissionColorMap;
            lightIntensityMap = data.lightIntensityMap;
            isDielectric = data.isDielectric;
            if (isDielectric){
                imaginaryRefractiveIndex = 0;
            }
            else{
                imaginaryRefractiveIndex = data.imaginaryRefractiveIndex;
            }
            isLightSource = data.isLightSource;
        }

    virtual vec3 eval(const double u, const double v){
        throw VirtualMethodNotAllowedException("this is a pure virtual method and should not be called.");
        vec3 vec;
        return vec;
    }

    virtual brdfData sample(const Hit& hit, Object** objectPtrList, const int numberOfObjects, const double u, const double v){
        throw VirtualMethodNotAllowedException("this is a pure virtual method and should not be called.");
        brdfData data;
        return data;
    }

    vec3 getLightEmittance(const double u, const double v){
        return emmissionColorMap -> get(u, v) * lightIntensityMap -> get(u, v);
    }
};


class DiffuseMaterial : public virtual Material{
    public:
        using Material::Material;

    vec3 eval(const double u, const double v) override{
        return albedoMap -> get(u, v) / M_PI;
    }

    brdfData sample(const Hit& hit, Object** objectPtrList, const int numberOfObjects, const double u, const double v) override{
        vec3 adjustedNormal = dotVectors(hit.incomingVector, hit.normalVector) < 0 ? hit.normalVector : -hit.normalVector;
        vec3 outgoingVector = sampleCosineHemisphere(adjustedNormal);
        vec3 brdfMultiplier = albedoMap -> get(u, v);
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

    vec3 eval(const double u, const double v) override{
        return BLACK;
    }

    brdfData sample(const Hit& hit, Object** objectPtrList, const int numberOfObjects, const double u, const double v) override{
        vec3 adjustedNormal = dotVectors(hit.incomingVector, hit.normalVector) < 0 ? hit.normalVector : -hit.normalVector;
        vec3 outgoingVector = reflectVector(hit.incomingVector, adjustedNormal);
        brdfData data;
        data.outgoingVector = outgoingVector;
        data.brdfMultiplier = isDielectric ? WHITE : albedoMap -> get(u, v);
        data.type = REFLECTED;
        return data;
    }
};


class TransparentMaterial : public virtual Material{
    public:
        using Material::Material;

    vec3 eval(const double u, const double v) override{
        return BLACK;
    }

    brdfData sample(const Hit& hit, Object** objectPtrList, const int numberOfObjects, const double u, const double v) override{
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
        if (isReflected){
            brdfMultiplier = isDielectric ? WHITE : albedoMap -> get(u, v);
            outgoingVector = reflectVector(hit.incomingVector, fresnelNormal);
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
                attenuationColor = expVector(log_attenuation);
            }
            else{
                attenuationColor = WHITE;
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
        ValueMap1D* roughnessMap;
        ValueMap1D* percentageSpecularMap;
        GlossyMaterial(MaterialData data)
        : Material(data){
            roughnessMap = data.roughnessMap;
            percentageSpecularMap = data.percentageSpecularMap;
        }
        ~GlossyMaterial(){
            delete roughnessMap;
            delete percentageSpecularMap;
        }

    vec3 eval(const double u, const double v) override{
        double randomNum = randomUniform(0, 1);
        bool doSpecular = randomNum <= percentageSpecularMap -> get(u, v);
        if (doSpecular){
            return ReflectiveMaterial::eval(u, v);
        }
        else{
            return DiffuseMaterial::eval(u, v);
        }
    }

    brdfData sample(const Hit& hit, Object** objectPtrList, const int numberOfObjects, const double u, const double v) override{
        double randomNum = randomUniform(0, 1);
        bool doSpecular = randomNum <= percentageSpecularMap -> get(u, v);
        if (doSpecular){
            brdfData data = ReflectiveMaterial::sample(hit, objectPtrList, numberOfObjects, u, v);
            vec3 randomRay = sampleCosineHemisphere(hit.normalVector);
            double roughness = roughnessMap -> get(u, v);
            data.outgoingVector = normalizeVector(data.outgoingVector * (1 - roughness) + randomRay * roughness);
            return data;
        }
        else{
            return DiffuseMaterial::sample(hit, objectPtrList, numberOfObjects, u, v);
        }
    }
};


class FrostyMaterial : public DiffuseMaterial, public TransparentMaterial{
    public:
        ValueMap1D* roughnessMap;
        ValueMap1D* percentageSpecularMap;
        FrostyMaterial(MaterialData data)
        : Material(data){
            roughnessMap = data.roughnessMap;
            percentageSpecularMap = data.percentageSpecularMap;
        }

        ~FrostyMaterial(){
            delete roughnessMap;
            delete percentageSpecularMap;
        }

    vec3 eval(const double u, const double v) override{
        double randomNum = randomUniform(0, 1);
        bool doSpecular = randomNum <= percentageSpecularMap -> get(u, v);
        if (doSpecular){
            return TransparentMaterial::eval(u, v);
        }
        else{
            return DiffuseMaterial::eval(u, v);
        }
    }

    brdfData sample(const Hit& hit, Object** objectPtrList, const int numberOfObjects, const double u, const double v) override{
        double randomNum = randomUniform(0, 1);
        bool doSpecular = randomNum <= percentageSpecularMap -> get(u, v);
        if (doSpecular){
            brdfData data = TransparentMaterial::sample(hit, objectPtrList, numberOfObjects, u, v);
            vec3 adjustedNormal;
            if (data.type == TRANSMITTED){
                adjustedNormal = -hit.normalVector;
            }   
            else{
                adjustedNormal = hit.normalVector;
            }        
            vec3 randomRay = sampleCosineHemisphere(adjustedNormal);
            double roughness = roughnessMap -> get(u, v);
            data.outgoingVector = normalizeVector(data.outgoingVector * (1 - roughness) + randomRay * roughness);
            return data;
        }
        else{
            return DiffuseMaterial::sample(hit, objectPtrList, numberOfObjects, u, v);
        }
    }
};
#endif