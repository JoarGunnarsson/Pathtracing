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
        if (isnan(u) or isnan(v)){
            return 0;
        }
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
        if (isnan(u) or isnan(v)){
            return vec3(0,0,0);
        }
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

    ValueMap3D* albedoMap = whiteMap;
    double refractiveIndex = 1;
    double attenuationCoefficient = 0;
    vec3 absorptionAlbedo = WHITE;
    ValueMap3D* emmissionColorMap = whiteMap;
    ValueMap1D* lightIntensityMap = zeroMap;
    bool isDielectric = true;
    double imaginaryRefractiveIndex = 0;
    ValueMap1D* roughnessMap = zeroMap;
    ValueMap1D* percentageDiffuseMap = zeroMap;
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
        ValueMap1D* roughnessMap;
        ValueMap1D* percentageDiffuseMap;
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

            roughnessMap = data.roughnessMap;
            percentageDiffuseMap = data.percentageDiffuseMap;
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
        
        brdfData data;
        if (isReflected){
            data.type = REFLECTED;
            data.brdfMultiplier = isDielectric ? WHITE : albedoMap -> get(u, v);
            data.outgoingVector = reflectVector(hit.incomingVector, fresnelNormal);
        }
        else{
            data.type = TRANSMITTED;
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

            data.brdfMultiplier = attenuationColor;
            data.outgoingVector = transmittedVector;
        }

        return data;
    }
};


class GlossyMaterial : public DiffuseMaterial, public ReflectiveMaterial{
    public:
        using DiffuseMaterial::DiffuseMaterial;

    vec3 eval(const double u, const double v) override{
        double randomNum = randomUniform(0, 1);
        bool doSpecular = randomNum <= 1.0 - percentageDiffuseMap -> get(u, v);
        if (doSpecular){
            return ReflectiveMaterial::eval(u, v);
        }
        else{
            return DiffuseMaterial::eval(u, v);
        }
    }

    brdfData sample(const Hit& hit, Object** objectPtrList, const int numberOfObjects, const double u, const double v) override{
        double randomNum = randomUniform(0, 1);
        bool doSpecular = randomNum <= 1.0 - percentageDiffuseMap -> get(u, v);
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
        using DiffuseMaterial::DiffuseMaterial;

    vec3 eval(const double u, const double v) override{
        double randomNum = randomUniform(0, 1);
        bool doSpecular = randomNum <= 1.0 - percentageDiffuseMap -> get(u, v);
        if (doSpecular){
            return TransparentMaterial::eval(u, v);
        }
        else{
            return DiffuseMaterial::eval(u, v);
        }
    }

    brdfData sample(const Hit& hit, Object** objectPtrList, const int numberOfObjects, const double u, const double v) override{
        double randomNum = randomUniform(0, 1);
        bool doSpecular = randomNum <= 1.0 - percentageDiffuseMap -> get(u, v);
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
            data.outgoingVector = normalizeVector(data.outgoingVector * (1.0 - roughness) + randomRay * roughness);
            return data;
        }
        else{
            return DiffuseMaterial::sample(hit, objectPtrList, numberOfObjects, u, v);
        }
    }
};


class MicrofacetMaterial : public Material{
    public:
        using Material::Material;

    vec3 eval(const double u, const double v) override{
        return BLACK;
    }

    double chi(const double x){
        return x > 0 ? 1 : 0;
    }

    double G1(const vec3& halfVector, const vec3& normalVector, const vec3& v, const double alpha){
        double cosTheta = dotVectors(halfVector, v);
        double tanTheta = sqrt((1.0 - cosTheta * cosTheta) / (cosTheta * cosTheta));
        double a = 1.0 / (alpha * tanTheta);

        return chi(cosTheta / dotVectors(v, normalVector)) * 2.0 / (1 + std::erf(a) + 1.0 / (sqrt(M_PI) * a) * std::exp(-a*a));
    }

    double G(const vec3& halfVector, const vec3& normalVector, const vec3& incidentVector, const vec3& outgoingVector, const double alpha){
        return G1(halfVector, normalVector, incidentVector, alpha) * G1(halfVector, normalVector, outgoingVector, alpha);
    }

    vec3 specularSample(const vec3& normalVector, const double alpha){
        double r1 = randomUniform(0, 1);
        double r2 = randomUniform(0, 1);
        double phi = 2 * M_PI * r2;
        double tanTheta2 = - alpha * alpha * std::log(1 - r1);
        double cosTheta2 = 1.0 / (1.0 + tanTheta2);
        
        double cosTheta = sqrt(cosTheta2);
        double sinTheta = sqrt(1 - cosTheta2);

        vec3 xHat;
        vec3 yHat;
        setPerpendicularVectors(normalVector, xHat, yHat);
        
        return xHat * sinTheta * cos(phi) + yHat * sinTheta * sin(phi) + normalVector * cosTheta;
    }

    brdfData sampleDiffuse(const vec3& normalVector, const double u, const double v){
        brdfData data;
        data.brdfMultiplier = albedoMap -> get(u, v);
        data.outgoingVector = sampleCosineHemisphere(normalVector);
        data.type = REFLECTED;
        return data;
    }

    brdfData sampleReflection(const vec3& sampledHalfVector, const vec3& normalVector, const vec3& incidentVector, const double cosineFactor, const double u, const double v, const double alpha){
        brdfData data;
        vec3 reflectionColor = isDielectric ? WHITE : albedoMap -> get(u, v);
        data.outgoingVector = sampledHalfVector * 2.0 * dotVectors(incidentVector, sampledHalfVector) - incidentVector;

        data.brdfMultiplier = reflectionColor * G(sampledHalfVector, normalVector, incidentVector, data.outgoingVector, alpha) * cosineFactor;
        data.type = REFLECTED;
        return data;
    }

    brdfData sampleTransmission(const vec3& sampledHalfVector, const vec3& normalVector, const vec3& incidentVector, const double cosineFactor, const double eta, const double alpha){
        //TODO: Add attenuation.
        brdfData data;
        vec3 attenuatedColor = vec3(1,1,1);
        double c = dotVectors(incidentVector, sampledHalfVector);
        double hFactor = eta * c - sign(dotVectors(incidentVector, normalVector)) * sqrt(1 + eta*eta * (c*c-1));
        data.outgoingVector = sampledHalfVector * hFactor - incidentVector * eta;
        data.brdfMultiplier = attenuatedColor * G(sampledHalfVector, normalVector, incidentVector, data.outgoingVector, alpha) * cosineFactor;
        data.type = TRANSMITTED;
        return data;
    }
    
    brdfData sample(const Hit& hit, Object** objectPtrList, const int numberOfObjects, const double u, const double v) override{
        double incomingDotNormal = dotVectors(hit.incomingVector, hit.normalVector);
        vec3 fresnelNormal;
        bool outside = incomingDotNormal <= 0.0;
        double n1;
        double k1;
        double n2;
        double k2;
        if (outside){
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

        double eta = n1 / n2;

        double alpha = roughnessMap -> get(u, v);

        vec3 sampledHalfVector = specularSample(-fresnelNormal, alpha);
                
        double iDotH = dotVectors(hit.incomingVector, sampledHalfVector);
        double iDotN = dotVectors(hit.incomingVector, fresnelNormal);
        double nDotH = dotVectors(sampledHalfVector, fresnelNormal);
        double cosineFactor = std::abs(iDotH / (iDotN * nDotH));

        double F_r = fresnelMultiplier(-iDotH, n1, k1, n2, k2, isDielectric);
        double randomNum = randomUniform(0, 1);
        bool reflectSpecular = randomNum < F_r || 1 + eta*eta * (iDotH*iDotH-1) < 0;
        
        if (reflectSpecular){
            return sampleReflection(sampledHalfVector, -fresnelNormal, -hit.incomingVector, cosineFactor, u, v, alpha);
        }
        else{
            double randomNum2 = randomUniform(0, 1);
            bool diffuse = randomNum2 < percentageDiffuseMap -> get(u, v);
            return diffuse ? sampleDiffuse(-fresnelNormal, u, v) : sampleTransmission(sampledHalfVector, -fresnelNormal, -hit.incomingVector, cosineFactor, eta, alpha);
        }
    }
};

#endif