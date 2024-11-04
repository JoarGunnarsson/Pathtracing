#ifndef materials
#define materials
#include "vec3.h"
#include "colors.h"
#include "utils.h"
#include "constants.h"
#include "valuemap.h"


class Object;


Hit findClosestHit(const Ray& ray, Object** objects, const int size);


struct brdfData{
    vec3 outgoingVector;
    vec3 brdfMultiplier;
    int type = DIFFUSE;
};


struct microfacetSampleArgs{
    vec3 sampledHalfVector;
    vec3 normalVector;
    vec3 incidentVector;
    vec3 intersectionPoint;
    double cosineFactor;
    double eta;
    double u;
    double v;
    double alpha;
    Object** objects;
    int numberOfObjects;
    bool outside;
};


struct microfacetData{
    bool outside;
    double alpha;
    double eta;
    vec3 normalIntoInterface;
    vec3 halfVector;
    double F_r;
};


struct MaterialData{
    ValueMap3D* albedoMap = nullptr;
    double refractiveIndex = 1;
    double extinctionCoefficient = 0;
    double attenuationCoefficient = 0;
    vec3 absorptionAlbedo = WHITE;
    ValueMap3D* emmissionColorMap = nullptr;
    ValueMap1D* lightIntensityMap = nullptr;
    bool isDielectric = true;
    ValueMap1D* roughnessMap = nullptr; 
    ValueMap1D* percentageDiffuseMap = nullptr;
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
        double extinctionCoefficient;
        bool isLightSource;
        ValueMap1D* roughnessMap;
        ValueMap1D* percentageDiffuseMap;
        Material(){}
        Material(MaterialData data){
            if (!data.albedoMap){
                data.albedoMap = new ValueMap3D(WHITE);
            }

            if (!data.emmissionColorMap){
                data.emmissionColorMap = new ValueMap3D(WHITE);
            }

            if (!data.lightIntensityMap){
                data.lightIntensityMap = new ValueMap1D(0);
            }

            if (!data.roughnessMap){
                data.roughnessMap = new ValueMap1D(0);
            }

            if (!data.percentageDiffuseMap){
                data.percentageDiffuseMap = new ValueMap1D(1);
            }

            albedoMap = data.albedoMap;
            refractiveIndex = data.refractiveIndex;
            attenuationCoefficient = data.attenuationCoefficient;
            absorptionAlbedo = data.absorptionAlbedo;
            emmissionColorMap = data.emmissionColorMap;
            lightIntensityMap = data.lightIntensityMap;
            isDielectric = data.isDielectric;
            if (isDielectric){
                extinctionCoefficient = 0;
            }
            else{
                extinctionCoefficient = data.extinctionCoefficient;
            }
            isLightSource = data.isLightSource;

            roughnessMap = data.roughnessMap;
            percentageDiffuseMap = data.percentageDiffuseMap;
        }

    void prepareDeletion(PointerManager* pm){
        pm -> addPointer(albedoMap);
        pm -> addPointer(emmissionColorMap);
        pm -> addPointer(lightIntensityMap);
        pm -> addPointer(roughnessMap);
        pm -> addPointer(percentageDiffuseMap);
        pm -> addPointer(this);
    }

    virtual vec3 eval(const Hit& hit, const double u, const double v){
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


class DiffuseMaterial : public Material{
    public:
        using Material::Material;

    vec3 eval(const Hit& hit, const double u, const double v) override{
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


class ReflectiveMaterial : public Material{
    public:
        using Material::Material;

    vec3 eval(const Hit& hit, const double u, const double v) override{
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


class TransparentMaterial : public Material{
    public:
        using Material::Material;

    vec3 eval(const Hit& hit, const double u, const double v) override{
        return BLACK;
    }

    brdfData sample(const Hit& hit, Object** objectPtrList, const int numberOfObjects, const double u, const double v) override{
        double incomingDotNormal = dotVectors(hit.incomingVector, hit.normalVector);
        vec3 normalIntoInterface;
        bool inside = incomingDotNormal > 0.0;
        double n1;
        double k1;
        double n2;
        double k2;
        if (!inside){
            normalIntoInterface = -hit.normalVector;
            n1 = constants::airRefractiveIndex;
            k1 = 0;
            n2 = refractiveIndex;
            k2 = extinctionCoefficient;
        }
        else{
            normalIntoInterface = hit.normalVector;
            n1 = refractiveIndex;
            k1 = extinctionCoefficient;
            n2 = constants::airRefractiveIndex;
            k2 = 0;
        }

        vec3 transmittedVector = refractVector(hit.incomingVector, normalIntoInterface, n1 / n2);

        double F_r = 1;
        if (transmittedVector.length_squared() != 0){
            double cosIncident = dotVectors(hit.incomingVector, normalIntoInterface);
            F_r = fresnelMultiplier(cosIncident, n1, k1, n2, k2, isDielectric);
        }

        double randomNum = randomUniform(0, 1);
        bool isReflected = randomNum <= F_r;
        
        brdfData data;
        if (isReflected){
            data.type = REFLECTED;
            data.brdfMultiplier = isDielectric ? WHITE : albedoMap -> get(u, v);
            data.outgoingVector = reflectVector(hit.incomingVector, normalIntoInterface);
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


class MicrofacetMaterial : public Material{
    public:
        using Material::Material;

    double chi(const double x){
        return x > 0 ? 1 : 0;
    }

    double G1(const vec3& halfVector, const vec3& normalVector, const vec3& v, const double alpha){
        double cosTheta = dotVectors(halfVector, v);
        double tanTheta = sqrt((1.0 - cosTheta * cosTheta) / (cosTheta * cosTheta));
        double a = 1.0 / (alpha * tanTheta);
        
        return chi(cosTheta / dotVectors(v, normalVector)) * a < 1.6 ? (3.535 * a  + 2.181 * a * a) / (1.0 + 2.276 * a + 2.577 * a * a) : 1.0;
    }

    double G(const vec3& halfVector, const vec3& normalVector, const vec3& incidentVector, const vec3& outgoingVector, const double alpha){
        return G1(halfVector, normalVector, incidentVector, alpha) * G1(halfVector, normalVector, outgoingVector, alpha);
    }

    microfacetData prepareMicrofacetData(const Hit& hit, const double u, const double v){
        double n1;
        double k1;
        double n2;
        double k2;
        vec3 normalIntoInterface;
        double incomingDotNormal = dotVectors(hit.incomingVector, hit.normalVector);
        bool outside = incomingDotNormal <= 0.0;
        if (outside){
            normalIntoInterface = -hit.normalVector;
            n1 = constants::airRefractiveIndex;
            k1 = 0;
            n2 = refractiveIndex;
            k2 = extinctionCoefficient;
        }
        else{
            normalIntoInterface = hit.normalVector;
            n1 = refractiveIndex;
            k1 = extinctionCoefficient;
            n2 = constants::airRefractiveIndex;
            k2 = 0;
        }

        double alpha = std::max(roughnessMap -> get(u, v), constants::EPSILON);

        vec3 sampledHalfVector = specularSample(-normalIntoInterface, alpha);
                
        double iDotH = dotVectors(hit.incomingVector, sampledHalfVector);

        double F_r = fresnelMultiplier(-iDotH, n1, k1, n2, k2, isDielectric);

        microfacetData data;
        data.outside = outside;
        data.alpha = alpha;
        data.eta = n1 / n2;
        data.normalIntoInterface = normalIntoInterface;
        data.halfVector = sampledHalfVector;
        data.F_r = F_r;
        return data;
    }

    vec3 eval(const Hit& hit, const double u, const double v) override{
        microfacetData data = prepareMicrofacetData(hit, u, v);
        
        double randomNum = randomUniform(0, 1);
        bool reflectSpecular = randomNum < data.F_r;
        if (reflectSpecular){
            return BLACK;
        }
        double transmit = randomUniform(0, 1) > percentageDiffuseMap -> get(u, v);

        if (transmit){
            return BLACK;
        }

        return albedoMap -> get(u, v) / M_PI;
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

    brdfData sampleDiffuse(const microfacetSampleArgs& args){
        brdfData data;
        data.brdfMultiplier = albedoMap -> get(args.u, args.v);
        data.outgoingVector = sampleCosineHemisphere(args.normalVector);
        data.type = DIFFUSE;
        return data;
    }

    brdfData sampleReflection(const microfacetSampleArgs& args){
        brdfData data;
        vec3 reflectionColor = isDielectric ? WHITE : albedoMap -> get(args.u, args.v);
        data.outgoingVector = reflectVector(-args.incidentVector, args.sampledHalfVector);
        data.brdfMultiplier = reflectionColor * G(args.sampledHalfVector, args.normalVector, args.incidentVector, data.outgoingVector, args.alpha) * args.cosineFactor;
        data.type = REFLECTED;
        return data;
    }

    vec3 computeAttenuatedColor(const microfacetSampleArgs& args, const vec3& outgoingVector){
        Ray transmissionRay;
        transmissionRay.directionVector = outgoingVector;
        transmissionRay.startingPosition = args.intersectionPoint;
        Hit transmissionHit = findClosestHit(transmissionRay, args.objects, args.numberOfObjects);
        vec3 attenuationColor;
        double distance = transmissionHit.distance;
        if (distance > 0 && args.outside){
            vec3 log_attenuation = absorptionAlbedo * attenuationCoefficient * (-distance);
            attenuationColor = expVector(log_attenuation);
        }
        else{
            attenuationColor = WHITE;
        }

        return attenuationColor;
    }

    brdfData sampleTransmission(const microfacetSampleArgs& args){
        vec3 refractedVector = refractVector(-args.incidentVector, -args.sampledHalfVector, args.eta);

        if (refractedVector.length() == 0){
            return sampleReflection(args);
        }

        vec3 attenuatedColor = computeAttenuatedColor(args, refractedVector);

        brdfData data;
        data.outgoingVector = refractedVector;

        data.brdfMultiplier = attenuatedColor * G(args.sampledHalfVector,args.normalVector, args.incidentVector, data.outgoingVector, args.alpha) * args.cosineFactor;
        data.type = TRANSMITTED;
        return data;
    }
    
    
    brdfData sample(const Hit& hit, Object** objects, const int numberOfObjects, const double u, const double v) override{
        microfacetData data = prepareMicrofacetData(hit, u, v);
                
        double iDotH = dotVectors(hit.incomingVector, data.halfVector);
        double iDotN = dotVectors(hit.incomingVector, data.normalIntoInterface);
        double nDotH = dotVectors(data.halfVector, data.normalIntoInterface);

        double cosineFactor = std::abs(iDotH / (iDotN * nDotH));

        double randomNum = randomUniform(0, 1);
        bool reflectSpecular = randomNum < data.F_r;
        
        microfacetSampleArgs args;
        args.sampledHalfVector = data.halfVector;
        args.normalVector = -data.normalIntoInterface;
        args.incidentVector = -hit.incomingVector;
        args.cosineFactor = cosineFactor;
        args.u = u;
        args.v = v;
        args.alpha = data.alpha;
        if (reflectSpecular){
            return sampleReflection(args);
        }
        else{
            args.intersectionPoint = hit.intersectionPoint;
            args.eta = data.eta;
            args.objects = objects;
            args.numberOfObjects = numberOfObjects;
            args.outside = data.outside;
            double randomNum2 = randomUniform(0, 1);
            bool diffuse = randomNum2 < percentageDiffuseMap -> get(u, v);
            return diffuse ? sampleDiffuse(args) : sampleTransmission(args);
        }
    }
};

#endif