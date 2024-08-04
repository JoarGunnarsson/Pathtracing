#ifndef materials
#define materials
#include "vec3.h"
#include "colors.h"
#include "utils.h"
#include "constants.h"


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
        vec3 absorptionColor;
        vec3 emmissionColor;
        double lightIntensity;
        Material(){
            albedo = WHITE;
            refractiveIndex = 1;
            attenuationCoefficient = 0;
            absorptionColor = WHITE;
            attenuationCoefficient = 0;
            absorptionColor = WHITE;
            emmissionColor = WHITE;
            lightIntensity = 0;
            // TODO: Combine emmissionColor and lightIntensity into one variable.
        }
        Material(vec3 _diffuseColor, double _diffuseCoefficient=0.8, double _refractiveIndex=1, double _attenuationCoefficient=0, vec3 _absorptionColor=WHITE,
        vec3 _emmissionColor=WHITE, double _lightIntensity=0){
            albedo = _diffuseColor * _diffuseCoefficient;
            refractiveIndex = _refractiveIndex;
            attenuationCoefficient = _attenuationCoefficient;
            absorptionColor = _absorptionColor;
            emmissionColor = _emmissionColor;
            lightIntensity = _lightIntensity;
        }

    virtual vec3 eval(){
        throw VirtualMethodNotAllowedException("this is a pure virtual method and should not be called.");
        vec3 vec;
        return vec;
    }

    virtual brdfData sample(Hit& hit){
        throw VirtualMethodNotAllowedException("this is a pure virtual method and should not be called.");
        brdfData data;
        return data;
    }
};


class DiffuseMaterial : public Material{
    public:
        using Material::Material;

    vec3 eval() override{
        return albedo / M_PI;
    }

    brdfData sample(Hit& hit) override{
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

    brdfData sample(Hit& hit) override{
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

    brdfData sample(Hit& hit) override{
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

        /*
        vec3 randomHemispherePoint = sampleHemisphere(fresnelNormal);
        double smoothness = 0.5;
        vec3 scaledTransmittedVector = multiplyVector(transmittedVector, smoothness);
        vec3 scaledHemispherePoint = multiplyVector(randomHemispherePoint, 1-smoothness);
        transmittedVector = addVectors(scaledTransmittedVector, scaledHemispherePoint);
        transmittedVector = normalizeVector(transmittedVector);
        */
        double F_r = 1;
        if (transmittedVector.length_squared() != 0){
            F_r = fresnelMultiplier(hit.incomingVector, transmittedVector, fresnelNormal, n1, n2);
        }
        
        double randomNum = randomUniform(0, 1);
        bool isReflected = randomNum <= F_r;

        /*
        # TODO: Add attenuation computation here.
        transmission_from_outside_indices = np.logical_and(transmitted_indices, outside_indices)
        attenuation_intersections = intersection_points[transmission_from_outside_indices]
        attenuation_outgoing_vectors = transmitted_vectors[transmission_from_outside_indices]
        _, transmitted_distances = utils.find_closest_intersected_object(attenuation_intersections, attenuation_outgoing_vectors, Scene.current_scene.objects)

        attenuation_factors = np.ones((transmission_size, 3), dtype=float)
        transmitted_distances = np.maximum(0, transmitted_distances)
        combined_attenuation_coefficient = self.attenuation_coefficient * self.absorption_color
        attenuation_factors[transmitted_into_an_object_indices] = np.exp(-combined_attenuation_coefficient * transmitted_distances)
        */
        vec3 reflectedVector = reflectVector(hit.incomingVector, fresnelNormal);
        // TODO: Ok vectors^?
        vec3 brdfMultiplier;
        vec3 outgoingVector;
        if (isReflected){
            brdfMultiplier = albedo;
            outgoingVector = reflectedVector;
        }
        else{
            Ray transmissionRay;
            transmissionRay.directionVector = transmittedVector;
            transmissionRay.startingPosition = hit.intersectionPoint;
            //Hit transmissionHit = findClosestHit(transmissionRay, objectPtrList);
            double attenuationFactor = 1;
            vec3 attenuationColor;
            double distance = 0;
            if (distance > 0 and !inside){
                vec3 combined_attenuation_coefficient = absorptionColor * attenuationCoefficient;
                vec3 log_attenuation = combined_attenuation_coefficient * -distance;
                attenuationColor = expVector(log_attenuation);
            }
            else{
                attenuationColor = albedo;
            }

            double refractionIntensityFactor = pow(n2, 2) / pow(n1, 2);
            brdfMultiplier = attenuationColor * refractionIntensityFactor; // White -> attenuationFactors.
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