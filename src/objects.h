#ifndef OBJECTS_H
#define OBJECTS_H
#include "vec3.h"
#include "utils.h"
#include "materials.h"
#include "constants.h"
#include "colors.h"


class Object{
    public:
        vec3 position;
        Material* material;
        double area;
        Object(){}
        Object(vec3 _position, Material* _material){
            position = _position;
            material = _material;
        }
        virtual ~Object(){}

        virtual Hit findClosestHit(const Ray& ray){
            Hit hit;
            return hit;
        }
        virtual vec3 getNormalVector(const vec3& intersectionPoint){
            vec3 vec;
            return vec;
        }
        virtual vec3 generateRandomSurfacePoint(){
            vec3 point;
            return point;
        }

        virtual vec3 randomLightPoint(const vec3& referencePoint, double& inversePDF){
            vec3 point;
            return point;
        }
        double areaToAnglePDFFactor(const vec3& surfacePoint, const vec3& referencePoint){
            vec3 normalVector = getNormalVector(surfacePoint);
            vec3 differenceVector = referencePoint - surfacePoint;
            vec3 vectorToPoint = normalizeVector(differenceVector);
            double PDF = dotVectors(normalVector, vectorToPoint) / differenceVector.length_squared();
            return std::max(0.0, PDF);
        }
};


class Sphere: public Object{
    public:
        double radius;
        Sphere(){}
        Sphere(vec3 _position, double _radius, Material* _material){
            position = _position;
            radius = _radius;
            material = _material;
            area = 4 * M_PI * pow(radius, 2);
        }
        Hit findClosestHit(const Ray& ray) override{

            double dotProduct = dotVectors(ray.directionVector, ray.startingPosition);
            double b = 2 * (dotProduct - dotVectors(ray.directionVector, position));
            vec3 difference_in_positions = position - ray.startingPosition;
            double c = difference_in_positions.length_squared() - pow(radius, 2);
            double distance = solveQuadratic(b, c);
            Hit hit;
            hit.objectID = 0;
            hit.distance = distance;
            return hit;
        }

        vec3 getNormalVector(const vec3& intersectionPoint) override{
            vec3 differenceVector = intersectionPoint - position;
            return normalizeVector(differenceVector);
        }

        vec3 generateRandomSurfacePoint() override{
            return sampleSpherical() * radius + position;
        }

        vec3 randomLightPoint(const vec3& referencePoint, double& inversePDF) override{
            double distanceToPoint = (referencePoint - position).length();
            if (distanceToPoint <= radius){
                vec3 randomPoint = generateRandomSurfacePoint();
                inversePDF = areaToAnglePDFFactor(randomPoint, referencePoint) * area;
                return randomPoint;
            }
            
            double thetaMax = asin(radius / distanceToPoint);
            double cosMax = cos(M_PI / 2 - thetaMax);
            vec3 randomPoint = sampleAngledHemisphere(getNormalVector(referencePoint), cosMax) * radius + position;

            inversePDF = 2 * M_PI * pow(radius, 2) * (1 - cosMax);
            inversePDF *= areaToAnglePDFFactor(randomPoint, referencePoint);
            return randomPoint;
        }
};


class Plane: public Object{
    public:
        vec3 v1;
        vec3 v2;
        vec3 normalVector;
        bool transparentBack;
        Plane(){}
        Plane(vec3 _position, vec3 _v1, vec3 _v2, Material* _material){
            position = _position;
            v1 = normalizeVector(_v1);
            v2 = normalizeVector(_v2);
            vec3 _normalVector = crossVectors(v1, v2);
            normalVector = normalizeVector(_normalVector);
            material = _material;
        }

        double computeDistanceInCenteredSystem(const vec3& startingPoint, const vec3& directionVector){
            double directionDotNormal = -dotVectors(directionVector, normalVector);
            if (std::abs(directionDotNormal) < constants::EPSILON){
                return -1;
            }

            double distancesToStart = dotVectors(startingPoint, normalVector);
            double distances = distancesToStart / directionDotNormal;
            return distances;
        }

        Hit findClosestHit(const Ray& ray) override{
            vec3 shiftedPoint = ray.startingPosition - position;
            double distance = computeDistanceInCenteredSystem(shiftedPoint, ray.directionVector);
            Hit hit;
            hit.objectID = 0;
            hit.distance = distance;
            return hit;
        }

        vec3 getNormalVector(const vec3& intersectionPoint) override{
            return normalVector;
        }

};


class Rectangle: public Plane{
    public:
        double L1;
        double L2;
        Rectangle(){}
        Rectangle(vec3 _position, vec3 _v1, vec3 _v2, double _L1, double _L2, Material* _material){
            position = _position;
            v1 = normalizeVector(_v1);
            v2 = normalizeVector(_v2);
            L1 = _L1;
            L2 = _L2;
            area = L1 * L2;
            vec3 _normalVector = crossVectors(v1, v2);
            normalVector = normalizeVector(_normalVector);
            material = _material;
        }

        Hit findClosestHit(const Ray& ray) override{
            Hit hit;
            hit.objectID = 0;

            vec3 shiftedPoint = ray.startingPosition - position;
            double distance = Plane::computeDistanceInCenteredSystem(shiftedPoint, ray.directionVector);
            if (distance < 0){
                hit.distance = distance;
                return hit;
            }
            double directionDotV1 = dotVectors(ray.directionVector, v1);
            double directionDotV2 = dotVectors(ray.directionVector, v2);
            double startDotV1 = dotVectors(shiftedPoint, v1);
            double startDotV2 = dotVectors(shiftedPoint, v2);

            if (std::abs(startDotV1 + directionDotV1 * distance) > L1 / 2.0 + constants::EPSILON || std::abs(startDotV2 + directionDotV2 * distance) > L2 / 2.0 + constants::EPSILON){
                distance = -1;
            }
            hit.distance = distance;
            return hit;
        }

        vec3 generateRandomSurfacePoint() override{
            double r1 = randomUniform(-L1/2, L1/2);
            double r2 = randomUniform(-L2/2, L2/2);
            vec3 randomPoint = v1 * r1 + v2 * r2 + position;
            return randomPoint;
        }

        vec3 randomLightPoint(const vec3& referencePoint, double& inversePDF) override{
            vec3 randomPoint = generateRandomSurfacePoint();
            inversePDF = area * areaToAnglePDFFactor(randomPoint, referencePoint);
            return randomPoint;
        }
};

Hit findClosestHit(Ray& ray, Object* objects[], int size){
    Hit closestHit;
    closestHit.distance = -1;

    for (int i = 0; i < size; i++){
        Hit hit = objects[i] -> findClosestHit(ray);
        if (hit.distance > constants::EPSILON && (hit.distance < closestHit.distance || closestHit.distance == -1)){
            hit.intersectedObjectIndex = i;
            closestHit = hit;
        }
    }
    if (closestHit.distance < constants::EPSILON){
        return closestHit;
    }

    closestHit.intersectionPoint = ray.startingPosition + ray.directionVector * closestHit.distance;
    closestHit.normalVector = objects[closestHit.intersectedObjectIndex] -> getNormalVector(closestHit.intersectionPoint);
    closestHit.incomingVector = ray.directionVector;
    return closestHit;
 }


#endif