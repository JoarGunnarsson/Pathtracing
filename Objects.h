#ifndef OBJECTS_H
#define OBJECTS_H
#include "vec3.h"
#include "utils.h"
#include "materials.h"
#include "constants.h"


std::random_device this_rand_dev;
std::minstd_rand uniform_generator(this_rand_dev());
std::uniform_real_distribution<double> uniform_dist(0, 1);


class Object{
    public:
        vec3 position;
        Material material;
        double area;
        Object(){}
        Object(vec3 _position, Material _material){
            position = _position;
            material = _material;
        }
        virtual ~Object(){}

        virtual Hit findClosestHit(Ray& ray){
            Hit hit;
            return hit;
        }
        virtual vec3 getNormalVector(Hit& hit){
            vec3 vec;
            return vec;
        }
        virtual vec3 generateRandomSurfacePoint(){
            vec3 point;
            return point;
        }
};


class Sphere: public Object{
    public:
        double radius;
        Sphere(){}
        Sphere(vec3 _position, double _radius, Material _material){
            position = _position;
            radius = _radius;
            material = _material;
            area = 4 * M_PI * pow(radius, 2);
        }
        Hit findClosestHit(Ray& ray) override{

            double dotProduct = dotVectors(ray.directionVector, ray.startingPosition);
            double b = 2 * (dotProduct - dotVectors(ray.directionVector, position));
            vec3 difference_in_positions = subtractVectors(position, ray.startingPosition);
            double c = difference_in_positions.length_squared() - pow(radius, 2);
            double distance = solveQuadratic(b, c);
            Hit hit;
            hit.objectID = 0;
            hit.distance = distance;
            return hit;
        }

        vec3 getNormalVector(Hit& hit) override{
            vec3 differenceVector = subtractVectors(hit.intersectionPoint, position);
            return normalizeVector(differenceVector);
        }

        vec3 generateRandomSurfacePoint() override{
            vec3 randomPoint = sampleSpherical();
            randomPoint = multiplyVector(randomPoint, radius);
            vec3 point = addVectors(randomPoint, position);
            return point;
        }
};


class Plane: public Object{
    public:
        vec3 v1;
        vec3 v2;
        vec3 normalVector;
        Plane(){}
        Plane(vec3 _position, vec3 _v1=vec3(1,0,0), vec3 _v2=vec3(0,1,0), Material _material=Material()){
            position = _position;
            v1 = normalizeVector(_v1);
            v2 = normalizeVector(_v2);
            vec3 _normalVector = crossVectors(v1, v2);
            normalVector = normalizeVector(_normalVector);
            material = _material;
        }

        double computeDistanceInCenteredSystem(vec3& startingPoint, vec3& directionVector){
            double directionDotNormal = -dotVectors(directionVector, normalVector);
            if (directionDotNormal < constants::EPSILON){
                return -1;
            }

            double distancesToStart = dotVectors(startingPoint, normalVector);
            double distances = distancesToStart / directionDotNormal;
            return distances;
        }

        Hit findClosestHit(Ray& ray) override{
            vec3 shiftedPoint = subtractVectors(ray.startingPosition, position);
            double distance = computeDistanceInCenteredSystem(shiftedPoint, ray.directionVector);
            Hit hit;
            hit.objectID = 0;
            hit.distance = distance;
            return hit;
        }

        vec3 getNormalVector(Hit& hit) override{
            return normalVector;
        }

};


class Rectangle: public Plane{
    public:
        double L1;
        double L2;
        Rectangle(){}
        Rectangle(vec3 _position, vec3 _v1=vec3(1,0,0), vec3 _v2=vec3(0,1,0), double _L1=1, double _L2=1, Material _material=Material()){
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

        double computeDistanceInCenteredSystem(vec3& startingPoint, vec3& directionVector){
            double directionDotNormal = -dotVectors(directionVector, normalVector);
            if (directionDotNormal < constants::EPSILON){
                return -1;
            }

            double distanceToStart = dotVectors(startingPoint, normalVector);
            double distance = distanceToStart / directionDotNormal;
            return distance;
        }

        Hit findClosestHit(Ray& ray) override{
            vec3 shiftedPoint = subtractVectors(ray.startingPosition, position);
            double distance = Plane::computeDistanceInCenteredSystem(shiftedPoint, ray.directionVector);
            double directionDotV1 = dotVectors(ray.directionVector, v1);
            double directionDotV2 = dotVectors(ray.directionVector, v2);
            double startDotV1 = dotVectors(shiftedPoint, v1);
            double startDotV2 = dotVectors(shiftedPoint, v2);

            if (std::abs(startDotV1 + directionDotV1 * distance) > L1 / 2.0 + constants::EPSILON || std::abs(startDotV2 + directionDotV2 * distance) > L2 / 2.0 + constants::EPSILON){
                distance = -1;
            }
            Hit hit;
            hit.objectID = 0;
            hit.distance = distance;
            return hit;
        }

        vec3 generateRandomSurfacePoint() override{
            double r1 = (uniform_dist(uniform_generator) - 0.5) * L1;
            double r2 = (uniform_dist(uniform_generator) - 0.5) * L2;
            vec3 randomInV1 = multiplyVector(v1, r1);
            vec3 randomInV2 = multiplyVector(v2, r2);
            vec3 randomInLocal = addVectors(randomInV1, randomInV2);
            vec3 randomPoint = addVectors(randomInLocal, position);
            return randomPoint;
        }
};



#endif