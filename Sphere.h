#ifndef OBJECTS_H
#define OBJECTS_H
#include "vec3.h"
#include "utils.h"
#include "materials.h"


struct Hit{
    int intersectedObjectIndex;
    int objectID;
    float distance;
    vec3 intersectionPoint;
    vec3 incomingVector;
    vec3 outgoingVector;
    vec3 normalVector;
};


class Object{
    public:
        vec3 position;
        Material material;
        Object(){}
        Object(vec3 _position, Material _material){
            position = _position;
            material = _material;
        }
        ~Object(){}

        virtual Hit findClosestHit(Ray ray){
            Hit hit;
            return hit;
        }
        virtual vec3 getNormalVector(Hit hit){
            vec3 vec;
            return vec;
        }
};


class Sphere: public Object{
    public:
        float radius;
        Sphere(){}
        Sphere(vec3 _position, float _radius, Material _material){
            position = _position;
            radius = _radius;
            material = _material;
        }
        ~Sphere(){}
        Hit findClosestHit(Ray ray) override{

            float dotProduct = dotVectors(ray.directionVector, ray.startingPosition);
            float b = 2 * (dotProduct - dotVectors(ray.directionVector, position));
            vec3 difference_in_positions = subtractVectors(position, ray.startingPosition);
            float c = dotVectors(difference_in_positions, difference_in_positions) - pow(radius, 2);
            float distance = solveQuadratic(b, c);
            Hit hit;
            hit.objectID = 0;
            hit.distance = distance;
            return hit;
        }

        vec3 getNormalVector(Hit hit) override{
            return subtractVectors(position, hit.intersectionPoint);
        }
};


#endif