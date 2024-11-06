#ifndef OBJECTS_H
#define OBJECTS_H
#include "vec3.h"
#include "utils.h"
#include "materials.h"
#include "constants.h"
#include "colors.h"


Hit findClosestHit(const Ray& ray, Object** objects, const int size);


class Object{
    public:
        Material* material;
        double area;
        int objectID = 0;
        Object(){}
        Object(Material* _material){
            material = _material;
        }

        virtual vec3 maxAxisPoint(){
            throw VirtualMethodNotAllowedException("maxAxisPoint is a pure virtual method and should not be called.");
            vec3 point;
            return point;
        }

        virtual vec3 minAxisPoint(){
            throw VirtualMethodNotAllowedException("minAxisPoint is a pure virtual method and should not be called.");
            vec3 point;
            return point;
        }

        virtual vec3 computeCentroid(){
            throw VirtualMethodNotAllowedException("computeCentroid is a pure virtual method and should not be called.");
            vec3 centroid;
            return centroid;
        }

        virtual vec3 getUV(const Hit& hit){
            throw VirtualMethodNotAllowedException("this is a pure virtual method and should not be called.");
            vec3 vec;
            return vec;
        }

        virtual bool isLightSource(){
            return material -> isLightSource;
        }

        virtual vec3 eval(const Hit& hit){
            vec3 UV = getUV(hit);
            return material -> eval(hit, UV[0], UV[1]);
        }

        virtual brdfData sample(const Hit& hit, Object** objectPtrList, const int numberOfObjects){
            vec3 UV = getUV(hit);
            return material -> sample(hit, objectPtrList, numberOfObjects, UV[0], UV[1]);
        }

        virtual vec3 getLightEmittance(const Hit& hit){
            vec3 UV = getUV(hit);
            return material -> getLightEmittance(UV[0], UV[1]);
        }

        virtual Hit findClosestObjectHit(const Ray& ray){
            throw VirtualMethodNotAllowedException("this is a pure virtual method and should not be called.");
            Hit hit;
            return hit;
        }
        
        virtual vec3 getNormalVector(const vec3& surfacePoint, const int objectID){
            throw VirtualMethodNotAllowedException("this is a pure virtual method and should not be called.");
            vec3 vec;
            return vec;
        }

        virtual vec3 generateRandomSurfacePoint(){
            throw VirtualMethodNotAllowedException("this is a pure virtual method and should not be called.");
            vec3 point;
            return point;
        }

        double areaToAnglePDFFactor(const vec3& surfacePoint, const vec3& intersectionPoint, const int objectID){
            vec3 normalVector = getNormalVector(surfacePoint, objectID);
            vec3 differenceVector = intersectionPoint - surfacePoint;
            vec3 vectorToPoint = normalize_vector(differenceVector);
            double inversePDF = dot_vectors(normalVector, vectorToPoint) / differenceVector.length_squared();
            return std::max(0.0, inversePDF);
        }

        virtual vec3 randomLightPoint(const vec3& intersectionPoint, double& inversePDF){
            vec3 randomPoint = generateRandomSurfacePoint();
            inversePDF = area * areaToAnglePDFFactor(randomPoint, intersectionPoint, 0);
            return randomPoint;
        }
};


class Sphere: public Object{
    public:
        vec3 position;
        double radius;
        double radiusSquared;
        Sphere(){}
        Sphere(const vec3& _position, const double _radius) : Object(){
            position = _position;
            radius = _radius;
            area = 4 * M_PI * radius * radius;
            radiusSquared = radius * radius;
        }
        Sphere(const vec3& _position, const double _radius, Material*_material) : Object(_material){
            position = _position;
            radius = _radius;
            area = 4 * M_PI * radius * radius;
            radiusSquared = radius * radius;
        }

        vec3 getUV(const Hit& hit) override{
            vec3 unitSpherePoint = (hit.intersection_point - position) / radius;
            double x = -unitSpherePoint[0];
            double y = -unitSpherePoint[1];
            double z = -unitSpherePoint[2];
            double u = 0.5 + atan2(z, x) / (2 * M_PI);
            double v = 0.5 + asin(y) / (M_PI);
            return vec3(u, v, 0);
        }

        Hit findClosestObjectHit(const Ray& ray) override{

            double dotProduct = dot_vectors(ray.direction_vector, ray.starting_position);
            double b = 2 * (dotProduct - dot_vectors(ray.direction_vector, position));
            vec3 difference_in_positions = position - ray.starting_position;
            double c = difference_in_positions.length_squared() - radiusSquared;
            double distance = solve_quadratic(b, c);
            Hit hit;
            hit.object_ID = objectID;
            hit.distance = distance;
            return hit;
        }

        vec3 getNormalVector(const vec3& surfacePoint, const int objectID) override{
            vec3 differenceVector = surfacePoint - position;
            return normalize_vector(differenceVector);
        }

        vec3 generateRandomSurfacePoint() override{
            return sample_spherical() * radius + position;
        }

        vec3 randomLightPoint(const vec3& intersectionPoint, double& inversePDF) override{
            double distance = (intersectionPoint - position).length();
            if (distance <= radius){
                vec3 randomPoint = generateRandomSurfacePoint();
                inversePDF = area * areaToAnglePDFFactor(randomPoint, intersectionPoint, 0);
                return randomPoint;
            }
        
            double cosThetaMax = sqrt(1 - pow(radius / distance, 2));
            inversePDF = 2 * M_PI * (1 - (cosThetaMax));

            double rand = random_uniform(0, 1);
            double cosTheta = 1 + rand * (cosThetaMax-1);
            double sinTheta = sqrt(1 - cosTheta * cosTheta);
            double cosAlpha = (radiusSquared + distance * distance - pow(distance * cosTheta - sqrt(radiusSquared - pow(distance*sinTheta, 2)), 2)) / (2.0 * distance * radius);
            double sinAlpha = sqrt(1.0 - cosAlpha * cosAlpha);
            
            vec3 xHat;
            vec3 yHat;
            vec3 zHat = getNormalVector(intersectionPoint, 0);
            set_perpendicular_vectors(zHat, xHat, yHat);
            double phi = random_uniform(0, 2.0*M_PI);
            vec3 randomPoint = xHat * sinAlpha * cos(phi) + yHat * sinAlpha * sin(phi) + zHat * cosAlpha;
            return randomPoint * radius + position;
        }
};


class Plane: public Object{
    public:
        vec3 position;
        vec3 v1;
        vec3 v2;
        vec3 normalVector;
        bool transparentBack;
        Plane(){}
        Plane(const vec3& _position, const vec3& _v1, const vec3& _v2, Material*_material) : Object(_material){
            position = _position;
            v1 = normalize_vector(_v1);
            v2 = normalize_vector(_v2);
            vec3 _normalVector = cross_vectors(v1, v2);
            normalVector = normalize_vector(_normalVector);
        }

        vec3 getUV(const Hit& hit) override{
            vec3 shiftedPoint = hit.intersection_point - position;
            double u = 1 - dot_vectors(shiftedPoint, v1) - 0.5;
            double v = 1 - dot_vectors(shiftedPoint, v2) - 0.5;
            return vec3(u, v, 0);
        }

        double computeDistanceInCenteredSystem(const vec3& startingPoint, const vec3& directionVector){
            double directionDotNormal = -dot_vectors(directionVector, normalVector);
            if (std::abs(directionDotNormal) < constants::EPSILON){
                return -1;
            }

            double distancesToStart = dot_vectors(startingPoint, normalVector);
            return distancesToStart / directionDotNormal;
        }

        Hit findClosestObjectHit(const Ray& ray) override{
            vec3 shiftedPoint = ray.starting_position - position;
            double distance = computeDistanceInCenteredSystem(shiftedPoint, ray.direction_vector);
            Hit hit;
            hit.object_ID = objectID;
            hit.distance = distance;
            return hit;
        }

        vec3 getNormalVector(const vec3& surfacePoint, const int objectID) override{
            return normalVector;
        }

};


class Rectangle: public Plane{
    public:
        double L1;
        double L2;
        Rectangle(){}
        Rectangle(const vec3& _position, const vec3& _v1, const vec3& _v2, const double _L1, const double _L2, Material*_material) : Plane(_position, _v1, _v2, _material){
            L1 = _L1;
            L2 = _L2;
            area = L1 * L2;
        }
        vec3 getUV(const Hit& hit) override{
            vec3 shiftedPoint = hit.intersection_point - position;
            double u = 1 - dot_vectors(shiftedPoint, v1) / L1 - 0.5;
            double v = 1 - dot_vectors(shiftedPoint, v2) / L2 - 0.5;
            return vec3(u, v, 0);
        }

        Hit findClosestObjectHit(const Ray& ray) override{
            Hit hit;
            hit.object_ID = objectID;
            hit.distance = -1;

            vec3 shiftedPoint = ray.starting_position - position;
            double distance = Plane::computeDistanceInCenteredSystem(shiftedPoint, ray.direction_vector);
            if (distance < 0){
                return hit;
            }
            double directionDotV1 = dot_vectors(ray.direction_vector, v1);
            double directionDotV2 = dot_vectors(ray.direction_vector, v2);
            double startDotV1 = dot_vectors(shiftedPoint, v1);
            double startDotV2 = dot_vectors(shiftedPoint, v2);

            if (std::abs(startDotV1 + directionDotV1 * distance) > L1 / 2.0 + constants::EPSILON || std::abs(startDotV2 + directionDotV2 * distance) > L2 / 2.0 + constants::EPSILON){
                return hit;
            }
            hit.distance = distance;
            return hit;
        }

        vec3 generateRandomSurfacePoint() override{
            double r1 = random_uniform(-L1/2, L1/2);
            double r2 = random_uniform(-L2/2, L2/2);
            return v1 * r1 + v2 * r2 + position;
        }
};


class Triangle: public Object{
    public:
        vec3 position;
        vec3 normalVector;
        vec3 p1;
        vec3 p2;
        vec3 p3;
        vec3 v1;
        vec3 v2;
        double x1;
        double y1;
        double x2;
        double y2;
        double x3;
        double y3;
        double detT;

        vec3 uv1;
        vec3 uv2;
        vec3 uv3;
        vec3 n1;
        vec3 n2;
        vec3 n3;

        bool smoothShaded = false;
        Triangle(){}
        Triangle(const vec3& _p1, const vec3& _p2, const vec3& _p3, Material*_material) : Object(_material){
            p1 = _p1;
            p2 = _p2;
            p3 = _p3;

            position = p1;

            v1 = p2 - p1;
            v2 = p3 - p1;
            normalVector = normalize_vector(cross_vectors(v1, v2));
            v1 = normalize_vector(v1);
            v2 = normalize_vector(cross_vectors(normalVector, v1));

            x1 = dot_vectors(p1, v1);
            y1 = dot_vectors(p1, v2);
            x2 = dot_vectors(p2, v1);
            y2 = dot_vectors(p2, v2);
            x3 = dot_vectors(p3, v1);
            y3 = dot_vectors(p3, v2);
            detT = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3);

            area = 0.5 * std::abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));

            uv1 = vec3(0, 0, 0);
            uv2 = vec3(0, 0, 0);
            uv3 = vec3(0, 0, 0);

            n1 = normalVector;
            n2 = normalVector;
            n3 = normalVector;
        }

        vec3 maxAxisPoint() override{
            vec3 point;
            for (int i = 0; i < 3; i++){
                point.e[i] = std::max(std::max(p1[i], p2[i]), p3[i]);
            }
            return point;
        }

        vec3 minAxisPoint() override{
            vec3 point;
            for (int i = 0; i < 3; i++){
                point.e[i] = std::min(std::min(p1[i], p2[i]), p3[i]);
            }
            return point;
        }


        vec3 computeCentroid() override{
            return (p1 + p2 + p3) / 3.0;
        }

        void setVertexUV(const vec3& _uv1, const vec3& _uv2, const vec3& _uv3){
            uv1 = _uv1;
            uv2 = _uv2;
            uv3 = _uv3;
        }

        void setVertexNormals(const vec3& _n1, const vec3& _n2, const vec3& _n3){
            n1 = _n1;
            n2 = _n2;
            n3 = _n3;
            smoothShaded = true;
        }

        vec3 getNormalVector(const vec3& surfacePoint, const int objectID) override{
            if (smoothShaded){
                return getNormalVectorSmoothed(surfacePoint, objectID);
            }
            return normalVector;
        }

        vec3 getNormalVectorSmoothed(const vec3& surfacePoint, const int objectID){
            vec3 barycentricVector = computeBarycentric(surfacePoint);
            return normalize_vector(n1 * barycentricVector[0] + n2 * barycentricVector[1] + n3 * barycentricVector[2]);
        }

        vec3 computeBarycentric(const vec3& point){
            double x = dot_vectors(point, v1);
            double y = dot_vectors(point, v2);

            double lambda1 = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / detT;
            double lambda2 = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / detT;
            return vec3(lambda1, lambda2, 1.0 - lambda1 - lambda2);
        }

        vec3 getUV(const Hit& hit) override{
            vec3 barycentricVector = computeBarycentric(hit.intersection_point);
            return uv1 * barycentricVector[0] + uv2 * barycentricVector[1] + uv3 * barycentricVector[2];
        }

        Hit findClosestObjectHit(const Ray& ray) override{
            Hit hit;
            hit.object_ID = objectID;
            hit.distance = -1;
            vec3 shiftedPoint = ray.starting_position - position;

            double directionDotNormal = -dot_vectors(ray.direction_vector, normalVector);
            if (std::abs(directionDotNormal) < constants::EPSILON){
                return hit;
            }

            double distancesToStart = dot_vectors(shiftedPoint, normalVector);
            double distance = distancesToStart / directionDotNormal;

            if (distance < constants::EPSILON){
                return hit;
            }

            vec3 inPlanePoint = ray.starting_position + ray.direction_vector * distance;

            vec3 barycentricVector = computeBarycentric(inPlanePoint);
            if (barycentricVector[0] < 0 || barycentricVector[1] < 0 || barycentricVector[2] < 0){
                return hit;
            }
            hit.distance = distance;
            return hit;
        }

        vec3 generateRandomSurfacePoint() override{
            double r1 = random_uniform(0, 1);
            double r2 = random_uniform(0, 1);
            return p1 * (1.0 - sqrt(r1)) + p2 * (sqrt(r1) * (1.0 - r2)) + p3 * (sqrt(r1) * r2);
        }
};


Hit findClosestHit(const Ray& ray, Object** objects, const int size){
    Hit closestHit;
    closestHit.distance = -1;
    for (int i = 0; i < size; i++){
        Hit hit = objects[i] -> findClosestObjectHit(ray);
        if (hit.distance > constants::EPSILON && (hit.distance < closestHit.distance || closestHit.distance == -1)){
            hit.intersected_object_index = i;
            closestHit = hit;
        }
    }
    if (closestHit.distance < constants::EPSILON){
        return closestHit;
    }

    closestHit.intersection_point = ray.starting_position + ray.direction_vector * closestHit.distance;
    closestHit.normal_vector = objects[closestHit.intersected_object_index] -> getNormalVector(closestHit.intersection_point, closestHit.object_ID);
    closestHit.incoming_vector = ray.direction_vector;
    return closestHit;
 }


#endif