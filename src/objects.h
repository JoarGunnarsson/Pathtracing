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
        Object(){}
        Object(Material* _material){
            material = _material;
        }

        virtual vec3 getUV(const vec3& hitPoint){
            throw VirtualMethodNotAllowedException("this is a pure virtual method and should not be called.");
            vec3 vec;
            return vec;
        }

        virtual bool isLightSource(const Hit& hit){
            return material -> isLightSource;
        }

        virtual vec3 eval(const Hit& hit){
            vec3 UV = getUV(hit.intersectionPoint);
            return material -> eval(hit, UV[0], UV[1]);
        }

        virtual brdfData sample(const Hit& hit, Object** objectPtrList, const int numberOfObjects){
            vec3 UV = getUV(hit.intersectionPoint);
            return material -> sample(hit, objectPtrList, numberOfObjects, UV[0], UV[1]);
        }

        virtual vec3 getLightEmittance(const Hit& hit){
            vec3 UV = getUV(hit.intersectionPoint);
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

        double areaToAnglePDFFactor(const vec3& surfacePoint, const Hit& hit){
            vec3 normalVector = getNormalVector(surfacePoint, hit.objectID);
            vec3 differenceVector = hit.intersectionPoint - surfacePoint;
            vec3 vectorToPoint = normalizeVector(differenceVector);
            double inversePDF = dotVectors(normalVector, vectorToPoint) / differenceVector.length_squared();
            return std::max(0.0, inversePDF);
        }

        virtual vec3 randomLightPoint(const Hit& hit, double& inversePDF){
            vec3 randomPoint = generateRandomSurfacePoint();
            inversePDF = area * areaToAnglePDFFactor(randomPoint, hit);
            return randomPoint;
        }
};


class Sphere: public Object{
    public:
        vec3 position;
        double radius;
        double radiusSquared;
        Sphere(){}
        Sphere(vec3 _position, double _radius, Material* _material) : Object(_material){
            position = _position;
            radius = _radius;
            area = 4 * M_PI * radius * radius;
            radiusSquared = radius * radius;
        }

        vec3 getUV(const vec3& hitPoint) override{
            vec3 unitSpherePoint = (hitPoint - position) / radius;
            double x = -unitSpherePoint[0];
            double y = -unitSpherePoint[1];
            double z = -unitSpherePoint[2];
            double u = 0.5 + atan2(z, x) / (2 * M_PI);
            double v = 0.5 + asin(y) / (M_PI);
            return vec3(u, v, 0);
        }

        Hit findClosestObjectHit(const Ray& ray) override{

            double dotProduct = dotVectors(ray.directionVector, ray.startingPosition);
            double b = 2 * (dotProduct - dotVectors(ray.directionVector, position));
            vec3 difference_in_positions = position - ray.startingPosition;
            double c = difference_in_positions.length_squared() - radiusSquared;
            double distance = solveQuadratic(b, c);
            Hit hit;
            hit.objectID = 0;
            hit.distance = distance;
            return hit;
        }

        vec3 getNormalVector(const vec3& surfacePoint, const int objectID) override{
            vec3 differenceVector = surfacePoint - position;
            return normalizeVector(differenceVector);
        }

        vec3 generateRandomSurfacePoint() override{
            return sampleSpherical() * radius + position;
        }

        vec3 randomLightPoint(const Hit& hit, double& inversePDF) override{
            double distance = (hit.intersectionPoint - position).length();
            if (distance <= radius){
                vec3 randomPoint = generateRandomSurfacePoint();
                inversePDF = areaToAnglePDFFactor(randomPoint, hit) * area;
                return randomPoint;
            }
        
            double cosThetaMax = sqrt(1 - pow(radius / distance, 2));
            inversePDF = 2 * M_PI * (1 - (cosThetaMax));

            double rand = randomUniform(0, 1);
            double cosTheta = 1 + rand * (cosThetaMax-1);
            double sinTheta = sqrt(1 - cosTheta * cosTheta);
            double cosAlpha = (radiusSquared + distance * distance - pow(distance * cosTheta - sqrt(radiusSquared - pow(distance*sinTheta, 2)), 2)) / (2.0 * distance * radius);
            double sinAlpha = sqrt(1.0 - cosAlpha * cosAlpha);
            
            vec3 xHat;
            vec3 yHat;
            vec3 zHat = getNormalVector(hit.intersectionPoint, hit.objectID);
            setPerpendicularVectors(zHat, xHat, yHat);
            double phi = randomUniform(0, 2.0*M_PI);
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
        Plane(vec3 _position, vec3 _v1, vec3 _v2, Material* _material) : Object(_material){
            position = _position;
            v1 = normalizeVector(_v1);
            v2 = normalizeVector(_v2);
            vec3 _normalVector = crossVectors(v1, v2);
            normalVector = normalizeVector(_normalVector);
        }

        vec3 getUV(const vec3& hitPoint) override{
            vec3 shiftedPoint = hitPoint - position;
            double u = 1 - dotVectors(shiftedPoint, v1) - 0.5;
            double v = 1 - dotVectors(shiftedPoint, v2) - 0.5;
            return vec3(u, v, 0);
        }

        double computeDistanceInCenteredSystem(const vec3& startingPoint, const vec3& directionVector){
            double directionDotNormal = -dotVectors(directionVector, normalVector);
            if (std::abs(directionDotNormal) < constants::EPSILON){
                return -1;
            }

            double distancesToStart = dotVectors(startingPoint, normalVector);
            return distancesToStart / directionDotNormal;
        }

        Hit findClosestObjectHit(const Ray& ray) override{
            vec3 shiftedPoint = ray.startingPosition - position;
            double distance = computeDistanceInCenteredSystem(shiftedPoint, ray.directionVector);
            Hit hit;
            hit.objectID = 0;
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
        Rectangle(vec3 _position, vec3 _v1, vec3 _v2, double _L1, double _L2, Material* _material) : Plane(_position, _v1, _v2, _material){
            L1 = _L1;
            L2 = _L2;
            area = L1 * L2;
        }
        vec3 getUV(const vec3& hitPoint) override{
            vec3 shiftedPoint = hitPoint - position;
            double u = 1 - dotVectors(shiftedPoint, v1) / L1 - 0.5;
            double v = 1 - dotVectors(shiftedPoint, v2) / L2 - 0.5;
            return vec3(u, v, 0);
        }

        Hit findClosestObjectHit(const Ray& ray) override{
            Hit hit;
            hit.objectID = 0;
            hit.distance = -1;

            vec3 shiftedPoint = ray.startingPosition - position;
            double distance = Plane::computeDistanceInCenteredSystem(shiftedPoint, ray.directionVector);
            if (distance < 0){
                return hit;
            }
            double directionDotV1 = dotVectors(ray.directionVector, v1);
            double directionDotV2 = dotVectors(ray.directionVector, v2);
            double startDotV1 = dotVectors(shiftedPoint, v1);
            double startDotV2 = dotVectors(shiftedPoint, v2);

            if (std::abs(startDotV1 + directionDotV1 * distance) > L1 / 2.0 + constants::EPSILON || std::abs(startDotV2 + directionDotV2 * distance) > L2 / 2.0 + constants::EPSILON){
                return hit;
            }
            hit.distance = distance;
            return hit;
        }

        vec3 generateRandomSurfacePoint() override{
            double r1 = randomUniform(-L1/2, L1/2);
            double r2 = randomUniform(-L2/2, L2/2);
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

        Triangle(vec3 _p1, vec3 _p2, vec3 _p3, Material* _material) : Object(_material){
            p1 = _p1;
            p2 = _p2;
            p3 = _p3;

            position = p1;

            v1 = normalizeVector(p2 - p1);
            v2 = normalizeVector(p3 - p1);
            normalVector = normalizeVector(crossVectors(v1, v2));
            v2 = normalizeVector(crossVectors(normalVector, v1));

            x1 = dotVectors(p1, v1);
            y1 = dotVectors(p1, v2);
            x2 = dotVectors(p2, v1);
            y2 = dotVectors(p2, v2);
            x3 = dotVectors(p3, v1);
            y3 = dotVectors(p3, v2);
            detT = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3);

            area = 0.5 * std::abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));

            uv1 = vec3(0, 0, 0);
            uv2 = vec3(1, 0, 0);
            uv3 = vec3(0, 1, 0);

            n1 = normalVector;
            n2 = normalVector;
            n3 = normalVector;
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
            return normalizeVector(n1 * barycentricVector[0] + n2 * barycentricVector[1] + n3 * barycentricVector[2]);
        }

        vec3 computeBarycentric(const vec3& point){
            double x = dotVectors(point, v1);
            double y = dotVectors(point, v2);

            double lambda1 = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / detT;
            double lambda2 = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / detT;
            return vec3(lambda1, lambda2, 1.0 - lambda1 - lambda2);
        }

        vec3 getUV(const vec3& hitPoint) override{
            vec3 barycentricVector = computeBarycentric(hitPoint);
            return uv1 * barycentricVector[0] + uv2 * barycentricVector[1] + uv3 * barycentricVector[2];
        }

        Hit findClosestObjectHit(const Ray& ray) override{
            Hit hit;
            hit.objectID = 0;
            hit.distance = -1;
            vec3 shiftedPoint = ray.startingPosition - position;

            double directionDotNormal = -dotVectors(ray.directionVector, normalVector);
            if (std::abs(directionDotNormal) < constants::EPSILON){
                return hit;
            }

            double distancesToStart = dotVectors(shiftedPoint, normalVector);
            double distance = distancesToStart / directionDotNormal;

            if (distance < constants::EPSILON){
                return hit;
            }

            vec3 inPlanePoint = ray.startingPosition + ray.directionVector * distance;

            vec3 barycentricVector = computeBarycentric(inPlanePoint);
            if (barycentricVector[0] < 0 || barycentricVector[1] < 0 || barycentricVector[2] < 0){
                return hit;
            }
            hit.distance = distance;
            return hit;
        }

        vec3 generateRandomSurfacePoint() override{
            double r1 = randomUniform(0, 1);
            double r2 = randomUniform(0, 1);
            return p1 * (1.0 - sqrt(r1)) + p2 * (sqrt(r1) * (1.0 - r2)) + p3 * (sqrt(r1) * r2);
        }
};


class ObjectUnion : public Object{
    public:
        Object** objects;
        int numberOfObjects;
        vec3 position;
        double size;
        Sphere boundingSphere;
        ObjectUnion(Object** _objects, int _numberOfObjects, vec3 _position, double _size) : Object(){
            objects = _objects;
            numberOfObjects = _numberOfObjects;
            position = _position;
            size = _size;
            area = 0;
            for (int i = 0; i < numberOfObjects; i++){
                area += objects[i] -> area;
            }
            Material material;
            boundingSphere = Sphere(position, _size, &material);
        }

        ~ObjectUnion(){
            for (int i = 0; i < numberOfObjects; i++){
                delete objects[i];
            }
        }

        bool isLightSource(const Hit& hit) override{
            return objects[hit.objectID] -> isLightSource(hit);
        }

        vec3 eval(const Hit& hit) override{
            return objects[hit.objectID]  -> eval(hit);
        }

        brdfData sample(const Hit& hit, Object** objectPtrList, const int numberOfObjects) override{
            return objects[hit.objectID]  -> sample(hit, objectPtrList, numberOfObjects);
        }

        vec3 getLightEmittance(const Hit& hit) override{
            return objects[hit.objectID] -> getLightEmittance(hit);
        }

        Hit findClosestObjectHit(const Ray& ray) override{
            Hit boundingHit = boundingSphere.findClosestObjectHit(ray);
            if (boundingHit.distance < constants::EPSILON){
                return boundingHit;
            }
            Hit hit = findClosestHit(ray, objects, numberOfObjects);
            hit.objectID = hit.intersectedObjectIndex;
            hit.intersectedObjectIndex = -1;
            return hit;
        }
        
        vec3 getNormalVector(const vec3& surfacePoint, const int objectID) override{
            return objects[objectID] -> getNormalVector(surfacePoint, objectID);
        }

        int sampleRandomObjectIndex(){
            double randomAreaSplit = randomUniform(0, 1) * area;
            double summedArea = 0;
            int sampledIndex = 0;
            for (int i = 0; i < numberOfObjects; i++){
                summedArea += objects[i] -> area;
                if (summedArea >= randomAreaSplit){
                    sampledIndex = i;
                    break;
                }
            }
            return sampledIndex;
        }

        vec3 generateRandomSurfacePoint() override{
            return objects[sampleRandomObjectIndex()] -> generateRandomSurfacePoint();
        }

        vec3 randomLightPoint(const Hit& hit, double& inversePDF) override{
            vec3 randomPoint = objects[sampleRandomObjectIndex()] -> generateRandomSurfacePoint();
            inversePDF = area * areaToAnglePDFFactor(randomPoint, hit);
            return randomPoint;
        }
};


Hit findClosestHit(const Ray& ray, Object** objects, const int size){
    Hit closestHit;
    closestHit.distance = -1;
    for (int i = 0; i < size; i++){
        Hit hit = objects[i] -> findClosestObjectHit(ray);
        if (hit.distance > constants::EPSILON && (hit.distance < closestHit.distance || closestHit.distance == -1)){
            hit.intersectedObjectIndex = i;
            closestHit = hit;
        }
    }
    if (closestHit.distance < constants::EPSILON){
        return closestHit;
    }

    closestHit.intersectionPoint = ray.startingPosition + ray.directionVector * closestHit.distance;
    closestHit.normalVector = objects[closestHit.intersectedObjectIndex] -> getNormalVector(closestHit.intersectionPoint, closestHit.objectID);
    closestHit.incomingVector = ray.directionVector;
    return closestHit;
 }



struct dataSizes{
    int numVertices = 0;
    int numVertexUVs = 0;
    int numVertexNormals = 0;
    int numTriangles = 0;
};


std::string getNthWord(std::string line, std::string delimiter, int n){
    int endLocation = 0;
    int startLocation = 0;
    for (int i = 0; i < n + 1; i++){
        int newEndLocation = line.find(delimiter, endLocation)+1;
        if (newEndLocation <= endLocation && !(newEndLocation == 0 && endLocation == 0)){
            return line.substr(endLocation, line.size() - endLocation);
        }
        startLocation = endLocation;
        endLocation = newEndLocation;
    }
    return line.substr(startLocation, endLocation - startLocation - 1);
}


dataSizes getVertexDataSizes(std::string fileName){
    std::ifstream modelFile(fileName);
    std::string line;
    dataSizes nums;
    
    while(std::getline(modelFile, line)){
        std::string firstWord = getNthWord(line, " ", 0);

        bool isVertex = firstWord == "v";
        bool isVertexUV = firstWord == "vt";
        bool isVertexNormal = firstWord == "vn";
        bool isShape = firstWord == "f";
        
        if (isVertex){
            nums.numVertices++;
        }
        else if (isVertexUV){
            nums.numVertexUVs++;
        }
        else if (isVertexNormal){
            nums.numVertexNormals++;
        }
        else if (isShape){
            int numberOfSpaces = 0;
            for (int i = 0; i < line.size(); i++){
                if (line.substr(i, 1) == " "){
                    numberOfSpaces++;
                }
            }
            
            bool isTriangle = numberOfSpaces == 3;
            bool isQuad = numberOfSpaces == 4;
            if (isTriangle){
                nums.numTriangles++;
            }

            else if (isQuad){
                nums.numTriangles += 2;
            }
        }
    }
    return nums;
}


void populateVertexArrays(std::string fileName, vec3* vertexArray, vec3* vertexUVArray, vec3* vertexNormalArray){
    int vertexIdx = 0;
    int vertexUVIdx = 0;
    int vertexNormalIdx = 0;

    std::ifstream modelFile(fileName);
    std::string line;
    while(std::getline(modelFile, line)){
        std::string firstWord = getNthWord(line, " ", 0);
        bool isVertex = firstWord == "v";
        bool isVertexUV = firstWord == "vt";
        bool isVertexNormal = firstWord == "vn";

        if (isVertex){
            double v1 = std::stod(getNthWord(line, " ", 1));
            double v2 = std::stod(getNthWord(line, " ", 2));
            double v3 = std::stod(getNthWord(line, " ", 3));
            vertexArray[vertexIdx] = vec3(v1, v2, v3);
            vertexIdx++;
        }
        else if (isVertexUV){
            double u = std::stod(getNthWord(line, " ", 1));
            double v = std::stod(getNthWord(line, " ", 2));
            vertexUVArray[vertexUVIdx] = vec3(u, v, 0);
            vertexUVIdx++;
        }
        else if (isVertexNormal){
            double n1 = std::stod(getNthWord(line, " ", 1));
            double n2 = std::stod(getNthWord(line, " ", 2));
            double n3 = std::stod(getNthWord(line, " ", 3));
            vertexNormalArray[vertexNormalIdx] = vec3(n1, n2, n3);
            vertexNormalIdx++;
        }
    }
}


void populateVertexVectors(const std::string vertexData, vec3& v, vec3& uv, vec3& n, const vec3* vertexArray, const vec3* vertexUVArray, const vec3* vertexNormalArray){
    std::string vIdx = getNthWord(vertexData, "/", 0);
    std::string UVIdx = getNthWord(vertexData, "/", 1);
    std::string nIdx = getNthWord(vertexData, "/", 2);

    v = vertexArray[std::stoi(vIdx)-1];
    uv = vertexUVArray[std::stoi(UVIdx)-1];
    n = vertexNormalArray[std::stoi(nIdx)-1];
}


Triangle* constructTriangle(std::string triangleData, const int idx1, const int idx2, const int idx3, Material* material, const vec3* vertexArray, const vec3* vertexUVArray, const vec3* vertexNormalArray, const bool allowSmoothShading){
    std::string v1Data = getNthWord(triangleData, " ", idx1);

    vec3 v1;
    vec3 uv1;
    vec3 n1;
    populateVertexVectors(v1Data, v1, uv1, n1, vertexArray, vertexUVArray, vertexNormalArray);

    std::string v2Data = getNthWord(triangleData, " ", idx2);
    vec3 v2;
    vec3 uv2;
    vec3 n2;
    populateVertexVectors(v2Data, v2, uv2, n2, vertexArray, vertexUVArray, vertexNormalArray);

    std::string v3Data = getNthWord(triangleData, " ", idx3);
    vec3 v3;
    vec3 uv3;
    vec3 n3;
    populateVertexVectors(v3Data, v3, uv3, n3, vertexArray, vertexUVArray, vertexNormalArray);

    Triangle* triangle = new Triangle(v1, v2, v3, material);
    triangle -> setVertexUV(uv1, uv2, uv3);
    if (allowSmoothShading){
        triangle -> setVertexNormals(n1, n2, n3);
    }
    return triangle;
}


void populateTriangleArray(std::string fileName, vec3* vertexArray, vec3* vertexUVArray, vec3* vertexNormalArray, Object** triangleArray, Material* material, const bool allowSmoothShading){
    std::ifstream modelFile(fileName);
    std::string line;
    int shapeIdx = 0;
    while(std::getline(modelFile, line)){
        std::string firstWord = getNthWord(line, " ", 0);
        bool isShape = firstWord == "f";
        if (!isShape){
            continue;
        }
        int numberOfSpaces = 0;
        for (int i = 0; i < line.size(); i++){
            if (line.substr(i, 1) == " "){
                numberOfSpaces++;
            }
        }
        bool isTriangle = numberOfSpaces == 3;
        bool isQuad = numberOfSpaces == 4;
        if (isTriangle){
            
            Triangle* triangle = constructTriangle(line, 1, 2, 3, material, vertexArray, vertexUVArray, vertexNormalArray, allowSmoothShading);
            triangleArray[shapeIdx] = triangle;
            shapeIdx++;
        }

        else if (isQuad){
            Triangle* triangle1 = constructTriangle(line, 1, 2, 3, material, vertexArray, vertexUVArray, vertexNormalArray, allowSmoothShading);
            triangleArray[shapeIdx] = triangle1;

            Triangle* triangle2 = constructTriangle(line, 1, 3, 4, material, vertexArray, vertexUVArray, vertexNormalArray, allowSmoothShading);
            triangleArray[shapeIdx+1] = triangle2;
            shapeIdx += 2;
        }
    }
}


vec3 computeAveragePosition(vec3* vertexArray, const int numberOfVertices){
    vec3 avg = vec3(0,0,0);
    for (int i = 0; i < numberOfVertices; i++){
        avg += vertexArray[i];
    }
    return avg / numberOfVertices;
}


double maximumDistance(vec3 center, vec3* vertexArray, const int numberOfVertices){
    double maxDistance = 0;
    for (int i = 0; i < numberOfVertices; i++){
        double distance = (vertexArray[i] - center).length();
        if (distance > maxDistance){
            maxDistance = distance;
        }
    }
    return maxDistance;
}

void changeVectors(vec3 desiredCenter, double desiredSize, vec3* vertexArray, const int numberOfVertices){
    vec3 averagePosition = computeAveragePosition(vertexArray, numberOfVertices);
    double maxDistance = maximumDistance(averagePosition, vertexArray,numberOfVertices);

    for (int i = 0; i < numberOfVertices; i++){
        vertexArray[i] = (vertexArray[i] - averagePosition + desiredCenter) / maxDistance * desiredSize;
    }
}


ObjectUnion* loadObjectModel(std::string fileName, Material* material, const bool allowSmoothShading){
    dataSizes nums = getVertexDataSizes(fileName);

    vec3 vertexArray[nums.numVertices];
    vec3 vertexUVArray[nums.numVertexUVs];
    vec3 vertexNormalArray[nums.numVertexNormals];
    populateVertexArrays(fileName, vertexArray, vertexUVArray, vertexNormalArray);

    double desiredSize = 0.7;
    vec3 desiredCenter = vec3(0, 0.2, 1);
    changeVectors(desiredCenter, desiredSize, vertexArray, nums.numVertices);

    Object** triangles  = new Object*[nums.numTriangles];
    populateTriangleArray(fileName, vertexArray, vertexUVArray, vertexNormalArray, triangles, material, allowSmoothShading);
    ObjectUnion* loadedObject = new ObjectUnion(triangles, nums.numTriangles, desiredCenter, desiredSize);
    return loadedObject;
}

#endif