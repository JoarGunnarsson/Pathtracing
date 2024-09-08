#ifndef OBJECTUNION_H
#define OBJECTUNION_H

#include <fstream>
#include "constants.h"
#include "vec3.h"
#include "utils.h"
#include "objects.h"


vec3 getMaxPoint(Object** triangles, int numberOfTriangles){
    vec3 maxPoint = triangles[0] -> maxAxisPoint();
    for (int i = 1; i < numberOfTriangles; i++){
        vec3 thisMaxPoint = triangles[i] -> maxAxisPoint();
        for (int j = 0; j < 3; j++){
            maxPoint.e[j] = std::max(thisMaxPoint[j], maxPoint[j]);
        }
    }
    return maxPoint;
}


vec3 getMinPoint(Object** triangles, int numberOfTriangles){
    vec3 minPoint = triangles[0] -> minAxisPoint();
    for (int i = 1; i < numberOfTriangles; i++){
        vec3 thisMinPoint = triangles[i] -> minAxisPoint();
        for (int j = 0; j < 3; j++){
            minPoint.e[j] = std::min(thisMinPoint[j], minPoint[j]);
        }
    }
    return minPoint;
}


class BoundingBox{
    public:
        Object** rects;
        int numberOfObjects = 6;
        BoundingBox(){}
        BoundingBox(Object** _triangles, int numberOfTriangles){
            vec3 maxPoint = getMaxPoint(_triangles, numberOfTriangles);
            vec3 minPoint = getMinPoint(_triangles, numberOfTriangles);
            
            double minX = minPoint[0];
            double minY = minPoint[1];
            double minZ = minPoint[2];

            double maxX = maxPoint[0];
            double maxY = maxPoint[1];
            double maxZ = maxPoint[2];

            vec3 p1 = minPoint;
            vec3 p2 = vec3(maxX, minY, minZ);
            vec3 p3 = vec3(maxX, minY, maxZ);
            vec3 p4 = vec3(minX, minY, maxZ);
            vec3 p5 = vec3(minX, maxY, minZ);
            vec3 p6 = vec3(maxX, maxY, minZ);
            vec3 p7 = maxPoint;
            vec3 p8 = vec3(minX, maxY, maxZ);

            vec3 v1 = vec3(1, 0, 0);
            vec3 v2 = vec3(0, 1, 0);
            vec3 v3 = vec3(0, 0, 1);
            double width = maxX - minX;
            double height = maxZ - minZ;
            double length = maxY - minY;

            Rectangle* rect1 = new Rectangle(p1 + v1 * width / 2.0 + v3 * height / 2.0, v1, v3, width, height);
            Rectangle* rect2 = new Rectangle(p7 - v1 * width / 2.0 - v3 * height / 2.0, -v1, v3, width, height);

            Rectangle* rect3 = new Rectangle(p1 + v2 * length / 2.0 + v3 * height / 2.0, -v2, v3, length, height);
            Rectangle* rect4 = new Rectangle(p7 - v2 * length / 2.0 - v3 * height / 2.0, v2, v3, length, height);

            Rectangle* rect5 = new Rectangle(p1 + v2 * length / 2.0 + v1 * width / 2.0, v2, v1, length, width);
            Rectangle* rect6 = new Rectangle(p7 - v2 * length / 2.0 - v1 * width / 2.0, -v2, v1, length, width);
            rects = new Object*[numberOfObjects]{rect1, rect2, rect3, rect4, rect5, rect6};
        }

    double intersect(const Ray& ray){
        Hit hit = findClosestHit(ray, rects, numberOfObjects);
        return hit.distance;
    }    
};


void sortByAxis(Object** triangles, int numberOfTriangles, int axis){
    for (int i = 1; i < numberOfTriangles; i++){
        Object* x = triangles[i];
        int j;
        for (j = i; j >= 1; j--){
            if ((x -> computeCentroid())[axis] < (triangles[j-1] -> computeCentroid())[axis]){
                triangles[j] = triangles[j-1];
            }
            else{
                break;
            }
        }
        triangles[j] = x;
    }
}


class Node{
    public:
        int leafSize;
        bool isLeafNode;
        Node* node1;
        Node* node2;
        int depth;
        Object** triangles;
        int numberOfTriangles;
        BoundingBox boundingBox;
        Node(){}
        Node(Object** _triangles, int _numberOfTriangles, int _leafSize=12, int depth=0){
            leafSize = _leafSize;
            boundingBox = BoundingBox(_triangles, _numberOfTriangles);
            if (_numberOfTriangles <= leafSize){
                triangles = _triangles;
                numberOfTriangles = _numberOfTriangles;
                isLeafNode = true;
                return;
            }
            isLeafNode = false;
            int axis = depth % 3;
            sortByAxis(_triangles, _numberOfTriangles, axis);
            int numberOfTriangles1 = _numberOfTriangles / 2;
            int numberOfTriangles2 = _numberOfTriangles - numberOfTriangles1;

            Object** node1Triangles = new Object*[numberOfTriangles1];
            Object** node2Triangles = new Object*[numberOfTriangles2];

            for (int i = 0; i < numberOfTriangles1; i++){
                node1Triangles[i] = _triangles[i];
            }
            for (int i = 0; i < numberOfTriangles2; i++){
                node2Triangles[i] = _triangles[i + numberOfTriangles1];
            }

            node1 = new Node(node1Triangles, numberOfTriangles1, _leafSize, depth+1);
            node2 = new Node(node2Triangles, numberOfTriangles2, _leafSize, depth+1);
        }


    void intersect(const Ray& ray, Hit& hit){
        if (isLeafNode){
            if (numberOfTriangles == 0){
                return;
            }
            Hit closestHit = findClosestHit(ray, triangles, numberOfTriangles);
            if (closestHit.distance > constants::EPSILON && (closestHit.distance < hit.distance || hit.distance == -1)){
                hit.distance = closestHit.distance;
                hit.objectID = closestHit.objectID;
            }
            return;
        }

        double d1 = node1 -> boundingBox.intersect(ray);
        double d2 = node2 -> boundingBox.intersect(ray);

        bool node1Hit = d1 > constants::EPSILON && (d1 < hit.distance || hit.distance == -1);
        bool node2Hit = d2 > constants::EPSILON && (d2 < hit.distance || hit.distance == -1);
        if (node1Hit && node2Hit){
            if (d1 > d2){
                node1 -> intersect(ray, hit);
                node2 -> intersect(ray, hit);
            }
            else{
                node2 -> intersect(ray, hit);
                node1 -> intersect(ray, hit);
            }

        }
        else if (node1Hit){
            node1 -> intersect(ray, hit);
        }
        else if (node2Hit){
            node2 -> intersect(ray, hit);
        }
    }
};


class BoundingVolumeHierarchy{
    public:
        Node* rootNode;
        BoundingVolumeHierarchy(){}
        BoundingVolumeHierarchy(Object** triangles, int numberOfTriangles, int leafSize){
            rootNode = new Node(triangles, numberOfTriangles, leafSize);
        }

        Hit intersect(const Ray& ray){
            Hit hit;
            hit.distance = -1;
            hit.objectID = -1;
            double distanceToBoundingBox = rootNode -> boundingBox.intersect(ray);

            if (distanceToBoundingBox > constants::EPSILON){
                rootNode -> intersect(ray, hit);
            }
            return hit;
        }
};


class ObjectUnion : public Object{
    public:
        Object** objects;
        int numberOfObjects;
        vec3 position;
        double size;
        Sphere boundingSphere;
        BoundingVolumeHierarchy bvh;
        bool useBVH;
        bool containsLight = false;
        ObjectUnion(Object** _objects, int _numberOfObjects, vec3 _position, double _size, bool constructBVH=false) : Object(){
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
            useBVH = constructBVH;
            if (constructBVH){
                bvh = BoundingVolumeHierarchy(_objects, _numberOfObjects, 12);
            }

            for (int i = 0; i < numberOfObjects; i++){
                objects[i] -> objectID = i;
                if (objects[i] -> isLightSource()){
                  containsLight = true;
                }
            }
        }

        ~ObjectUnion(){
            for (int i = 0; i < numberOfObjects; i++){
                delete objects[i];
            }
        }

        bool isLightSource() override{
            return containsLight;
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
            if (useBVH){
                return bvh.intersect(ray);
            }

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
            // TODO: Only sample the index of the objects that are lightsources.
            vec3 randomPoint = objects[sampleRandomObjectIndex()] -> generateRandomSurfacePoint();
            inversePDF = area * areaToAnglePDFFactor(randomPoint, hit);
            return randomPoint;
        }
};


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
    double maxDistance = maximumDistance(averagePosition, vertexArray, numberOfVertices);

    for (int i = 0; i < numberOfVertices; i++){
        vertexArray[i] = ((vertexArray[i] - averagePosition) / maxDistance + desiredCenter) * desiredSize;
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


ObjectUnion* loadObjectModel(std::string fileName, Material* material, const bool allowSmoothShading){
    dataSizes nums = getVertexDataSizes(fileName);

    vec3 vertexArray[nums.numVertices];
    vec3 vertexUVArray[nums.numVertexUVs];
    vec3 vertexNormalArray[nums.numVertexNormals];
    populateVertexArrays(fileName, vertexArray, vertexUVArray, vertexNormalArray);

    double desiredSize = 0.7;
    vec3 desiredCenter = vec3(0, 0.2, 1);
    changeVectors(desiredCenter, desiredSize, vertexArray, nums.numVertices);

    Object** triangles = new Object*[nums.numTriangles];
    populateTriangleArray(fileName, vertexArray, vertexUVArray, vertexNormalArray, triangles, material, allowSmoothShading);
    ObjectUnion* loadedObject = new ObjectUnion(triangles, nums.numTriangles, desiredCenter, desiredSize, true);
    std::cout << nums.numTriangles << "\n";
    return loadedObject;
}

#endif