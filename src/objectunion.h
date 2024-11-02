#ifndef OBJECTUNION_H
#define OBJECTUNION_H

#include <fstream>
#include "constants.h"
#include "vec3.h"
#include "utils.h"
#include "objects.h"


vec3 getMaxPoint(Object** triangles, int numberOfTriangles){
    if (numberOfTriangles == 0){
        return vec3(0,0,0);
    }
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
    if (numberOfTriangles == 0){
        return vec3(0,0,0);
    }
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
        vec3 p1;
        vec3 p2;
        double width;
        double height;
        double length;
        double axisLength[3];
        BoundingBox(){}
        BoundingBox(Object** _triangles, int numberOfTriangles){
            p1 = getMinPoint(_triangles, numberOfTriangles);
            p2 = getMaxPoint(_triangles, numberOfTriangles);

            width = p2[0] - p1[0];
            length = p2[1] - p1[1];
            height = p2[2] - p1[2];
            axisLength[0] = width;
            axisLength[1] = length;
            axisLength[2] = height;
        }
    
    inline bool isWithinBounds(const double x, const double lower, const double higher){
        return lower <= x && x <= higher;
    }

    double intersect(const Ray& ray){
        double t[6];
        bool insideBounds[6];
        for (int i = 0; i < 3; i++){
            if (std::abs(ray.directionVector[i]) < constants::EPSILON){
                t[i] = -1;
                insideBounds[i] = false;
                continue;
            }
            t[i] = (p1[i] - ray.startingPosition[i]) / ray.directionVector[i];
            vec3 hitPoint = ray.directionVector * t[i] + ray.startingPosition;
            vec3 differenceVector = hitPoint - p1;
            if (i == 0){
                insideBounds[i] = isWithinBounds(differenceVector[1], 0, length) && isWithinBounds(differenceVector[2], 0, height);
            }
            else if (i ==  1){
                insideBounds[i] = isWithinBounds(differenceVector[0], 0, width) && isWithinBounds(differenceVector[2], 0, height);
            }
            else if (i == 2){
                insideBounds[i] = isWithinBounds(differenceVector[0], 0, width) && isWithinBounds(differenceVector[1], 0, length);
            }
        }

        for (int i = 0; i < 3; i++){
            if (std::abs(ray.directionVector[i]) < constants::EPSILON){
                t[i+3] = -1;
                insideBounds[i+3] = false;
                continue;
            }
            t[i+3] = (p2[i] - ray.startingPosition[i]) / ray.directionVector[i];
            vec3 hitPoint = ray.directionVector * t[i+3] + ray.startingPosition;
            vec3 differenceVector = hitPoint - p2;
            if (i == 0){
                insideBounds[i+3] = isWithinBounds(differenceVector[1], -length, 0) && isWithinBounds(differenceVector[2], -height, 0);
            }
            else if (i ==  1){
                insideBounds[i+3] = isWithinBounds(differenceVector[0], -width, 0) && isWithinBounds(differenceVector[2], -height, 0);
            }
            else if (i == 2){
                insideBounds[i+3] = isWithinBounds(differenceVector[0], -width, 0) && isWithinBounds(differenceVector[1], -length, 0);
            }
        }

        double minT = -1;
        for (int i = 0; i < 6; i++){
            if (insideBounds[i] && (minT == -1 || minT > t[i]) && t[i] > constants::EPSILON){
                minT = t[i];
            }
        }

        return minT;
    }    
};


void sortByAxis(Object** triangles, int numberOfTriangles, int axis){
    std::sort(triangles, triangles + numberOfTriangles, [axis](Object* obj1, Object* obj2){ 
        return (obj1 -> computeCentroid())[axis] < (obj2 -> computeCentroid())[axis]; 
        });
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
            int axis = getSplitAxis();

            sortByAxis(_triangles, _numberOfTriangles, axis);
            int splitIndex = _numberOfTriangles / 2;

            Object** node1Triangles = new Object*[splitIndex];
            Object** node2Triangles = new Object*[_numberOfTriangles - splitIndex];

            for (int i = 0; i < splitIndex; i++){
                node1Triangles[i] = _triangles[i];
            }
            for (int i = splitIndex; i < _numberOfTriangles; i++){
                node2Triangles[i - splitIndex] = _triangles[i];
            }

            node1 = new Node(node1Triangles, splitIndex, _leafSize, depth+1);
            node2 = new Node(node2Triangles, _numberOfTriangles - splitIndex, _leafSize, depth+1);
        }

    int getSplitAxis(){
        int axis;
        double maxLength = 0;
        for (int i = 0; i < 3; i++){
            if (boundingBox.axisLength[i] >= maxLength){
                axis = i;
                maxLength = boundingBox.axisLength[i];
            }
        }
        return axis;
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
                if (d2 < hit.distance || hit.distance == -1){
                    node2 -> intersect(ray, hit);
                }
            }
            else{
                node2 -> intersect(ray, hit);
                if (d1 < hit.distance || hit.distance == -1){
                    node1 -> intersect(ray, hit);
                }
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
            double distanceToBoundingBox = rootNode -> boundingBox.intersect(ray);

            Hit hit;
            hit.distance = -1;
            hit.objectID = -1;
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
        double* cumulativeArea;
        int* lightSourceConversionIndices;
        int numberOfLightSources;
        BoundingVolumeHierarchy bvh;
        bool useBVH;
        bool containsLightSource = false;
        ObjectUnion(Object** _objects, const int _numberOfObjects, const bool constructBVH=false) : Object(){
            objects = _objects;
            numberOfObjects = _numberOfObjects;

            area = 0;
            for (int i = 0; i < numberOfObjects; i++){
                area += objects[i] -> area;
                if (objects[i] -> isLightSource()){
                    numberOfLightSources++;
                }
            }

            cumulativeArea = new double[numberOfObjects];
            lightSourceConversionIndices = new int[numberOfLightSources];
            int j = 0;
            for (int i = 0; i < numberOfObjects; i++){
                if (!objects[i] -> isLightSource()){
                    continue;
                }
                cumulativeArea[j] = objects[i] -> area;
                if (j != 0){
                    cumulativeArea[j] += cumulativeArea[j-1];
                }
                lightSourceConversionIndices[j] = i;
                j++;
            }

            useBVH = constructBVH;
            if (constructBVH){
                bvh = BoundingVolumeHierarchy(_objects, _numberOfObjects, 12);
            }

            for (int i = 0; i < numberOfObjects; i++){
                objects[i] -> objectID = i;
                if (objects[i] -> isLightSource()){
                  containsLightSource = true;
                }
            }
        }

        ~ObjectUnion() override{
            for (int i = 0; i < numberOfObjects; i++){
                if (objects[i] -> alive){delete objects[i];}
            }

            delete[] objects;
            delete[] cumulativeArea;
            delete[] lightSourceConversionIndices;
        }

        bool isLightSource() override{
            return containsLightSource;
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

            Hit hit = findClosestHit(ray, objects, numberOfObjects);
            hit.objectID = hit.intersectedObjectIndex;
            hit.intersectedObjectIndex = -1;
            return hit;
        }
        
        vec3 getNormalVector(const vec3& surfacePoint, const int objectID) override{
            return objects[objectID] -> getNormalVector(surfacePoint, objectID);
        }

        int sampleRandomObjectIndex(){
            double randomAreaSplit = randomUniform(0, area);
            int max = numberOfLightSources - 1;
            int min = 0;
            int index;

            if (cumulativeArea[0] >= randomAreaSplit){
                return lightSourceConversionIndices[0];
            }

            while (min <= max){
                index = (max - min) / 2 + min;

                if (cumulativeArea[index] < randomAreaSplit){
                    min = index + 1;
                }
                else if (cumulativeArea[index] == randomAreaSplit || (cumulativeArea[index] >= randomAreaSplit && cumulativeArea[index-1] < randomAreaSplit)){
                    break;
                }
                else{
                    max = index - 1;
                }
            }

            return lightSourceConversionIndices[index];
        }

        vec3 generateRandomSurfacePoint() override{
            return objects[sampleRandomObjectIndex()] -> generateRandomSurfacePoint();
        }

        vec3 randomLightPoint(const vec3& intersectionPoint, double& inversePDF) override{
            int randomIndex = sampleRandomObjectIndex();
            vec3 randomPoint = objects[randomIndex] -> generateRandomSurfacePoint();
            inversePDF = cumulativeArea[numberOfLightSources-1] * areaToAnglePDFFactor(randomPoint, intersectionPoint, randomIndex);
            return randomPoint;
        }
};


struct dataSizes{
    int numVertices = 0;
    int numVertexUVs = 0;
    int numVertexNormals = 0;
    int numTriangles = 0;
};


int numberOfCharOccurances(const std::string& line, const char character){
    int count = 0;
    for (int i = 0; i < line.length(); i++){
        if (line[i] == character){
            count++;
        }
    }
    return count;
}


std::string getNthWord(const std::string& line, const char delimiter, const int n){
    int numberOfWords = numberOfCharOccurances(line, delimiter);

    if (n > numberOfWords){
        return "";
    }

    int start = 0;
    int end = 0;
    for (int i = 0; i < n+1; i++){
        start = end;
        end = line.find(delimiter, start)+1;
    }
    return line.substr(start, end - start - 1);
}


dataSizes getVertexDataSizes(const std::string& fileName){
    std::ifstream modelFile(fileName);
    std::string line;
    dataSizes nums;
    
    while(std::getline(modelFile, line)){
        std::string firstWord = getNthWord(line, ' ', 0);

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


void populateVertexArrays(const std::string& fileName, vec3* vertexArray, vec3* vertexUVArray, vec3* vertexNormalArray){
    int vertexIdx = 0;
    int vertexUVIdx = 0;
    int vertexNormalIdx = 0;

    std::ifstream modelFile(fileName);
    std::string line;
    while(std::getline(modelFile, line)){
        std::string firstWord = getNthWord(line, ' ', 0);
        bool isVertex = firstWord == "v";
        bool isVertexUV = firstWord == "vt";
        bool isVertexNormal = firstWord == "vn";

        if (isVertex){
            double v1 = std::stod(getNthWord(line, ' ', 1));
            double v2 = std::stod(getNthWord(line, ' ', 2));
            double v3 = std::stod(getNthWord(line, ' ', 3));
            vertexArray[vertexIdx] = vec3(v1, v2, v3);
            vertexIdx++;
        }
        else if (isVertexUV){
            double u = std::stod(getNthWord(line, ' ', 1));
            double v = std::stod(getNthWord(line, ' ', 2));
            vertexUVArray[vertexUVIdx] = vec3(u, v, 0);
            vertexUVIdx++;
        }
        else if (isVertexNormal){
            double n1 = std::stod(getNthWord(line, ' ', 1));
            double n2 = std::stod(getNthWord(line, ' ', 2));
            double n3 = std::stod(getNthWord(line, ' ', 3));
            vertexNormalArray[vertexNormalIdx] = vec3(n1, n2, n3);
            vertexNormalIdx++;
        }
    }
}


vec3 computeAveragePosition(const vec3* vertexArray, const int numberOfVertices){
    vec3 avg = vec3(0,0,0);
    for (int i = 0; i < numberOfVertices; i++){
        avg += vertexArray[i];
    }
    return avg / numberOfVertices;
}


double maximumDistance(const vec3& center, const vec3* vertexArray, const int numberOfVertices){
    double maxDistance = 0;
    for (int i = 0; i < numberOfVertices; i++){
        double distance = (vertexArray[i] - center).length();
        if (distance > maxDistance){
            maxDistance = distance;
        }
    }
    return maxDistance;
}

void changeVectors(const vec3& desiredCenter, const double desiredSize, vec3* vertexArray, const int numberOfVertices){
    vec3 averagePosition = computeAveragePosition(vertexArray, numberOfVertices);
    double maxDistance = maximumDistance(averagePosition, vertexArray, numberOfVertices);

    for (int i = 0; i < numberOfVertices; i++){
        vertexArray[i] = ((vertexArray[i] - averagePosition) / maxDistance + desiredCenter) * desiredSize;
    }
}


struct populateVertexVectorData{
    const std::string vertexData;
    vec3 v;
    bool vSuccess = false;
    vec3 uv;
    bool uvSuccess = false;
    vec3 n;
    bool nSuccess = false;
    const vec3* vertexArray;
    const vec3* vertexUVArray;
    const vec3* vertexNormalArray;

    populateVertexVectorData(const std::string& data, const vec3* vertexArray, const vec3* vertexUVArray, const vec3* vertexNormalArray)
        : vertexData(data), vertexArray(vertexArray), vertexUVArray(vertexUVArray), vertexNormalArray(vertexNormalArray) {}
};


void populateVertexVectors(populateVertexVectorData& args){
    std::string vIdx = getNthWord(args.vertexData, '/', 0);
    std::string UVIdx = getNthWord(args.vertexData, '/', 1);
    std::string nIdx = getNthWord(args.vertexData, '/', 2);
    
    if (vIdx != ""){
        args.v = args.vertexArray[std::stoi(vIdx)-1];
        args.vSuccess = true;
    }

    if (UVIdx != ""){
        args.uv = args.vertexUVArray[std::stoi(UVIdx)-1];
        args.uvSuccess = true;
    }

    if (nIdx != ""){
        args.n = args.vertexNormalArray[std::stoi(nIdx)-1];
        args.nSuccess = true;
    }
}


struct TriangleConstructionArgs{
    const std::string triangleData;
    const int idx1;
    const int idx2;
    const int idx3;
    Material* material;
    const vec3* vertexArray;
    const vec3* vertexUVArray;
    const vec3* vertexNormalArray;
    const bool enableSmoothShading;

    TriangleConstructionArgs(const std::string& data, const int idx1, const int idx2, const int idx3, Material* material, const vec3* vertexArray,
    const vec3* vertexUVArray, const vec3* vertexNormalArray, const bool enableSmoothShading) : triangleData(data), idx1(idx1), idx2(idx2),
    idx3(idx3), material(material), vertexArray(vertexArray), vertexUVArray(vertexUVArray), vertexNormalArray(vertexNormalArray), enableSmoothShading(enableSmoothShading) {}
};


struct TriangleCreationResult{
    Triangle* triangle;
    bool success = false;
};


TriangleCreationResult constructTriangle(TriangleConstructionArgs& args){
    std::string v1Data = getNthWord(args.triangleData, ' ', args.idx1);
    populateVertexVectorData data1 = populateVertexVectorData(v1Data, args.vertexArray, args.vertexUVArray, args.vertexNormalArray);
    populateVertexVectors(data1);

    std::string v2Data = getNthWord(args.triangleData, ' ', args.idx2);
    populateVertexVectorData data2 = populateVertexVectorData(v2Data, args.vertexArray, args.vertexUVArray, args.vertexNormalArray);
    populateVertexVectors(data2);

    std::string v3Data = getNthWord(args.triangleData, ' ', args.idx3);
    populateVertexVectorData data3 = populateVertexVectorData(v3Data, args.vertexArray, args.vertexUVArray, args.vertexNormalArray);
    populateVertexVectors(data3);

    bool loadedVerticesSuccessfully = data1.vSuccess && data2.vSuccess && data3.vSuccess;
    
    TriangleCreationResult result;

    if (!loadedVerticesSuccessfully){
        return result;
    }

    Triangle* triangle = new Triangle(data1.v, data2.v, data3.v, args.material);

    bool loadedUVSuccessfully = data1.uvSuccess && data2.uvSuccess && data3.uvSuccess;
    if (loadedUVSuccessfully){triangle -> setVertexUV(data1.uv, data2.uv, data3.uv);}

    bool loadedNormalsSuccessfully  = data1.nSuccess && data2.nSuccess && data3.nSuccess;
    if (loadedNormalsSuccessfully && args.enableSmoothShading){
        triangle -> setVertexNormals(data1.n, data2.n, data3.n);
    }
    result.success = true;
    result.triangle = triangle;
    return result;
}


int populateTriangleArray(std::string fileName, vec3* vertexArray, vec3* vertexUVArray, vec3* vertexNormalArray, Object** triangleArray, Material* material, const bool enableSmoothShading){
    std::ifstream modelFile(fileName);
    std::string line;
    int shapeIdx = 0;
    while(std::getline(modelFile, line)){
        std::string firstWord = getNthWord(line, ' ', 0);
        
        bool isShape = firstWord == "f";
        if (!isShape){
            continue;
        }
        int numberOfSpaces = numberOfCharOccurances(line, ' ');
        bool isTriangle = numberOfSpaces == 3;
        bool isQuad = numberOfSpaces == 4;
        if (isTriangle){
            TriangleConstructionArgs args = TriangleConstructionArgs(line, 1, 2, 3, material, vertexArray, vertexUVArray, vertexNormalArray, enableSmoothShading);
            TriangleCreationResult result = constructTriangle(args);
            if (!result.success){
                continue;
            }
            triangleArray[shapeIdx] = result.triangle;
            shapeIdx++;
        }

        else if (isQuad){
            TriangleConstructionArgs args1 = TriangleConstructionArgs(line, 1, 2, 3, material, vertexArray, vertexUVArray, vertexNormalArray, enableSmoothShading);
            TriangleCreationResult result1 = constructTriangle(args1);
            if (result1.success){
                triangleArray[shapeIdx] = result1.triangle;
                shapeIdx++;
            }
            
            TriangleConstructionArgs args2 = TriangleConstructionArgs(line, 1, 3, 4, material, vertexArray, vertexUVArray, vertexNormalArray, enableSmoothShading);
            TriangleCreationResult result2 = constructTriangle(args2);
            if (result1.success){
                triangleArray[shapeIdx] = result2.triangle;
                shapeIdx++;
            }
        }
    }
    return shapeIdx;
}


ObjectUnion* loadObjectModel(std::string fileName, Material* material, const bool enableSmoothShading){
    dataSizes nums = getVertexDataSizes(fileName);

    vec3 vertexArray[nums.numVertices];
    vec3 vertexUVArray[nums.numVertexUVs];
    vec3 vertexNormalArray[nums.numVertexNormals];
    populateVertexArrays(fileName, vertexArray, vertexUVArray, vertexNormalArray);

    double desiredSize = 0.7;
    vec3 desiredCenter = vec3(0, 0.2, 1);
    changeVectors(desiredCenter, desiredSize, vertexArray, nums.numVertices);

    Object** triangles = new Object*[nums.numTriangles];
    int validTriangles = populateTriangleArray(fileName, vertexArray, vertexUVArray, vertexNormalArray, triangles, material, enableSmoothShading);
    ObjectUnion* loadedObject = new ObjectUnion(triangles, validTriangles, true);
    return loadedObject;
}

#endif