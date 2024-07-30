#include <iostream>
#include <fstream>
#include <cstdint>
#include <vector>
#include <cmath>
#include <chrono>
#include "vec3.h"
#include "Sphere.h"
#include "utils.h"

const double EPSILON = 0.0001;
int mcIterations = 100;
int WIDTH = 1000;
int HEIGHT = 1000;
std::shared_ptr<Sphere> thisFloor = std::make_shared<Sphere>(vec3(0,-0.1-100000,0), 100000, Material(WHITE));
std::shared_ptr<Sphere> roof = std::make_shared<Sphere>(vec3(0,0.35 + 100000,0), 100000, Material(WHITE));
std::shared_ptr<Sphere> leftWall = std::make_shared<Sphere>(vec3(0,0,-0.35 -100000), 100000, Material(RED));
std::shared_ptr<Sphere> rightWall = std::make_shared<Sphere>(vec3(0,0,0.35 + 100000), 100000, Material(GREEN));
std::shared_ptr<Sphere> backWall = std::make_shared<Sphere>(vec3(1 + 100000,0,0), 100000, Material(WHITE));
std::shared_ptr<Sphere> frontWall = std::make_shared<Sphere>(vec3(-2 -100000,0,0), 100000, Material(WHITE));
std::shared_ptr<Sphere> blueBall = std::make_shared<Sphere>(vec3(1,0.1,0), 0.1, Material(WHITE));
std::shared_ptr<Sphere> lightSource = std::make_shared<Sphere>(vec3(1, 0.35, 0), 0.1, Material(WHITE, 1, 1, 1, WHITE, WHITE, 10));
std::vector<std::shared_ptr<Object>> objectPtrList = {thisFloor, roof, leftWall, rightWall, backWall, frontWall, blueBall, lightSource};
vec3 cameraPosition = vec3(0, 0.1, 0);
vec3 screenPosition = vec3(1, 0.1, 0);
vec3 screenNormalVector = vec3(-1, 0, 0);
vec3 screenYVector = vec3(0, 1, 0);
// TODO: Add classes etc for camera, screen etc. 


vec3 raytrace(Ray& ray, int depth, vec3 throughput);


vec3 indexToPosition(int x, int y){
    // TODO: Move to a camera/screen class
    double screenWidth = 1;
    double screenHeight = screenWidth * (double) HEIGHT / (double) WIDTH; 
    vec3 screenXVector = crossVectors(screenNormalVector, screenYVector);

    double localXCoordinate = (double) x * screenWidth / (double) WIDTH - (double) screenWidth / (double) 2;
    vec3 localX = multiplyVector(screenXVector, localXCoordinate);

    double localYCoordinate = (double) y * screenHeight / (double) HEIGHT - (double) screenHeight / (double) 2.0;

    vec3 localY = multiplyVector(screenYVector, localYCoordinate);
    vec3 xPlusY = addVectors(localX, localY);
    return addVectors(xPlusY, screenPosition);
}


 vec3 getStartingDirections(int x, int y){
    vec3 pixelVector = indexToPosition(x, y);
    vec3  directionVector = subtractVectors(pixelVector, cameraPosition);
    return normalizeVector(directionVector);
 }


Hit findClosestDistance(Ray& ray, std::vector<std::shared_ptr<Object>>& objects){
    // TODO: Ask object for the actual intersected object. Is better for object unions.
    Hit closestHit;
    closestHit.distance = -1;

    for (int i = 0; i < objects.size(); i++){
        Hit hit = (*objects[i]).findClosestHit(ray);
        if (hit.distance > EPSILON && (hit.distance < closestHit.distance || closestHit.distance == -1)){
            hit.intersectedObjectIndex = i;
            closestHit = hit;
        }
    }
    if (closestHit.distance < EPSILON){
        return closestHit;
    }

    vec3 scaledDirectionVector = multiplyVector(ray.directionVector, closestHit.distance);
    closestHit.intersectionPoint = addVectors(ray.startingPosition, scaledDirectionVector);
    closestHit.normalVector = (*(objects[closestHit.intersectedObjectIndex])).getNormalVector(closestHit);
    return closestHit;
 }

vec3 indirectLighting(Ray& ray, Hit& hit, int depth, vec3 throughput, double randomThreshold){
    Material objectMaterial = (*objectPtrList[hit.intersectedObjectIndex]).material;
    brdfData brdfResult = objectMaterial.sample(hit);

    throughput = multiplyVectorElementwise(throughput, brdfResult.brdfMultiplier);
    throughput = divideVector(throughput, randomThreshold);
    Ray newRay;
    newRay.startingPosition = hit.intersectionPoint;
    newRay.directionVector = brdfResult.outgoingVector;
    vec3 recursiveColor = raytrace(newRay, depth+1, throughput);
    return multiplyVectorElementwise(recursiveColor, brdfResult.brdfMultiplier);
}

vec3 raytrace(Ray& ray, int depth, vec3 throughput){
    int forceRecusionLimit = 3;
    Hit rayHit = findClosestDistance(ray, objectPtrList);
    if (rayHit.distance <= EPSILON){
        return BLACK;
    }

    Material objectMaterial = (*objectPtrList[rayHit.intersectedObjectIndex]).material;

    vec3 color = multiplyVector(objectMaterial.emmissionColor, objectMaterial.lightIntensity);

    bool allowRecusion;
    double randomThreshold;

    if (depth < forceRecusionLimit){
        randomThreshold = 1;
        allowRecusion = true;
    }
    else{
        randomThreshold = std::min(throughput.max(), 0.9);
        double randomValue = ((double) rand() / (RAND_MAX));
        allowRecusion = randomValue < randomThreshold;
    }
    
    if (!allowRecusion){
        return color;
    }   

    vec3 indirect = indirectLighting(ray, rayHit, depth, throughput, randomThreshold);
    vec3 tempColor = divideVector(indirect, randomThreshold);
    color = addVectors(color, tempColor);
    return color;

 }

vec3 computePixelColor(int x, int y){
    vec3 pixelColor = BLACK;
    Ray ray;
    ray.directionVector = getStartingDirections(x, y);
    ray.startingPosition = cameraPosition;
    for (int iter = 0; iter < mcIterations; iter++){
        vec3 throughput = vec3(1,1,1);
        vec3 sampledColor = raytrace(ray, 0, throughput);
        pixelColor = addVectors(pixelColor, sampledColor);
    }
    return divideVector(pixelColor, (double)mcIterations);
}


void printPixelColor(vec3 rgb){
    int r = int(rgb.x() * (double) 255);
    int g = int(rgb.y() * (double) 255);
    int b = int(rgb.z() * (double) 255);
    std::cout << r << ' ' << g << ' ' << b << '\n';
}


int main() {
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    std::cout << "P3\n" << WIDTH << ' ' << HEIGHT << "\n255\n";
    std::vector<double> colorArray(WIDTH * HEIGHT * 3, 0);
    // Write pixel data (BGR format)
    for (int y = HEIGHT-1; y >= 0; y--) {
        for (int x = WIDTH-1; x >= 0; x--) {
            vec3 rgb = computePixelColor(x, y);

            printPixelColor(colorClip(rgb));
        }
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    std::clog << "Time taken: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

    return 0;
}