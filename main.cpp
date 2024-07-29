#include <iostream>
#include <fstream>
#include <cstdint>
#include <vector>
#include <cmath>
#include <chrono>
#include "vec3.h"
#include "Sphere.h"
#include "utils.h"


int WIDTH = 1000;
int HEIGHT = 1000;
std::shared_ptr<Sphere> s1 = std::make_shared<Sphere>(vec3(5,1,-1), 1, Material(BLUE));
std::shared_ptr<Sphere> s2 = std::make_shared<Sphere>(vec3(4,1,1), 0.5, Material(RED));
std::shared_ptr<Sphere> s3 = std::make_shared<Sphere>(vec3(4, 4, 4), 1, Material(GREEN));
std::vector<std::shared_ptr<Object>> objectPtrList = {s1, s2, s3};

// TODO: Add classes etc for camera, screen etc. 


vec3 indexToPosition(int x, int y){
    // TODO: Move to a camera/screen class
    float screenWidth = 1;
    float screenHeight = screenWidth * (float) HEIGHT / (float) WIDTH; 
    vec3 screenNormalVector = vec3(-1, 0, 0);
    vec3 screenYVector = vec3(0, 1, 0);
    vec3 screenXVector = crossVectors(screenNormalVector, screenYVector);

    vec3 screenPosition = vec3(1, 0, 0);

    float localXCoordinate = (float) x * screenWidth / (float) WIDTH - (float) screenWidth / (float) 2;
    vec3 localX = multiplyVector(screenXVector, localXCoordinate);

    float localYCoordinate = (float) y * screenHeight / (float) HEIGHT - (float) screenHeight / (float) 2.0;
    vec3 localY = multiplyVector(screenYVector, localYCoordinate);
    return addVectors(addVectors(localX, localY), screenPosition);
}


 vec3 getStartingDirections(int x, int y){
    vec3 pixelVector = indexToPosition(x, y);
    vec3 cameraPosition = vec3(0, 0, 0);
    vec3  directionVector = subtractVectors(pixelVector, cameraPosition);
    return normalizeVector(directionVector);
 }


Hit findClosestDistance(Ray ray, std::vector<std::shared_ptr<Object>> objects){
    // TODO: Ask object for the actual intersected object. Is better for object unions.
    Hit closestHit;
    closestHit.distance = -1;

    for (int i = 0; i < objects.size(); i++){
        Hit hit = (*objects[i]).findClosestHit(ray);
        if (hit.distance > 0 && (hit.distance < closestHit.distance || closestHit.distance == -1)){
            hit.intersectedObjectIndex = i;
            closestHit = hit;
        }
    }
    if (closestHit.distance < 0){
        return closestHit;
    }

    closestHit.intersectionPoint = addVectors(ray.startingPosition, multiplyVector(ray.directionVector, closestHit.distance));
    closestHit.normalVector = (*(objects[closestHit.intersectedObjectIndex])).getNormalVector(closestHit);
    return closestHit;
 }

vec3 raytrace(Ray ray){
    Hit rayHit = findClosestDistance(ray, objectPtrList);
    if (rayHit.distance <= 0){
        return BLACK;
    }

    //displayVector((*objectPtrList[0]).material.diffuseColor);
    return (*objectPtrList[rayHit.intersectedObjectIndex]).material.diffuseColor;

 }

vec3 computePixelColor(int x, int y){

    int mcIterations = 1;
    vec3 pixelColor = vec3(0, 0, 0);
    vec3 cameraPosition = vec3(0, 0, 0);
    Ray ray;
    ray.directionVector = getStartingDirections(x, y);
    ray.startingPosition = cameraPosition;

    for (int iter = 0; iter < mcIterations; iter++){
        vec3 sampledColor = raytrace(ray);
        pixelColor = addVectors(pixelColor, sampledColor);
    }
    return divideVector(pixelColor, mcIterations);
}


void printPixelColor(vec3 rgb){
    int r = int(rgb.x() * (float) 255);
    int g = int(rgb.y() * (float) 255);
    int b = int(rgb.z() * (float) 255);
    std::cout << r << ' ' << g << ' ' << b << '\n';
}


int main() {
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    std::cout << "P3\n" << WIDTH << ' ' << HEIGHT << "\n255\n";
    std::vector<float> colorArray(WIDTH * HEIGHT * 3, 0);
    // Write pixel data (BGR format)
    for (int y = HEIGHT-1; y >= 0; y--) {
        for (int x = WIDTH-1; x >= 0; x--) {
            vec3 rgb = computePixelColor(x, y);
            printPixelColor(rgb);
        }
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    std::clog << "Time taken: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

    return 0;
}