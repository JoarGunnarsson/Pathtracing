#include <iostream>
#include <fstream>
#include <cstdint>
#include <vector>
#include <cmath>
#include <chrono>
#include <random>
#include "vec3.h"
#include "objects.h"
#include "camera.h"
#include "utils.h"
#include "constants.h"


DiffuseMaterial whiteDiffuseMaterial = DiffuseMaterial(WHITE);
DiffuseMaterial redDiffuseMaterial = DiffuseMaterial(RED);
DiffuseMaterial greenDiffuseMaterial = DiffuseMaterial(GREEN);
ReflectiveMaterial blueReflectiveMaterial = ReflectiveMaterial(BLUE);
DiffuseMaterial whiteTransparentMaterial = DiffuseMaterial(WHITE, 0.8, 2);
DiffuseMaterial lightSourceMaterial = DiffuseMaterial(WHITE, 0.8, 1, 1, WHITE, WARM_WHITE, 10);

Plane thisFloor = Plane(vec3(0,-0.35,0), vec3(1,0,0), vec3(0,0,-1), &whiteDiffuseMaterial);
Plane  frontWall = Plane(vec3(0,0,-0.35), vec3(1,0,0), vec3(0,1,0), &whiteDiffuseMaterial);
Plane leftWall = Plane(vec3(1,0,0), vec3(0,0,1), vec3(0,1,0), &redDiffuseMaterial);
Plane rightWall = Plane(vec3(-1,0,0), vec3(0,1,0), vec3(0,0,1), &greenDiffuseMaterial);
Plane roof = Plane(vec3(0,1.2,0), vec3(0,0,-1), vec3(1,0,0), &whiteDiffuseMaterial);
Plane backWall = Plane(vec3(0,0,3.5), vec3(0,1,0), vec3(1,0,0), &whiteDiffuseMaterial);

Sphere ball1 = Sphere(vec3(0.35,0,0), 0.35, &blueReflectiveMaterial);
Sphere ball2 = Sphere(vec3(-0.45,0,0.6), 0.35, &whiteTransparentMaterial);

Rectangle lightSource = Rectangle(vec3(0, 1.199, 1), vec3(0,0,-1), vec3(1,0,0), 1, 1, &lightSourceMaterial);

Object* objectPtrList[] = {&thisFloor, &roof, &frontWall, &backWall, &rightWall, &leftWall, &ball1, &ball2, &lightSource};
Object* lightsourcePtrList[] = {&lightSource};

vec3 cameraPosition = vec3(0, 1, 3);
vec3 viewingDirection = vec3(0.0, -0.3, -1);
vec3 screenYVector = vec3(0, 1, 0);
Camera camera = Camera(cameraPosition, viewingDirection, screenYVector);


vec3 raytrace(Ray& ray, int depth, vec3& throughput);


vec3 indirectLighting(Hit& hit, int depth, vec3 throughput, double randomThreshold){
    brdfData brdfResult = (*objectPtrList[hit.intersectedObjectIndex]).material -> sample(hit);

    throughput *= brdfResult.brdfMultiplier / randomThreshold;
    Ray newRay;
    newRay.startingPosition = hit.intersectionPoint;
    newRay.directionVector = brdfResult.outgoingVector;
    newRay.specular = brdfResult.specular;
    vec3 recursiveColor = raytrace(newRay, depth+1, throughput);
    return recursiveColor * brdfResult.brdfMultiplier;
}


vec3 directLighting(Hit& hit){
    int numLightSources = sizeof(lightsourcePtrList) / sizeof(Object*); // lightsourcePtrList.size();
    int randomIndex = randomInt(0, numLightSources);
    int lightIndex;
    for (int i = 0; i < sizeof(objectPtrList) / sizeof(Object*); i++){
        if (objectPtrList[i] == lightsourcePtrList[randomIndex]){
            lightIndex = i;
            break;
        }
    }
    vec3 randomPoint = (*lightsourcePtrList[randomIndex]).generateRandomSurfacePoint();
    Material lightMaterial = *(*lightsourcePtrList[randomIndex]).material; // TODO: Perhaps an issue here.
    double area = (*lightsourcePtrList[randomIndex]).area;

    vec3 vectorTowardsLight = randomPoint - hit.intersectionPoint;
    double distanceToLight = vectorTowardsLight.length();
    vectorTowardsLight = normalizeVector(vectorTowardsLight);
    Ray lightRay;
    lightRay.startingPosition = hit.intersectionPoint;
    lightRay.directionVector =  vectorTowardsLight;
    int numberOfObjects = sizeof(objectPtrList) / sizeof(Object*);
    Hit lightHit = findClosestHit(lightRay, objectPtrList, numberOfObjects);

    if (lightHit.intersectedObjectIndex != lightIndex || hit.intersectedObjectIndex == lightIndex){
        return BLACK;
    }

    vec3 lightVector = -vectorTowardsLight;

    double P = dotVectors(lightHit.normalVector, lightVector) / pow(lightHit.distance, 2);
    P = std::max(0.0, P);
    vec3 brdfMultiplier = (*objectPtrList[hit.intersectedObjectIndex]).material -> eval();
    vec3 lightEmmitance = lightMaterial.emmissionColor * lightMaterial.lightIntensity;
    double cosine = dotVectors(hit.normalVector, vectorTowardsLight);
    cosine = std::max(0.0, cosine);
    vec3 color = brdfMultiplier * lightEmmitance * P * cosine * area * (double) numLightSources;
    return color;
}


vec3 raytrace(Ray& ray, int depth, vec3& throughput){
    if (depth > constants::maxRecursionDepth){
        return BLACK;
    }
    int forceRecusionLimit = 3;
    int numberOfObjects = sizeof(objectPtrList) / sizeof(Object*);
    Hit rayHit = findClosestHit(ray, objectPtrList, numberOfObjects);
    vec3 color;
    if (rayHit.distance <= constants::EPSILON){
        return BLACK;
    }

    if (!constants::enableNextEventEstimation || depth == 0 || ray.specular){
        Material objectMaterial = *(*objectPtrList[rayHit.intersectedObjectIndex]).material;
        color = objectMaterial.emmissionColor * objectMaterial.lightIntensity;
    }

    bool allowRecursion;
    double randomThreshold;

    if (depth < forceRecusionLimit){
        randomThreshold = 1;
        allowRecursion = true;
    }
    else{
        randomThreshold = std::min(throughput.max(), 0.9);
        double randomValue = randomUniform(0, 1);
        allowRecursion = randomValue < randomThreshold;
    }
    
    if (!allowRecursion){
        return BLACK;
        
    }   

    vec3 indirect = indirectLighting(rayHit, depth, throughput, randomThreshold);

    vec3 direct;

    if (constants::enableNextEventEstimation){
        direct = directLighting(rayHit);
    }
    else{
        direct = BLACK;
    }

    color += (indirect + direct) / randomThreshold;

    return color;

 }


vec3 computePixelColor(int x, int y){
    vec3 pixelColor = BLACK;
    Ray ray;
    ray.startingPosition = cameraPosition;
    for (int i = 0; i < constants::mcIterations; i++){
        double x_offset = randomNormal();
        double y_offset = randomNormal();
        ray.directionVector = camera.getStartingDirections(x + x_offset, y + y_offset);
        vec3 throughput = vec3(1,1,1);
        vec3 sampledColor = raytrace(ray, 0, throughput);
        pixelColor = pixelColor + sampledColor;
            
    }
    return pixelColor / (double) constants::mcIterations;
}


void printPixelColor(vec3 rgb){
    int r = int(rgb[0] * (double) 255);
    int g = int(rgb[1] * (double) 255);
    int b = int(rgb[2] * (double) 255);
    std::cout << r << ' ' << g << ' ' << b << '\n';
}


void printProgress(double progress){
    if (progress < 1.0) {
        int barWidth = 70;

        std::clog << "[";
        int pos = barWidth * progress;
        for (int i = 0; i < barWidth; ++i) {
            if (i < pos) std::clog << "=";
            else if (i == pos) std::clog << ">";
            else std::clog << " ";
        }
        std::clog << "] " << int(progress * 100.0) << " %\r";
    }
}


int main() {
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    std::vector<double> colorArray(constants::WIDTH * constants::HEIGHT * 3, 0);
    std::cout << "SIZE:" << constants::WIDTH << ' ' << constants::HEIGHT << "\n";
    for (int y = constants::HEIGHT-1; y >= 0; y--) {
        double progress = double(constants::HEIGHT - y) / (double) constants::HEIGHT;
        printProgress(progress);
        for (int x = constants::WIDTH-1; x >= 0; x--) {
            vec3 rgb = computePixelColor(x, y);
            printPixelColor(colorClip(rgb));
        }
    }
    std::clog << std::endl;
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    std::clog << "Time taken: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

    return 0;
}