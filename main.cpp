#include <iostream>
#include <fstream>
#include <cstdint>
#include <vector>
#include <cmath>
#include <chrono>
#include <random>
#include "vec3.h"
#include "Objects.h"
#include "utils.h"
#include "constants.h"


std::shared_ptr<Plane> thisFloor = std::make_shared<Plane>(vec3(0,-0.35,0), vec3(1,0,0), vec3(0,0,-1), std::make_shared<DiffuseMaterial>(WHITE));
std::shared_ptr<Plane> frontWall = std::make_shared<Plane>(vec3(0,0,-0.35), vec3(1,0,0), vec3(0,1,0), std::make_shared<DiffuseMaterial>(WHITE));
std::shared_ptr<Plane> leftWall = std::make_shared<Plane>(vec3(1,0,0), vec3(0,0,1), vec3(0,1,0), std::make_shared<DiffuseMaterial>(RED));
std::shared_ptr<Plane> rightWall = std::make_shared<Plane>(vec3(-1,0,0), vec3(0,1,0), vec3(0,0,1), std::make_shared<DiffuseMaterial>(GREEN));
std::shared_ptr<Plane> roof = std::make_shared<Plane>(vec3(0,1.2,0), vec3(0,0,-1), vec3(1,0,0), std::make_shared<DiffuseMaterial>(WHITE));
std::shared_ptr<Plane> backWall = std::make_shared<Plane>(vec3(0,0,3.5), vec3(0,1,0), vec3(1,0,0), std::make_shared<DiffuseMaterial>(WHITE));

std::shared_ptr<Sphere> ball1 = std::make_shared<Sphere>(vec3(0.35,0,0), 0.35, std::make_shared<ReflectiveMaterial>(BLUE));
std::shared_ptr<Sphere> ball2 = std::make_shared<Sphere>(vec3(-0.45,0,0.6), 0.35, std::make_shared<TransparentMaterial>(WHITE, 0.8, 2));

std::shared_ptr<Rectangle> lightSource = std::make_shared<Rectangle>(vec3(0, 1.199, 1), vec3(0,0,-1), vec3(1,0,0), 1, 1, std::make_shared<DiffuseMaterial>(WHITE, 0.8, 1, 1, WHITE, WARM_WHITE, 10));

std::vector<std::shared_ptr<Object>> objectPtrList = {thisFloor, roof, frontWall, backWall, rightWall, leftWall, ball1, ball2, lightSource};
std::vector<std::shared_ptr<Object>> lightsourcePtrList = {lightSource};
vec3 cameraPosition = vec3(0, 1, 3);
vec3 viewingDirection = vec3(0.0, -0.28734788556, -0.95782628522);
vec3 screenPosition = addVectors(cameraPosition,viewingDirection);
vec3 screenNormalVector = -viewingDirection;
vec3 screenYVector = vec3(0, 0.95782628522, -0.28734788556);
std::random_device rand_dev;
std::minstd_rand generator(rand_dev());
std::uniform_int_distribution<int> distr(0, lightsourcePtrList.size()-1);
// TODO: Add classes etc for camera, screen etc. 


vec3 raytrace(Ray& ray, int depth, vec3& throughput);


vec3 indexToPosition(int x, int y){
    // TODO: Move to a camera/screen class
    double screenWidth = 1.0;
    double screenHeight = screenWidth * (double) constants::HEIGHT / (double) constants::WIDTH; 
    vec3 screenXVector = -crossVectors(screenNormalVector, screenYVector);

    double localXCoordinate = (double) x * screenWidth / (double) constants::WIDTH - (double) screenWidth / (double) 2;
    vec3 localX = multiplyVector(screenXVector, localXCoordinate);

    double localYCoordinate = (double) y * screenHeight / (double) constants::HEIGHT - (double) screenHeight / (double) 2.0;

    vec3 localY = multiplyVector(screenYVector, localYCoordinate);
    vec3 xPlusY = addVectors(localX, localY);
    return addVectors(xPlusY, screenPosition);
}


 vec3 getStartingDirections(int x, int y){
    vec3 pixelVector = indexToPosition(x, y);
    vec3  directionVector = subtractVectors(pixelVector, cameraPosition);
    return normalizeVector(directionVector);
 }


Hit findClosestHit(Ray& ray, std::vector<std::shared_ptr<Object>>& objects){
    Hit closestHit;
    closestHit.distance = -1;

    for (int i = 0; i < objects.size(); i++){
        Hit hit = (*objects[i]).findClosestHit(ray);
        if (hit.distance > constants::EPSILON && (hit.distance < closestHit.distance || closestHit.distance == -1)){
            hit.intersectedObjectIndex = i;
            closestHit = hit;
        }
    }
    if (closestHit.distance < constants::EPSILON){
        return closestHit;
    }

    vec3 scaledDirectionVector = multiplyVector(ray.directionVector, closestHit.distance);
    closestHit.intersectionPoint = addVectors(ray.startingPosition, scaledDirectionVector);
    closestHit.normalVector = (*(objects[closestHit.intersectedObjectIndex])).getNormalVector(closestHit);
    closestHit.incomingVector = ray.directionVector;
    return closestHit;
 }

vec3 indirectLighting(Hit& hit, int depth, vec3 throughput, double randomThreshold){
    brdfData brdfResult = (*objectPtrList[hit.intersectedObjectIndex]).material -> sample(hit);

    throughput = multiplyVectorElementwise(throughput, brdfResult.brdfMultiplier);
    throughput = divideVector(throughput, randomThreshold);
    Ray newRay;
    newRay.startingPosition = hit.intersectionPoint;
    newRay.directionVector = brdfResult.outgoingVector;
    newRay.specular = brdfResult.specular;
    vec3 recursiveColor = raytrace(newRay, depth+1, throughput);
    return multiplyVectorElementwise(recursiveColor, brdfResult.brdfMultiplier);
}

inline int randomLightsourceIndex(){
    return distr(generator);
}

vec3 directLighting(Hit& hit){
    int numLightSources = lightsourcePtrList.size();
    int randomIndex = randomLightsourceIndex();
    int lightIndex;
    for (int i = 0; i < objectPtrList.size(); i++){
        if (objectPtrList[i] == lightsourcePtrList[randomIndex]){
            lightIndex = i;
            break;
        }
    }
    vec3 randomPoint = (*lightsourcePtrList[randomIndex]).generateRandomSurfacePoint();
    Material lightMaterial = *(*lightsourcePtrList[randomIndex]).material; // TODO: Perhaps an issue here.
    double area = (*lightsourcePtrList[randomIndex]).area;

    vec3 vectorTowardsLight = subtractVectors(randomPoint, hit.intersectionPoint);
    double distanceToLight = vectorTowardsLight.length();
    vectorTowardsLight = normalizeVector(vectorTowardsLight);
    Ray lightRay;
    lightRay.startingPosition = hit.intersectionPoint;
    lightRay.directionVector =  vectorTowardsLight;
    Hit lightHit = findClosestHit(lightRay, objectPtrList);

    if (lightHit.intersectedObjectIndex != lightIndex || hit.intersectedObjectIndex == lightIndex){
        return BLACK;
    }

    vec3 lightVector = -vectorTowardsLight;

    double P = dotVectors(lightHit.normalVector, lightVector) / pow(lightHit.distance, 2);
    P = std::max(0.0, P);
    vec3 brdfMultiplier = (*objectPtrList[hit.intersectedObjectIndex]).material -> eval();
    vec3 lightEmmitance = multiplyVector(lightMaterial.emmissionColor, lightMaterial.lightIntensity);
    double cosine = dotVectors(hit.normalVector, vectorTowardsLight);
    cosine = std::max(0.0, cosine);
    vec3 color = multiplyVector(lightEmmitance, P * cosine * area * (double) numLightSources);
    color =  multiplyVectorElementwise(color, brdfMultiplier);
    return color;
}

vec3 raytrace(Ray& ray, int depth, vec3& throughput){
    int forceRecusionLimit = 3;
    Hit rayHit = findClosestHit(ray, objectPtrList);
    vec3 color;
    if (rayHit.distance <= constants::EPSILON){
        return BLACK;
    }

    if (!constants::enableNextEventEstimation || depth == 0 || ray.specular){
        Material objectMaterial = *(*objectPtrList[rayHit.intersectedObjectIndex]).material;
        color = multiplyVector(objectMaterial.emmissionColor, objectMaterial.lightIntensity);
    }

    bool allowRecursion;
    double randomThreshold;

    if (depth < forceRecusionLimit){
        randomThreshold = 1;
        allowRecursion = true;
    }
    else{
        randomThreshold = std::min(throughput.max(), 0.9);
        double randomValue = ((double) rand() / (RAND_MAX));
        allowRecursion = randomValue < randomThreshold;
    }
    
    if (!allowRecursion){
        return color;
    }   

    vec3 indirect = indirectLighting(rayHit, depth, throughput, randomThreshold);
    vec3 direct;

    if (constants::enableNextEventEstimation){
        direct = directLighting(rayHit);
    }
    else{
        direct = BLACK;
    }

    vec3 tempColor = addVectors(indirect, direct);
    tempColor = divideVector(tempColor, randomThreshold);
    color = addVectors(color, tempColor);

    return color;

 }

vec3 computePixelColor(int x, int y){
    vec3 pixelColor = BLACK;
    Ray ray;
    ray.directionVector = getStartingDirections(x, y);
    ray.startingPosition = cameraPosition;
    for (int iter = 0; iter < constants::mcIterations; iter++){
        vec3 throughput = vec3(1,1,1);
        vec3 sampledColor = raytrace(ray, 0, throughput);
        pixelColor = addVectors(pixelColor, sampledColor);
    }
    return divideVector(pixelColor, (double)constants::mcIterations);
}


void printPixelColor(vec3 rgb){
    int r = int(rgb.x() * (double) 255);
    int g = int(rgb.y() * (double) 255);
    int b = int(rgb.z() * (double) 255);
    std::cout << r << ' ' << g << ' ' << b << '\n';
}

int main() {
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    std::cout << "P3\n" << constants::WIDTH << ' ' << constants::HEIGHT << "\n255\n";
    std::vector<double> colorArray(constants::WIDTH * constants::HEIGHT * 3, 0);
    // Write pixel data (BGR format)
    for (int y = constants::HEIGHT-1; y >= 0; y--) {
        for (int x = constants::WIDTH-1; x >= 0; x--) {
            vec3 rgb = computePixelColor(x, y);

            printPixelColor(colorClip(rgb));
        }
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    std::clog << "Time taken: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

    return 0;
}