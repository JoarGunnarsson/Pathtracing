#include <iostream>
#include <fstream>
#include <chrono>
#include "vec3.h"
#include "objects.h"
#include "camera.h"
#include "utils.h"
#include "constants.h"


struct Scene{
    Object** objects;
    Camera* camera;
    int numberOfObjects;
};


ValueMap1D* createValueMap1D(std::string fileName, double uMax=1, double vMax=1){
    std::fstream mapFile(fileName, std::ios_base::in); 
    int width;
    int height;
    int dimension;
    mapFile >> width;
    mapFile >> height;
    mapFile >> dimension;
    int N = width * height * dimension;
    double valueHolder;
    double* dataArray = new double[N];
    for (int i = 0; i < N; i++){
        mapFile >> dataArray[i];
    }

    return new ValueMap1D(dataArray, width, height, uMax, vMax);;
}



ValueMap3D* createValueMap3D(std::string fileName, double uMax=1, double vMax=1){
    std::fstream mapFile(fileName, std::ios_base::in); 
    int width;
    int height;
    int dimension;
    mapFile >> width;
    mapFile >> height;
    mapFile >> dimension;
    int N = width * height * dimension;
    double valueHolder;
    double* dataArray = new double[N];
    for (int i = 0; i < N; i++){
        mapFile >> dataArray[i];
    }
    return new ValueMap3D(dataArray, width, height, uMax, vMax);;
}


vec3 raytrace(Ray ray, Object** objects, const int numberOfObjects);


vec3 directLighting(const Hit& hit, Object** objects, const int numberOfObjects){
    int lightSourceIdxArray[numberOfObjects];

    int numberOfLightSources = 0;
    for (int i = 0; i < numberOfObjects; i++){
        if (objects[i] -> material -> isLightSource){
            lightSourceIdxArray[numberOfLightSources] = i;
            numberOfLightSources++;
        }
    }
    int randomIndex = randomInt(0, numberOfLightSources);
    int lightIndex = lightSourceIdxArray[randomIndex];

    double inversePDF;
    vec3 randomPoint = objects[lightIndex] -> randomLightPoint(hit.intersectionPoint, inversePDF);

    vec3 vectorTowardsLight = randomPoint - hit.intersectionPoint;
    double distanceToLight = vectorTowardsLight.length();
    vectorTowardsLight = normalizeVector(vectorTowardsLight);

    Ray lightRay;
    lightRay.startingPosition = hit.intersectionPoint;
    lightRay.directionVector =  vectorTowardsLight;

    Hit lightHit = findClosestHit(lightRay, objects, numberOfObjects);
    bool sameDistance = std::abs(distanceToLight - lightHit.distance) <= constants::EPSILON;
    bool hitFromBehind = dotVectors(vectorTowardsLight, hit.normalVector) < 0.0;
    if (lightHit.intersectedObjectIndex != lightIndex || hit.intersectedObjectIndex == lightIndex || !sameDistance || hitFromBehind){
        return BLACK;
    }

    vec3 brdfMultiplier = objects[hit.intersectedObjectIndex] -> eval(hit.intersectionPoint);
    vec3 lightEmitance = objects[lightIndex] -> getLightEmittance(lightHit.intersectionPoint);

    double cosine = dotVectors(hit.normalVector, vectorTowardsLight);
    cosine = std::max(0.0, cosine);

    return brdfMultiplier * cosine * lightEmitance * inversePDF * (double) numberOfLightSources;
}


vec3 raytrace(Ray ray, Object** objects, const int numberOfObjects){
    vec3 color = vec3(0,0,0);
    vec3 brdfAccumulator = vec3(1,1,1);
    vec3 throughput = vec3(1,1,1);
    int forceRecusionLimit = 3;
    double randomThreshold = 1;
    for (int depth = 0; depth <= constants::maxRecursionDepth; depth++){
        Hit rayHit = findClosestHit(ray, objects, numberOfObjects);
        if (rayHit.distance <= constants::EPSILON){
            break;
        }
        
        bool specularRay = ray.type == REFLECTED || ray.type == TRANSMITTED;
        if (!constants::enableNextEventEstimation || depth == 0 || specularRay){
            vec3 lightEmitance = objects[rayHit.intersectedObjectIndex] -> getLightEmittance(rayHit.intersectionPoint);
            color += lightEmitance * brdfAccumulator * (dotVectors(ray.directionVector, rayHit.normalVector) < 0);
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
            break;
        }   

        if (constants::enableNextEventEstimation){
            vec3 direct = directLighting(rayHit, objects, numberOfObjects);
            color += direct * brdfAccumulator / randomThreshold;
        }

        brdfData brdfResult = objects[rayHit.intersectedObjectIndex] -> sample(rayHit, objects, numberOfObjects);

        throughput *= brdfResult.brdfMultiplier / randomThreshold;
        ray.startingPosition = rayHit.intersectionPoint;
        ray.directionVector = brdfResult.outgoingVector;
        ray.type = brdfResult.type;

        brdfAccumulator *= brdfResult.brdfMultiplier / randomThreshold;
    }
    return color;

 }


vec3 computePixelColor(const int x, const int y, Scene scene){
    vec3 pixelColor = vec3(0,0,0);
    Ray ray;
    ray.startingPosition = scene.camera -> position;
    for (int i = 0; i < constants::samplesPerPixel; i++){
        double x_offset = randomNormal() / 2.0;
        double y_offset = randomNormal() / 2.0;
        ray.directionVector = scene.camera -> getStartingDirections(x + x_offset, y + y_offset);
        vec3 sampledColor = raytrace(ray, scene.objects, scene.numberOfObjects);
        pixelColor += sampledColor;    
    }

    return pixelColor / (double) constants::samplesPerPixel;
}


void printPixelColor(vec3 rgb, std::ofstream& file){
    int r = int(rgb[0] * (double) 255);
    int g = int(rgb[1] * (double) 255);
    int b = int(rgb[2] * (double) 255);
    file << r << ' ' << g << ' ' << b << '\n';
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


Scene createScene(){
    ValueMap3D* whiteMap = new ValueMap3D(WHITE*0.8);
    ValueMap3D* warmWhiteMap = new ValueMap3D(WARM_WHITE);
    ValueMap3D* pureWhiteMap = new ValueMap3D(WHITE);
    ValueMap3D* redMap = new ValueMap3D(RED*0.8);
    ValueMap3D* greenMap = new ValueMap3D(GREEN*0.8);
    ValueMap3D* blueMap = new ValueMap3D(BLUE*0.8);
    double* zero = new double(0);
    double* one = new double(1);
    double* ten = new double(6);
    ValueMap1D* zeroMap = new ValueMap1D(zero);
    ValueMap1D* oneMap = new ValueMap1D(one);
    ValueMap1D* tenMap = new ValueMap1D(ten);

    ValueMap3D* worldMap = createValueMap3D("./maps/world.map");

    ValueMap1D* worldRoughnessMap = createValueMap1D("./maps/world_roughness.map");
    
    MaterialData defaultMaterialData;
    defaultMaterialData.albedoMap = whiteMap;
    DiffuseMaterial* whiteDiffuseMaterial = new DiffuseMaterial(defaultMaterialData);
    ReflectiveMaterial* whiteReflectiveMaterial = new ReflectiveMaterial(defaultMaterialData);   

    MaterialData stripedData;
    double* stripes = new double[6]{0.8, 0.8, 0.8, 0, 0, 0};
    ValueMap3D* stripedMap = new ValueMap3D(stripes, 2, 1, 0.1, 1);
    stripedData.albedoMap = stripedMap;
    DiffuseMaterial* stripedMaterial = new DiffuseMaterial(stripedData);

    MaterialData redMaterialData;
    redMaterialData.albedoMap = redMap;
    DiffuseMaterial* redDiffuseMaterial = new DiffuseMaterial(redMaterialData);

    MaterialData greenMaterialData;
    greenMaterialData.albedoMap = greenMap;
    DiffuseMaterial* greenDiffuseMaterial = new DiffuseMaterial(greenMaterialData);

    MaterialData blueMaterialData;
    blueMaterialData.albedoMap = blueMap;
    blueMaterialData.isDielectric = false;
    ReflectiveMaterial* blueReflectiveMaterial = new ReflectiveMaterial(blueMaterialData);

    MaterialData ball2Data;
    ball2Data.albedoMap = blueMap;
    ball2Data.refractiveIndex = 2;
    ball2Data.attenuationCoefficient = 0.1;
    ball2Data.roughnessMap = worldRoughnessMap;
    ball2Data.absorptionAlbedo = WHITE - BLUE;
    ball2Data.isDielectric = true;
    FrostyMaterial* ball2Material = new FrostyMaterial(ball2Data);

    MaterialData lightMaterialData;
    lightMaterialData.albedoMap = blueMap;
    lightMaterialData.emmissionColorMap = warmWhiteMap;
    lightMaterialData.lightIntensityMap = tenMap;
    lightMaterialData.isLightSource = true;
    DiffuseMaterial* lightSourceMaterial = new DiffuseMaterial(lightMaterialData);

    MaterialData glassData;
    glassData.albedoMap = pureWhiteMap;
    glassData.refractiveIndex = 1.5;
    TransparentMaterial* pane1Material = new TransparentMaterial(glassData);

    MaterialData frostyGlassData;
    frostyGlassData.albedoMap = pureWhiteMap;
    frostyGlassData.refractiveIndex = 1.5;
    frostyGlassData.roughnessMap = worldRoughnessMap;
    FrostyMaterial* pane2Material = new FrostyMaterial(frostyGlassData);

    Plane* thisFloor = new Plane(vec3(0,-0.35,0), vec3(-1,0,0), vec3(0,0,1), whiteDiffuseMaterial);
    Plane*  frontWall = new Plane(vec3(0,0,-0.35), vec3(-1,0,0), vec3(0,-1,0), stripedMaterial);
    Plane* leftWall = new Plane(vec3(1,0,0), vec3(0,0,1), vec3(0,1,0), redDiffuseMaterial);
    Plane* rightWall = new Plane(vec3(-1,0,0), vec3(0,1,0), vec3(0,0,1), greenDiffuseMaterial);
    Plane* roof = new Plane(vec3(0,1.2,0), vec3(1,0,0), vec3(0,0,1), whiteDiffuseMaterial);
    Plane* backWall = new Plane(vec3(0,0,3.5), vec3(0,1,0), vec3(1,0,0), whiteDiffuseMaterial);

    Rectangle* frontPane1 = new Rectangle(vec3(0.25,0.5,1.2), vec3(1,0,0), vec3(0,1,0), 0.5, 0.5, pane1Material);
    Rectangle* BackPane1 = new Rectangle(vec3(0.25,0.5,1.15), vec3(-1,0,0), vec3(0,1,0), 0.5, 0.5, pane1Material);

    Rectangle* frontPane2 = new Rectangle(vec3(-0.25,0.5,1.2), vec3(1,0,0), vec3(0,1,0), 0.5, 0.5, pane2Material);
    Rectangle* BackPane2 = new Rectangle(vec3(-0.25,0.5,1.15), vec3(-1,0,0), vec3(0,1,0), 0.5, 0.5, pane2Material);

    Sphere* ball1 = new Sphere(vec3(0.35,0,0), 0.35, blueReflectiveMaterial);
    Sphere* ball2 = new Sphere(vec3(-0.45, 0, 0.6), 0.35, ball2Material);

    Rectangle* lightSource = new Rectangle(vec3(0, 1.199, 1), vec3(0,0,-1), vec3(1,0,0), 1, 1, lightSourceMaterial);

    Object** objects = new Object*[9]{thisFloor, roof, frontWall, backWall, rightWall, leftWall, ball1, ball2, lightSource};

    vec3 cameraPosition = vec3(0, 0.7, 3.4);
    vec3 viewingDirection = vec3(0.0, -0.2, -1);
    vec3 screenYVector = vec3(0, 1, 0);
    Camera* camera = new Camera(cameraPosition, viewingDirection, screenYVector);
    int numberOfObjects = 9;
    Scene scene;
    scene.objects = objects;
    scene.camera = camera;
    scene.numberOfObjects = numberOfObjects;
    return scene;
}


int main() {
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    std::ofstream dataFile;
    dataFile.open("./temp/result_data.txt");
    
    dataFile << "SIZE:" << constants::WIDTH << ' ' << constants::HEIGHT << "\n";

    Scene scene = createScene();

    int maxIters = 0;
    for (int y = constants::HEIGHT-1; y >= 0; y--) {
        double progress = double(constants::HEIGHT - y) / (double) constants::HEIGHT;
        printProgress(progress);
        for (int x = constants::WIDTH-1; x >= 0; x--) {
            vec3 rgb = computePixelColor(x, y, scene);
            printPixelColor(colorClip(rgb), dataFile);
        }
    }
    dataFile.close();
    std::clog << std::endl;
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    std::clog << "Time taken: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

    //TODO: remember to free the scene.
    return 0;
}