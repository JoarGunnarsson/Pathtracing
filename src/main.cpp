#include <iostream>
#include <fstream>
#include <chrono>
#include "vec3.h"
#include "objects.h"
#include "camera.h"
#include "utils.h"
#include "constants.h"
#include "objectunion.h"


struct Scene{
    Object** objects;
    Camera* camera;
    int numberOfObjects;
    MaterialManager* material_manager;
};


vec3 directLighting(Hit& hit, Object** objects, const int numberOfObjects){
    int lightSourceIdxArray[numberOfObjects];

    int numberOfLightSources = 0;
    
    for (int i = 0; i < numberOfObjects; i++){
        if (objects[i] -> isLightSource()){
            lightSourceIdxArray[numberOfLightSources] = i;
            numberOfLightSources++;
        }
    }
    if (numberOfLightSources == 0){
        return BLACK;
    }
    
    int randomIndex = randomInt(0, numberOfLightSources);
    int lightIndex = lightSourceIdxArray[randomIndex];

    double inversePDF;
    vec3 randomPoint = objects[lightIndex] -> randomLightPoint(hit.intersectionPoint, inversePDF);

    vec3 vectorTowardsLight = randomPoint - hit.intersectionPoint;
    double distanceToLight = vectorTowardsLight.length();
    vectorTowardsLight = normalizeVector(vectorTowardsLight);
    hit.outgoingVector = vectorTowardsLight;
    
    Ray lightRay;
    lightRay.startingPosition = hit.intersectionPoint;
    lightRay.directionVector =  vectorTowardsLight;

    Hit lightHit = findClosestHit(lightRay, objects, numberOfObjects);
    bool inShadow = lightHit.intersectedObjectIndex != lightIndex;
    bool sameDistance = std::abs(distanceToLight - lightHit.distance) <= constants::EPSILON;
    bool hitFromBehind = dotVectors(vectorTowardsLight, hit.normalVector) < 0.0;
    bool insideObject = dotVectors(hit.incomingVector, hit.normalVector) > 0.0;
    if ( inShadow || !sameDistance || hitFromBehind || insideObject){
        return BLACK;
    }

    vec3 brdfMultiplier = objects[hit.intersectedObjectIndex] -> eval(hit);
    vec3 lightEmitance = objects[lightIndex] -> getLightEmittance(lightHit);

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
        
        bool isSpecularRay = ray.type == REFLECTED || ray.type == TRANSMITTED;

        if (!constants::enableNextEventEstimation || depth == 0 || isSpecularRay){
            vec3 lightEmitance = objects[rayHit.intersectedObjectIndex] -> getLightEmittance(rayHit);
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
        int barWidth = 60;

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
    /*
    ValueMap3D* worldMap = createValueMap3D("./maps/world.map");
    ValueMap3D* sakuraMap = createValueMap3D("./maps/sakura.map");
    ValueMap3D* templeMap = createValueMap3D("./maps/temple.map");
    ValueMap3D* cobbleMap = createValueMap3D("./maps/cobblestone.map");
    ValueMap1D* worldRoughnessMap = createValueMap1D("./maps/world_roughness.map");
    
    */

    MaterialManager* manager = new MaterialManager();
    MaterialData whiteData;
    whiteData.albedoMap = new ValueMap3D(WHITE * 0.8);
    DiffuseMaterial* whiteDiffuseMaterial = new DiffuseMaterial(whiteData);
    manager -> add_material(whiteDiffuseMaterial);
    
    MaterialData whiteReflectiveData;
    whiteReflectiveData.albedoMap = new ValueMap3D(WHITE * 0.8);
    ReflectiveMaterial* whiteReflectiveMaterial = new ReflectiveMaterial(whiteReflectiveData);   
    manager -> add_material(whiteReflectiveMaterial);

    MaterialData stripedData;
    double* stripes = new double[6]{0.8, 0.8, 0.8, 0, 0, 0};
    ValueMap3D* stripedMap = new ValueMap3D(stripes, 2, 1, 0.1, 1);
    stripedData.albedoMap = stripedMap;
    DiffuseMaterial* stripedMaterial = new DiffuseMaterial(stripedData);
    manager -> add_material(stripedMaterial);

    MaterialData redMaterialData;
    redMaterialData.albedoMap = new ValueMap3D(RED);
    DiffuseMaterial* redDiffuseMaterial = new DiffuseMaterial(redMaterialData);
    manager -> add_material(redDiffuseMaterial);

    MaterialData greenMaterialData;
    greenMaterialData.albedoMap = new ValueMap3D(GREEN);
    DiffuseMaterial* greenDiffuseMaterial = new DiffuseMaterial(greenMaterialData);
    manager -> add_material(greenDiffuseMaterial);

    MaterialData goldData;
    goldData.albedoMap = new ValueMap3D(GOLD);
    goldData.roughnessMap = new ValueMap1D(0.1);
    goldData.refractiveIndex = 0.277;
    goldData.extinctionCoefficient = 2.92;
    goldData.isDielectric = false;
    MicrofacetMaterial* goldMaterial = new MicrofacetMaterial(goldData);
    manager -> add_material(goldMaterial);

    MaterialData lightMaterialData;
    lightMaterialData.albedoMap =  new ValueMap3D(WHITE * 0.8);
    lightMaterialData.emmissionColorMap =  new ValueMap3D(WARM_WHITE);
    lightMaterialData.lightIntensityMap = new ValueMap1D(10.0);
    lightMaterialData.isLightSource = true;
    DiffuseMaterial* lightSourceMaterial = new DiffuseMaterial(lightMaterialData);
    manager -> add_material(lightSourceMaterial);

    /*
    MaterialData glassData;
    glassData.refractiveIndex = 1.5;
    TransparentMaterial* pane1Material = new TransparentMaterial(glassData);

    MaterialData frostyGlassData;
    frostyGlassData.albedoMap = pureWhiteMap;
    frostyGlassData.refractiveIndex = 1.5;
    frostyGlassData.roughnessMap = new ValueMap1D(0.1);
    frostyGlassData.percentageDiffuseMap = new ValueMap1D(0.0);
    MicrofacetMaterial* pane2Material = new MicrofacetMaterial(frostyGlassData);

    MaterialData sakuraData;
    sakuraData.albedoMap = sakuraMap;
    DiffuseMaterial* sakuraMaterial = new DiffuseMaterial(sakuraData);

    MaterialData templeData;
    templeData.albedoMap = templeMap;
    DiffuseMaterial* templeMaterial = new DiffuseMaterial(templeData);

    MaterialData cobbleData;
    cobbleData.albedoMap = cobbleMap;
    DiffuseMaterial* cobbleMaterial = new DiffuseMaterial(cobbleData);
    */
   
    MaterialData modelData;
    modelData.albedoMap = createValueMap3D("./maps/bunny.map", 1, -1);
    DiffuseMaterial* modelMaterial = new DiffuseMaterial(modelData);
    manager -> add_material(modelMaterial);
    
    Plane* thisFloor = new Plane(vec3(0,-0.35,0), vec3(1,0,0), vec3(0,0,-1), whiteDiffuseMaterial);
    Rectangle* frontWall = new Rectangle(vec3(0,0.425,-0.35), vec3(1,0,0), vec3(0,1,0), 2, 1.55, whiteDiffuseMaterial);
    Rectangle* leftWall = new Rectangle(vec3(-1,0.425,1.575), vec3(0,0,-1), vec3(0,1,0), 3.85, 1.55, redDiffuseMaterial);
    Rectangle* rightWall = new Rectangle(vec3(1,0.425,1.575), vec3(0,0,-1), vec3(0,-1,0), 3.85, 1.55, greenDiffuseMaterial);
    Plane* roof = new Plane(vec3(0,1.2,0), vec3(1,0,0), vec3(0,0,1), whiteDiffuseMaterial);
    Rectangle* backWall = new Rectangle(vec3(0,0.425,3.5), vec3(0,1,0), vec3(1,0,0), 2, 3.85/2.0, whiteDiffuseMaterial);

    /*
    Rectangle* frontPane1 = new Rectangle(vec3(-0.25,0.5,1.2), vec3(1,0,0), vec3(0,1,0), 0.5, 0.5, pane1Material);
    Rectangle* backPane1 = new Rectangle(vec3(-0.25,0.5,1.15), vec3(-1,0,0), vec3(0,1,0), 0.5, 0.5, pane1Material);
    Object** pane1Objects = new Object*[2]{frontPane1, backPane1};
    ObjectUnion* pane1 = new ObjectUnion(pane1Objects, 2);

    Rectangle* frontPane2 = new Rectangle(vec3(0.25,0.5,1.2), vec3(1,0,0), vec3(0,1,0), 0.5, 0.5, pane2Material);
    Rectangle* backPane2 = new Rectangle(vec3(0.25,0.5,1.15), vec3(-1,0,0), vec3(0,1,0), 0.5, 0.5, pane2Material);
    Object** pane2Objects = new Object*[2]{frontPane2, backPane2};
    ObjectUnion* pane2 = new ObjectUnion(pane2Objects, 2);
    */

    //Sphere* ball1 = new Sphere(vec3(-0.35,0,0), 0.35, greenDiffuseMaterial);
    //Sphere* ball2 = new Sphere(vec3(0.45, 0, 0.6), 0.35, lightSourceMaterial);

    Rectangle* lightSource = new Rectangle(vec3(0, 1.199, 1), vec3(0,0,-1), vec3(1,0,0), 1, 1, lightSourceMaterial);
    
    ObjectUnion* loadedModel = loadObjectModel("./models/dragon.obj", goldMaterial, true);

    int numberOfObjects = 8;
    Object** objects = new Object*[numberOfObjects]{thisFloor, frontWall, leftWall, rightWall, roof, backWall, lightSource, loadedModel};

    vec3 cameraPosition = vec3(0, 1, 3);
    vec3 viewingDirection = vec3(0.0, -0.3, -1);
    vec3 screenYVector = vec3(0, 1, 0);
    Camera* camera = new Camera(cameraPosition, viewingDirection, screenYVector);

    Scene scene;
    scene.objects = objects;
    scene.camera = camera;
    scene.numberOfObjects = numberOfObjects;
    scene.material_manager = manager;
    return scene;
}


void clearScene(Scene& scene){
    for (int i = 0; i < scene.numberOfObjects; i++){
        delete scene.objects[i];
    }

    delete[] scene.objects;
    delete scene.material_manager;
    delete scene.camera;
}


int main() {

    std::ofstream dataFile;
    dataFile.open("./temp/result_data.txt");
    
    dataFile << "SIZE:" << constants::WIDTH << ' ' << constants::HEIGHT << "\n";

    Scene scene = createScene();
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    int maxIters = 0;
    for (int y = constants::HEIGHT-1; y >= 0; y--) {
        double progress = double(constants::HEIGHT - y) / (double) constants::HEIGHT;
        printProgress(progress);
        for (int x = 0; x < constants::WIDTH; x++) {
            vec3 rgb = computePixelColor(x, y, scene);
            printPixelColor(colorClip(rgb), dataFile);
        }
    }
    dataFile.close();
    std::clog << std::endl;
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    std::clog << "Time taken: " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]" << std::endl;

    //TODO: remember to free the scene.
    clearScene(scene);
    return 0;
}