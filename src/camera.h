#ifndef CAMERA_H
#define CAMERA_H
#include "vec3.h"
#include "constants.h"


class Camera{
    public:
        vec3 position;
        vec3 viewingDirection;
        vec3 yVector;
        vec3 screenPosition;
        vec3 screenXVector;
        double screenWidth;
        double screenHeight;
        Camera(){}
        Camera(vec3 _position=vec3(0,0,0), vec3 _viewingDirection=vec3(0,0,1), vec3 _yVector=vec3(0,1,0)){
            position = _position;
            viewingDirection = normalizeVector(_viewingDirection);
            if (dotVectors(viewingDirection, _yVector) != 0){
                vec3 perpendicularVector = crossVectors(viewingDirection, _yVector);
                _yVector = crossVectors(perpendicularVector, viewingDirection);
            }
            yVector = normalizeVector(_yVector);
            screenWidth = 1.0;
            screenHeight = screenWidth * (double) constants::HEIGHT / (double) constants::WIDTH; 
            screenXVector = crossVectors(viewingDirection, yVector);
            screenPosition = addVectors(position,viewingDirection);
        }

    vec3 indexToPosition(double x, double y){
        double localXCoordinate = x * screenWidth / (double) constants::WIDTH - (double) screenWidth / 2.0;
        vec3 localX = multiplyVector(screenXVector, localXCoordinate);

        double localYCoordinate = y * screenHeight / (double) constants::HEIGHT - (double) screenHeight / 2.0;

        vec3 localY = multiplyVector(yVector, localYCoordinate);
        vec3 xPlusY = addVectors(localX, localY);
        return addVectors(xPlusY, screenPosition);
    }


    vec3 getStartingDirections(double x, double y){
        vec3 pixelVector = indexToPosition(x, y);
        vec3  directionVector = subtractVectors(pixelVector, position);
        return normalizeVector(directionVector);
 }

}; 


#endif