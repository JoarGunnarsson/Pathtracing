#ifndef VALUEMAP_H
#define VALUEMAP_H
#include "vec3.h"
#include "utils.h"
#include <cstdio>  // For file handling functions

class ValueMap{
    public:
        double* data;
        int width;
        double uMax;
        int height;
        double vMax;
        ValueMap(){}

        ValueMap(const int _data, const int _width=1, const int _height=1, const double _uMax=1, const double _vMax=1){
            double* _data_ptr = new double[1]{double(_data)};
            initialise(_data_ptr, _width, _height, _uMax, _vMax);
        }

        ValueMap(const double _data, const int _width=1, const int _height=1, const double _uMax=1, const double _vMax=1){
            double* _data_ptr = new double[1]{_data};
            initialise(_data_ptr, _width, _height, _uMax, _vMax);
        }

        ValueMap(double* _data, const int _width=1, const int _height=1, const double _uMax=1, const double _vMax=1){
            initialise(_data, _width, _height, _uMax, _vMax);
        }
        
        ~ValueMap(){
            delete[] data;
        }
        
    private: 
        void initialise(double* _data, const int _width, const int _height, const double _uMax, const double _vMax){
            data = _data;
            width = _width;
            height = _height;
            uMax = _uMax;
            vMax = _vMax;
        }
};


class ValueMap1D : public ValueMap{
    public:
        using ValueMap::ValueMap;

    double get(const double u, const double v) {
        if (isnan(u) or isnan(v)){
            return 0;
        }
        int uIdx = int((double) width * posFmod(u / uMax, 1.0));
        int vIdx = int((double) height * (posFmod((vMax - v) / vMax, 1.0)));
        int index = (vIdx * width + uIdx);
        return data[index];
    }
};


class ValueMap3D : public ValueMap{
    public:
        using ValueMap::ValueMap;

    vec3 get(const double u, const double v){
        if (isnan(u) or isnan(v)){
            return vec3(0,0,0);
        }
        int uIdx = int((double) width * posFmod(u / uMax, 1.0));
        int vIdx = int((double) height * posFmod(v / vMax, 1.0));
        int startIndex = 3 * (vIdx * width + uIdx);
        return vec3(data[startIndex], data[startIndex + 1], data[startIndex + 2]);
    }
};


ValueMap1D* createValueMap1D(const char* fileName, double uMax = 1, double vMax = 1) {
    FILE* mapFile = fopen(fileName, "r");
    if (!mapFile) {
        return nullptr;
    }

    int width, height, dimension;
    if (fscanf(mapFile, "%d %d %d", &width, &height, &dimension) != 3) {
        fclose(mapFile);
        return nullptr;
    }

    int N = width * height * dimension;
    double* dataArray = new double[N];

    for (int i = 0; i < N; i++) {
        fscanf(mapFile, "%lf", &dataArray[i]);
    }

    fclose(mapFile);
    return new ValueMap1D(dataArray, width, height, uMax, vMax);
}


ValueMap3D* createValueMap3D(const char* fileName, double uMax = 1, double vMax = 1) {
    FILE* mapFile = fopen(fileName, "r");
    if (!mapFile) {
        return nullptr;
    }

    int width, height, dimension;
    if (fscanf(mapFile, "%d %d %d", &width, &height, &dimension) != 3) {
        fclose(mapFile);
        return nullptr;
    }

    int N = width * height * dimension;
    double* dataArray = new double[N];

    for (int i = 0; i < N; i++) {
        fscanf(mapFile, "%lf", &dataArray[i]);
    }

    fclose(mapFile);
    return new ValueMap3D(dataArray, width, height, uMax, vMax);
}


#endif