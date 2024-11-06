#ifndef VALUEMAP_H
#define VALUEMAP_H
#include "vec3.h"
#include "utils.h"
#include <cstdio>  // For file handling functions

class ValueMap{
    public:
        double* data;
        int width;
        double u_max;
        int height;
        double v_max;
        ValueMap(){}

        ValueMap(const int _data, const int _width=1, const int _height=1, const double _u_max=1, const double _v_max=1){
            double* _data_ptr = new double[1]{double(_data)};
            initialise(_data_ptr, _width, _height, _u_max, _v_max);
        }

        ValueMap(const double _data, const int _width=1, const int _height=1, const double _uMax=1, const double _v_max=1){
            double* _data_ptr = new double[1]{_data};
            initialise(_data_ptr, _width, _height, _uMax, _v_max);
        }

        ValueMap(double* _data, const int _width=1, const int _height=1, const double _u_max=1, const double _v_max=1){
            initialise(_data, _width, _height, _u_max, _v_max);
        }
        
        ~ValueMap(){
            delete[] data;
        }
        
    private: 
        void initialise(double* _data, const int _width, const int _height, const double _u_max, const double _v_max){
            data = _data;
            width = _width;
            height = _height;
            u_max = _u_max;
            v_max = _v_max;
        }
};


class ValueMap1D : public ValueMap{
    public:
        using ValueMap::ValueMap;

    double get(const double u, const double v) {
        if (isnan(u) or isnan(v)){
            return 0;
        }
        int u_idx = int((double) width * pos_fmod(u / u_max, 1.0));
        int v_idx = int((double) height * (pos_fmod((v_max - v) / v_max, 1.0)));
        int index = (v_idx * width + u_idx);
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
        int u_idx = int((double) width * pos_fmod(u / u_max, 1.0));
        int v_idx = int((double) height * pos_fmod(v / v_max, 1.0));
        int start_index = 3 * (v_idx * width + u_idx);
        return vec3(data[start_index], data[start_index + 1], data[start_index + 2]);
    }
};


ValueMap1D* create_value_map_1D(const char* fileName, double u_max = 1, double v_max = 1) {
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
    return new ValueMap1D(dataArray, width, height, u_max, v_max);
}


ValueMap3D* create_value_map_3D(const char* fileName, double u_max = 1, double v_max = 1) {
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
    return new ValueMap3D(dataArray, width, height, u_max, v_max);
}


#endif