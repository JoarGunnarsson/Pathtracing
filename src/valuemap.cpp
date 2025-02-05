#include "valuemap.h"


ValueMap::ValueMap(const int _data, const int _width, const int _height, const double _u_max, const double _v_max){
    double* _data_ptr = new double[1]{double(_data)};
    initialise(_data_ptr, _width, _height, _u_max, _v_max);
}

ValueMap::ValueMap(const double _data, const int _width, const int _height, const double _uMax, const double _v_max){
    double* _data_ptr = new double[1]{_data};
    initialise(_data_ptr, _width, _height, _uMax, _v_max);
}

ValueMap::ValueMap(double* _data, const int _width, const int _height, const double _u_max, const double _v_max){
    initialise(_data, _width, _height, _u_max, _v_max);
}

ValueMap::~ValueMap(){
    delete[] data;
}

void ValueMap::initialise(double* _data, const int _width, const int _height, const double _u_max, const double _v_max){
    data = _data;
    width = _width;
    height = _height;
    u_max = _u_max;
    v_max = _v_max;
}


double ValueMap1D::get(const double u, const double v) const {
    if (isnan(u) || isnan(v)){
        return 0;
    }
    int u_idx = int((double) width * pos_fmod(u / u_max, 1.0));
    int v_idx = int((double) height * (pos_fmod((v_max - v) / v_max, 1.0)));
    int index = (v_idx * width + u_idx);
    return data[index];
}

vec3 ValueMap3D::get(const double u, const double v) const {
    if (isnan(u) || isnan(v)){
        return vec3(0,0,0);
    }
    int u_idx = int((double) width * pos_fmod(u / u_max, 1.0));
    int v_idx = int((double) height * pos_fmod(v / v_max, 1.0));
    int start_index = 3 * (v_idx * width + u_idx);
    return vec3(data[start_index], data[start_index + 1], data[start_index + 2]);
}


ValueMap1D* create_value_map_1D(const char* file_name, double u_max, double v_max) {
    FILE* map_file = fopen(file_name, "r");
    if (!map_file) {
        return nullptr;
    }

    int width, height, dimension;
    if (fscanf(map_file, "%d %d %d", &width, &height, &dimension) != 3) {
        fclose(map_file);
        return nullptr;
    }

    int N = width * height * dimension;
    double* data_array = new double[N];

    for (int i = 0; i < N; i++) {
        fscanf(map_file, "%lf", &data_array[i]);
    }

    fclose(map_file);
    return new ValueMap1D(data_array, width, height, u_max, v_max);
}


ValueMap3D* create_value_map_3D(const char* file_name, double u_max, double v_max) {
    FILE* map_file = fopen(file_name, "r");
    if (!map_file) {
        return nullptr;
    }

    int width, height, dimension;
    if (fscanf(map_file, "%d %d %d", &width, &height, &dimension) != 3) {
        fclose(map_file);
        return nullptr;
    }

    int N = width * height * dimension;
    double* data_array = new double[N];

    for (int i = 0; i < N; i++) {
        fscanf(map_file, "%lf", &data_array[i]);
    }

    fclose(map_file);
    return new ValueMap3D(data_array, width, height, u_max, v_max);
}
