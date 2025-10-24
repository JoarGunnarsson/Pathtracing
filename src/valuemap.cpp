#include "valuemap.h"

ValueMap::ValueMap(const int _data, const int _width, const int _height, const double _u_max, const double _v_max) {
    double* _data_ptr = new double[1]{double(_data)};
    initialise(_data_ptr, _width, _height, _u_max, _v_max);
}

ValueMap::ValueMap(const double _data, const int _width, const int _height, const double _uMax, const double _v_max) {
    double* _data_ptr = new double[1]{_data};
    initialise(_data_ptr, _width, _height, _uMax, _v_max);
}

ValueMap::ValueMap(double* _data, const int _width, const int _height, const double _u_max, const double _v_max) {
    initialise(_data, _width, _height, _u_max, _v_max);
}

ValueMap::~ValueMap() {
    delete[] data;
}

void ValueMap::initialise(double* _data, const int _width, const int _height, const double _u_max,
                          const double _v_max) {
    data = _data;
    width = _width;
    height = _height;
    u_max = _u_max;
    v_max = _v_max;
}

double ValueMap1D::get(const double u, const double v) const {
    if (isnan(u) || isnan(v)) {
        return 0;
    }
    int u_idx = int((double) width * pos_fmod(u / u_max, 1.0));
    int v_idx = int((double) height * (pos_fmod((v_max - v) / v_max, 1.0)));
    int index = (v_idx * width + u_idx);
    return data[index];
}

vec3 ValueMap3D::get(const double u, const double v) const {
    if (isnan(u) || isnan(v)) {
        return vec3(0, 0, 0);
    }
    int u_idx = int((double) width * pos_fmod(u / u_max, 1.0));
    int v_idx = int((double) height * pos_fmod(v / v_max, 1.0));
    int start_index = 3 * (v_idx * width + u_idx);
    return vec3(data[start_index], data[start_index + 1], data[start_index + 2]);
}

template<typename ValueMapType>
ValueMapType* create_value_map(const std::string& file_name, double u_max, double v_max) {
    std::ifstream file(file_name, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Could not open file '" + file_name + "'");
    }

    file.seekg(0, std::ios::end);
    std::streamsize size = file.tellg();
    file.seekg(0, std::ios::beg);

    size_t N = size / sizeof(double);

    double* data = new double[N];
    if (!file.read(reinterpret_cast<char*>(data), size)) {
        throw std::runtime_error("Error reading file '" + file_name + "'");
    }

    int width = data[0];
    int height = data[1];
    int dimension = data[2];
    if (width * height * dimension != N - 3) {
        throw std::runtime_error("File '" + file_name + "' does not follow expected format and cannot be loaded");
    }
    double* data_array = new double[N - 3];
    for (size_t i = 3; i < N; ++i) {
        data_array[i - 3] = data[i];
    }

    delete[] data;
    return new ValueMapType(data_array, width, height, u_max, v_max);
}

ValueMap1D* create_value_map_1D(const std::string& file_name, double u_max, double v_max) {
    return create_value_map<ValueMap1D>(file_name, u_max, v_max);
}

ValueMap3D* create_value_map_3D(const std::string& file_name, double u_max, double v_max) {
    return create_value_map<ValueMap3D>(file_name, u_max, v_max);
}
