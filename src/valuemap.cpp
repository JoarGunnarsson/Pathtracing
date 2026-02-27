#include "valuemap.h"
#include "colors.h"

ValueMap::ValueMap(const int _data, const int _width, const int _height) {
    double* _data_ptr = new double[1]{(double) _data};
    initialise(_data_ptr, _width, _height);
}

ValueMap::ValueMap(const double _data, const int _width, const int _height) {
    double* _data_ptr = new double[1]{_data};
    initialise(_data_ptr, _width, _height);
}

ValueMap::ValueMap(const vec3& _data, const int _width, const int _height) {
    // TODO: I dont like using new here, maybe move to std::vector?
    double* _data_ptr = new double[3]{_data[0], _data[1], _data[2]};
    initialise(_data_ptr, _width, _height);
}

ValueMap::ValueMap(double* _data, const int _width, const int _height) {
    initialise(_data, _width, _height);
}

ValueMap::~ValueMap() {
    delete[] data;
}

void ValueMap::initialise(double* _data, const int _width, const int _height) {
    data = _data;
    width = _width;
    height = _height;
}

double ValueMap1D::get(const double u, const double v) const {
    if (isnan(u) || isnan(v)) {
        return 0;
    }
    int u_idx = static_cast<int>(static_cast<double>(width) * pos_fmod(u, 1.0));
    int v_idx = static_cast<int>(static_cast<double>(height) * pos_fmod(1 - v, 1.0));
    int index = (v_idx * width + u_idx);
    return data[index];
}

vec3 ValueMap3D::get(const double u, const double v) const {
    if (isnan(u) || isnan(v)) {
        return vec3(0, 0, 0);
    }
    int u_idx = static_cast<int>(static_cast<double>(width) * pos_fmod(u, 1.0));
    int v_idx = static_cast<int>(static_cast<double>(height) * pos_fmod(1 - v, 1.0));
    int start_index = 3 * (v_idx * width + u_idx);
    return vec3(data[start_index], data[start_index + 1], data[start_index + 2]);
}

template<typename ValueMapType>
ValueMapType* create_value_map(const std::string& file_name, const bool gamma_correct) {
    std::ifstream file(file_name, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Could not open file '" + file_name + "'");
    }

    file.seekg(0, std::ios::end);
    std::streamsize size = file.tellg();
    file.seekg(0, std::ios::beg);

    size_t N = static_cast<size_t>(size) / sizeof(double);

    double* data = new double[N];
    if (!file.read(reinterpret_cast<char*>(data), size)) {
        throw std::runtime_error("Error reading file '" + file_name + "'");
    }

    int width = static_cast<int>(data[0]);
    int height = static_cast<int>(data[1]);
    int dimension = static_cast<int>(data[2]);
    // TODO: should cast before, in case the map is very large.
    size_t num_values = static_cast<size_t>(width * height * dimension);
    if (num_values != N - 3) {
        throw std::runtime_error("File '" + file_name + "' does not follow the expected format and cannot be loaded");
    }
    double* data_array = new double[N - 3];
    for (size_t i = 3; i < N; ++i) {
        data_array[i - 3] = gamma_correct ? apply_gamma_correction(data[i]) : data[i];
    }

    delete[] data;
    return new ValueMapType(data_array, width, height);
}

ValueMap1D const* create_value_map_1D(const std::string& file_name) {
    // Gamma correction does not make any sense for ValueMap1Ds.
    return create_value_map<ValueMap1D>(file_name, false);
}

ValueMap3D const* create_value_map_3D(const std::string& file_name, const bool gamma_correct) {
    return create_value_map<ValueMap3D>(file_name, gamma_correct);
}
