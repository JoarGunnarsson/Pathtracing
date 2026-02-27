#ifndef VALUEMAP_H
#define VALUEMAP_H

#include <cstdio>
#include <fstream>

#include "vec3.h"
#include "utils.h"

class ValueMap {
  public:
    ValueMap() {}
    ValueMap(const int _data, const int _width = 1, const int _height = 1);
    ValueMap(const double _data, const int _width = 1, const int _height = 1);
    ValueMap(const vec3& _data, const int _width = 1, const int _height = 1);
    ValueMap(double* _data, const int _width = 1, const int _height = 1);
    ~ValueMap();

  protected:
    double* data;
    int width;
    int height;
    void initialise(double* _data, const int _width, const int _height);
};

class ValueMap1D : public ValueMap {
  public:
    using ValueMap::ValueMap;

    double get(const double u, const double v) const;
};

class ValueMap3D : public ValueMap {
  public:
    using ValueMap::ValueMap;

    vec3 get(const double u, const double v) const;
};

ValueMap1D const* create_value_map_1D(const std::string& file_name);
ValueMap3D const* create_value_map_3D(const std::string& file_name, const bool gamma_correct);
#endif
