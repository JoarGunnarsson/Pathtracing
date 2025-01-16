#ifndef VALUEMAP_H
#define VALUEMAP_H

#include "vec3.h"
#include "utils.h"
#include <cstdio>


class ValueMap{
    public:
        double* data;
        int width;
        double u_max;
        int height;
        double v_max;

        ValueMap(){}
        ValueMap(const int _data, const int _width=1, const int _height=1, const double _u_max=1, const double _v_max=1);
        ValueMap(const double _data, const int _width=1, const int _height=1, const double _uMax=1, const double _v_max=1);
        ValueMap(double* _data, const int _width=1, const int _height=1, const double _u_max=1, const double _v_max=1);
        ~ValueMap();
        
    private: 
        void initialise(double* _data, const int _width, const int _height, const double _u_max, const double _v_max);
};


class ValueMap1D : public ValueMap{
    public:
        using ValueMap::ValueMap;

        double get(const double u, const double v) const;
};


class ValueMap3D : public ValueMap{
    public:
        using ValueMap::ValueMap;

        vec3 get(const double u, const double v) const;
};


ValueMap1D* create_value_map_1D(const char* file_name, double u_max = 1, double v_max = 1);
ValueMap3D* create_value_map_3D(const char* file_name, double u_max = 1, double v_max = 1);
#endif
