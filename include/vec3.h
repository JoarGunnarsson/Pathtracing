#ifndef VEC3_H
#define VEC3_H

#include <cmath>
#include <iostream>

class vec3 {
  public:
    double e[3];

    vec3() : e{0, 0, 0} {}

    vec3(double e0, double e1, double e2) : e{e0, e1, e2} {}
    vec3(double e) : e{e, e, e} {}

    inline vec3 operator-() const {
        return vec3(-e[0], -e[1], -e[2]);
    }

    inline vec3 operator+(const vec3& v2) const {
        double e0 = e[0] + v2[0];
        double e1 = e[1] + v2[1];
        double e2 = e[2] + v2[2];
        return vec3(e0, e1, e2);
    }

    inline vec3& operator+=(const vec3& v2) {
        e[0] += v2[0];
        e[1] += v2[1];
        e[2] += v2[2];
        return *this;
    }

    inline vec3 operator-(const vec3& v2) const {
        double e0 = e[0] - v2[0];
        double e1 = e[1] - v2[1];
        double e2 = e[2] - v2[2];
        return vec3(e0, e1, e2);
    }

    inline vec3& operator-=(const vec3& v2) {
        e[0] -= v2[0];
        e[1] -= v2[1];
        e[2] -= v2[2];
        return *this;
    }

    inline vec3 operator*(const vec3& v2) const {
        double e0 = e[0] * v2[0];
        double e1 = e[1] * v2[1];
        double e2 = e[2] * v2[2];
        return vec3(e0, e1, e2);
    }

    inline vec3 operator*(const double value) const {
        double e0 = e[0] * value;
        double e1 = e[1] * value;
        double e2 = e[2] * value;
        return vec3(e0, e1, e2);
    }

    inline vec3& operator*=(const vec3& v2) {
        e[0] *= v2[0];
        e[1] *= v2[1];
        e[2] *= v2[2];
        return *this;
    }

    inline vec3& operator*=(const double value) {
        e[0] *= value;
        e[1] *= value;
        e[2] *= value;
        return *this;
    }

    inline vec3& operator/=(const double value) {
        e[0] /= value;
        e[1] /= value;
        e[2] /= value;
        return *this;
    }

    inline vec3 operator/(const double value) const {
        double e0 = e[0] / value;
        double e1 = e[1] / value;
        double e2 = e[2] / value;
        return vec3(e0, e1, e2);
    }

    inline vec3 operator/(const vec3& v) const {
        double e0 = e[0] / v[0];
        double e1 = e[1] / v[1];
        double e2 = e[2] / v[2];
        return vec3(e0, e1, e2);
    }

    inline double& operator[](int i) {
        return e[i];
    }

    inline double length() const {
        return std::sqrt(length_squared());
    }

    inline double length_squared() const {
        return e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
    }

    inline double max() const {
        return std::max(std::max(e[0], e[1]), e[2]);
    }

    inline double mean() const {
        return (e[0] + e[1] + e[2]) / 3.0;
    }

    operator double*() const {
        // Cast vec3 to double array, used for ValueMaps, which are responsible for deleting resource.
        // TODO: Migrate to std::vector instead
        double* result = new double[3];
        result[0] = e[0];
        result[1] = e[1];
        result[2] = e[2];
        return result;
    }
};

inline vec3 operator*(double value, const vec3& v) {
    double v0 = v[0] * value;
    double v1 = v[1] * value;
    double v2 = v[2] * value;
    return vec3(v0, v1, v2);
}

inline double dot_vectors(const vec3& v1, const vec3& v2) {
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

inline vec3 cross_vectors(const vec3& v1, const vec3& v2) {
    return vec3(v1[1] * v2[2] - v1[2] * v2[1], v1[2] * v2[0] - v1[0] * v2[2], v1[0] * v2[1] - v1[1] * v2[0]);
}

inline vec3 normalize_vector(const vec3& v) {
    return v / v.length();
}

inline vec3 exp_vector(const vec3& v) {
    return vec3(std::exp(v[0]), std::exp(v[1]), std::exp(v[2]));
}

inline vec3 abs(const vec3& v) {
    return vec3(std::abs(v[0]), std::abs(v[1]), std::abs(v[2]));
}

inline int argmax(const vec3& v) {
    int max_index = 0;
    if (v[1] > v[max_index]) {
        max_index = 1;
    }
    if (v[2] > v[max_index]) {
        max_index = 2;
    }
    return max_index;
}

inline vec3 permute(const vec3& v, const int idx1, const int idx2, const int idx3) {
    return vec3(v[idx1], v[idx2], v[idx3]);
}
void display_vector(const vec3& v);

#endif
