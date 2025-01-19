#ifndef UTILS_H
#define UTILS_H

#include "vec3.h"
#include "constants.h"
#include <random>
#include <complex>


double random_uniform(const double low, const double high);
int random_int(const int low, const int high);
double random_normal();

enum reflection_type{
    DIFFUSE = 0,
    REFLECTED = 1,
    TRANSMITTED = 2
};

struct Hit{
    int intersected_object_index = -1;
    int primitive_ID = -1;
    double distance = constants::max_ray_distance;
    vec3 intersection_point;
    vec3 incident_vector;
    vec3 normal_vector;
};

struct Ray{
    vec3 starting_position;
    vec3 direction_vector;
    int type = DIFFUSE;
    double t_max = constants::max_ray_distance;
    double debug;
    int kx;
    int ky;
    int kz;
    vec3 d;
    double Sx;
    double Sy;
    double Sz;

    void prepare(){
        kz = argmax(abs(direction_vector));
        kx = kz + 1;
        if (kx == 3){
            kx = 0;
        }

        ky = kx + 1;
        if (ky == 3){
            ky = 0;
        }

        d = permute(direction_vector, kx, ky, kz);

        Sx = -d[0]/ d[2];
        Sy = -d[1] / d[2];
        Sz = 1.0 / d[2];
    }
};


double pos_fmod(const double a, const double b);
double clamp(const double value, const double min, const double max);
double sign(const double x);
bool solve_quadratic(const double b, const double c, double& distance);

vec3 sample_spherical();
vec3 sample_hemisphere(const vec3& normal);
void set_perpendicular_vectors(const vec3& z_hat, vec3& x_hat, vec3& y_hat);
vec3 sample_angled_hemisphere(const vec3& normal_vector, const double cos_max);
vec3 sample_cosine_hemisphere(const vec3& normal_vector);
vec3 reflect_vector(const vec3& direction_vector, const vec3& normal_vector);
vec3 refract_vector(const vec3& incident_vector, const vec3& normal_vector, const double eta);

double fresnel_dielectric(const double cos_incident, const double n1, const double n2);
double fresnel_conductor(double cos_theta_real, const double n1, const double k1, const double n2, const double k2);
double fresnel_multiplier(const double cos_incident, const double n1, const double k1, const double n2, const double k2, const bool is_dielectric);
#endif
