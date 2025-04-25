#ifndef UTILS_H
#define UTILS_H
#include "vec3.h"
#include "constants.h"
#include <random>
#include <complex>

std::random_device rand_dev_normal;
std::minstd_rand normal_generator(rand_dev_normal());
std::normal_distribution<double> normal_distribution(0, 1);

std::random_device rand_dev_uniform;
std::minstd_rand uniform_generator(rand_dev_uniform());
std::uniform_real_distribution<double> uniform_dist(0, 1);

double random_uniform(const double low, const double high){
    return (high - low) * uniform_dist(uniform_generator) + low;
}

int random_int(const int low, const int high){
    // Returns a random int between low (inclusive) and high (exclusive). 
    return (int) random_uniform(low, high);
}

double random_normal(){
    return normal_distribution(normal_generator);
}

double pos_fmod(const double a, const double b){
    return fmod((fmod(a, b) + b), b);
}

double clamp(const double value, const double min, const double max){
    return std::max(std::min(value, max), min);
}

double sign(const double x){
    return x > 0 ? 1 : -1;
}

bool solve_quadratic(const double b, const double c, double& distance){
    double discriminant = b*b - 4.0 * c;
    if (discriminant <= 0){
        return false;
    }
    double root_discriminant = sqrt(discriminant);

    double minimum_solution = - 1.0 / 2.0 * (b + root_discriminant);
    if (minimum_solution > constants::EPSILON){
        distance = minimum_solution;
        return true;
    }

    double maximum_solution = - 1.0 / 2.0 * (b - root_discriminant);
    if (maximum_solution > constants::EPSILON){
        distance = maximum_solution;
        return true;
    }
    return false;

}

vec3 sample_spherical(){
    double r1 = random_normal();
    double r2 = random_normal();
    double r3 = random_normal();
    vec3 sample = vec3(r1, r2, r3);
    sample = normalize_vector(sample);
    return sample;
}

vec3 sample_hemisphere(const vec3& normal){
    vec3 sample = sample_spherical();
    if (dot_vectors(normal, sample) < 0){
        return -sample;
    }
    return sample;
}

void set_perpendicular_vectors(const vec3& z_hat, vec3& x_hat, vec3& y_hat){
    vec3 non_parallel_vector = vec3(1.0, 0.0, 0.0);
    if (std::abs(dot_vectors(non_parallel_vector, z_hat)) == 1.0){
        non_parallel_vector = vec3(0.0, 1.0, 0.0);
    }
    
    x_hat = cross_vectors(z_hat, non_parallel_vector);
    x_hat = normalize_vector(x_hat);
    y_hat = cross_vectors(z_hat, x_hat);
    y_hat = normalize_vector(y_hat);
}

vec3 sample_angled_hemisphere(const vec3& normal_vector, const double cos_max){
    vec3 x_hat;
    vec3 y_hat;
    set_perpendicular_vectors(normal_vector, x_hat, y_hat);
    double phi = random_uniform(0, 2.0 * M_PI);
    double cos_theta = random_uniform(cos_max, 1);
    double sin_theta = sqrt(1 - (cos_theta * cos_theta));
    double x = sin_theta * cos(phi);
    double y = sin_theta * sin(phi); 
    double z = cos_theta;
    return x_hat * x + y_hat * y + normal_vector * z;
}

vec3 sample_cosine_hemisphere(const vec3& normal_vector){
    vec3 x_hat;
    vec3 y_hat;
    set_perpendicular_vectors(normal_vector, x_hat, y_hat);

    double theta = ((double) rand() / (RAND_MAX)) * 2 * M_PI;
    double radius = sqrt((double) rand() / (RAND_MAX));
    double x = cos(theta) * radius;
    double y = sin(theta) * radius;
    double z = sqrt(1 - x * x - y * y);
    return x_hat * x + y_hat * y + normal_vector * z;
}

vec3 reflect_vector(const vec3& direction_vector, const vec3& normal_vector){
    return direction_vector - normal_vector * 2.0 * dot_vectors(normal_vector, direction_vector);
}

vec3 refract_vector(const vec3& incident_vector, const vec3& normal_vector, const double eta){
    double cos_incident = dot_vectors(normal_vector, incident_vector);
    double length_in_normal_direction_squared = 1 - eta * eta * (1 - cos_incident * cos_incident);
    if (length_in_normal_direction_squared < 0){
        return vec3(0,0,0);
    }
    vec3 perpendicular_vectors = incident_vector - normal_vector * cos_incident;
    return normal_vector * sqrt(length_in_normal_direction_squared) + perpendicular_vectors * eta;
}

double fresnel_dielectric(const double cos_incident, const double n1, const double n2){
    // Can result in nans - seems to be caused by rays entering through the surface of concave models.
    double sin_incident = sqrt(1 - cos_incident * cos_incident);
    double cos_transmitted = sqrt(1 - pow(n1 / n2 * sin_incident, 2));
    if (isnan(cos_transmitted)){
        return 0.0;
    }
    double n1_cos_incident = n1 * cos_incident;
    double n2_cos_transmitted = n2 * cos_transmitted;
    double n1_cos_transmitted = n1 * cos_transmitted;
    double n2_cos_incident = n2 * cos_incident;
    double R_s = pow((n1_cos_incident - n2_cos_transmitted) / (n1_cos_incident + n2_cos_transmitted), 2);
    double R_p = pow((n1_cos_transmitted - n2_cos_incident) / (n1_cos_transmitted + n2_cos_incident), 2);
    return 0.5 * (R_s + R_p);
}

double fresnel_conductor(double cos_theta_real, const double n1, const double k1, const double n2, const double k2){
    std::complex<double> cos_theta(cos_theta_real, 0);
    double eta;
    double k;
    std::complex<double> one(1,0);
    if (k1 == 0){
        eta = n2 / n1;
        k = k2 / n1;
    }
    else{
        eta = n1 / n2;
        k = k1 / n2;
        std::complex<double> sin_theta = sqrt(one - cos_theta * cos_theta);
        std::complex<double> n_tilde(n1, k1);
        std::complex<double> sin_theta_T = n_tilde / n2 * sin_theta;
        std::complex<double> sin_theta_T2 = sin_theta_T * sin_theta_T;
        std::complex<double> cos_theta = sqrt(one - sin_theta_T2);
    }

    std::complex<double> cos_theta2 = cos_theta * cos_theta;
    std::complex<double> sin_theta2 = one - cos_theta2;
    std::complex<double> f0 = sqrt(pow(eta * eta - k * k - sin_theta2, 2) + 4.0 * eta * eta * k * k);
    std::complex<double> a2b2 = f0;
    std::complex<double> a = sqrt(0.5 * f0 + eta * eta - k * k - sin_theta2);
    std::complex<double> f1 = a2b2 + cos_theta2;
    std::complex<double> f2 = 2.0 * a * cos_theta;
    std::complex<double> f3 = cos_theta2 * a2b2 + sin_theta2 * sin_theta2;
    std::complex<double> f4 = f2 * sin_theta2;

    std::complex<double> R_s = (f1 - f2) / (f1 + f2);
    std::complex<double> R_p = R_s * (f3 - f4) / (f3 + f4); 
    return 0.5 * real(R_p + R_s);
}

double fresnel_multiplier(const double cos_incident, const double n1, const double k1, const double n2, const double k2, const bool is_dielectric){
    double R_0_sqrt = (n1 - n2) / (n1 + n2);
    double R_0 = R_0_sqrt * R_0_sqrt;
    double x = 1 - cos_incident;
    double fresnel = R_0 + (1 - R_0) * (x * x * x * x * x);
    return fresnel;

    if (is_dielectric || (k1==0 && k2==0)){
        return fresnel_dielectric(cos_incident, n1, n2);
    }

    return fresnel_conductor(cos_incident, n1, k1, n2, k2);
}
#endif
