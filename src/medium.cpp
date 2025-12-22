#include "medium.h"

// TODO: Does media even need all existing objects? For heterogenous scattering maybe?
Medium::Medium(const vec3& _scattering_albedo, const vec3& _absorption_albedo, const vec3& _emission_coefficient,
               const double _refractive_index) :
    refractive_index(_refractive_index),
    absorption_albedo(_absorption_albedo),
    scattering_albedo(_scattering_albedo),
    emission_coefficient(_emission_coefficient) {
    extinction_albedo = absorption_albedo + scattering_albedo;
}

double Medium::sample_distance() const {
    return constants::max_ray_distance;
}

vec3 Medium::sample_direction(const vec3&) const {
    return sample_spherical();
}

double Medium::phase_function(const vec3&, const vec3&) const {
    return 1.0 / (4 * M_PI);
}

vec3 Medium::transmittance_albedo(const double distance) const {
    return exp_vector(-extinction_albedo * distance);
}

vec3 Medium::sample(Object**, const int, const double, const bool) const {
    return colors::WHITE;
}

vec3 Medium::sample_emission() const {
    return colors::BLACK;
}

BeersLawMedium::BeersLawMedium(const vec3&, const vec3& _absorption_albedo, const vec3& _emission_coefficient,
                               const double _refractive_index) :
    Medium(0, _absorption_albedo, _emission_coefficient, _refractive_index) {}

vec3 BeersLawMedium::sample(Object**, const int, const double distance, const bool) const {
    return transmittance_albedo(distance);
}

double HomogenousScatteringMedium::sample_distance() const {
    int channel = random_int(0, 3);
    if (extinction_albedo[channel] == 0) {
        return constants::max_ray_distance;
    }
    return -std::log(random_uniform(0, 1)) / extinction_albedo[channel];
}

vec3 HomogenousScatteringMedium::sample(Object**, const int, const double distance, const bool scatter) const {
    vec3 tr = transmittance_albedo(distance);
    vec3 density = scatter ? extinction_albedo * tr : tr;
    double pdf = 0;
    for (int i = 0; i < 3; i++) {
        pdf += density[i];
    }
    pdf *= 1.0 / 3.0;

    return scatter ? tr * scattering_albedo / pdf : tr / pdf;
}

vec3 HomogenousScatteringMedium::sample_emission() const {
    double pdf = 0;
    for (int i = 0; i < 3; i++) {
        pdf += extinction_albedo[i];
    }
    pdf *= 1.0 / 3.0;
    return emission_coefficient * absorption_albedo / pdf;
}
