#ifndef MEDIUM_H
#define MEDIUM_H

#include <stdexcept>
#include "utils.h"
#include "objects.h"

class Object;

class Medium {
  public:
    double refractive_index;
    Medium(const vec3& _scattering_albedo, const vec3& _absorption_albedo, const vec3& _emission_coefficient,
           const double _refractive_index);
    virtual ~Medium() {}

    virtual double sample_distance() const;
    virtual vec3 sample_direction(const vec3& incident_vector) const;
    virtual double phase_function(const vec3& incident_vector, const vec3& outgoing_vector) const;
    virtual vec3 transmittance_albedo(const double distance) const;
    virtual vec3 sample(Object** objects, const int number_of_objects, const double distance, const bool scatter) const;
    virtual vec3 sample_emission() const;

  protected:
    vec3 absorption_albedo;
    vec3 scattering_albedo;
    vec3 extinction_albedo;
    vec3 emission_coefficient;
};

class BeersLawMedium : public Medium {
  public:
    BeersLawMedium(const vec3& scattering_albedo, const vec3& _absorption_albedo, const vec3& _emission_coefficient,
                   const double _refractive_index);
    virtual vec3 sample(Object** objects, const int number_of_objects, const double distance,
                        const bool scatter) const override;
};

class HomogenousScatteringMedium : public Medium {
  public:
    using Medium::Medium;
    virtual double sample_distance() const override;
    virtual vec3 sample(Object** objects, const int number_of_objects, const double distance,
                        const bool scatter) const override;
    virtual vec3 sample_emission() const override;
};

#endif
