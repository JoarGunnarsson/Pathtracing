#ifndef MEDIUM_H
#define MEDIUM_H

#include <stdexcept>
#include "utils.h"
#include "objects.h"

class Object;


class Medium{
    public:
        int id;
        Medium(const vec3& _scattering_albedo, const vec3& _absorption_albedo);

        virtual double sample_distance() const;
        virtual vec3 sample_direction(const vec3& incident_vector) const;
        virtual vec3 transmittance_albedo(const double distance) const;
        virtual void Integrate(Object** objects, const int number_of_objects, Ray& incoming_ray, vec3& Lv, vec3& transmittance, vec3& weight, Ray& outgoing_ray) const;
        virtual vec3 sample(Object** objects, const int number_of_objects, const double distance, const bool scatter) const;
        virtual vec3 sample_direct(const vec3& scattering_point, Object** objects, const int number_of_objects) const;

    protected:
        vec3 scattering_albedo;
        vec3 absorption_albedo;
        vec3 extinction_albedo;
};


class BeersLawMedium: public Medium{
    public:
        BeersLawMedium(const vec3& scattering_albedo, const vec3& _absorption_albedo);

        virtual void Integrate(Object** objects, const int number_of_objects, Ray& incoming_ray, vec3& Lv, vec3& transmittance, vec3& weight, Ray& outgoing_ray) const override;
        virtual vec3 sample(Object** objects, const int number_of_objects, const double distance, const bool scatter) const override;
};


class SingleScatteringHomogenousMedium: public Medium{
    public:
        using Medium::Medium;

        virtual void Integrate(Object** objects, const int number_of_objects, Ray& incoming_ray, vec3& L, vec3& transmittance, vec3& weight, Ray& outgoing_ray) const override;

};


class ScatteringMediumHomogenous : public Medium{
    public:
        using Medium::Medium;
        virtual double sample_distance() const override;
        virtual vec3 sample(Object** objects, const int number_of_objects, const double distance, const bool scatter) const override;
        virtual vec3 sample_direct(const vec3& scattering_point, Object** objects, const int number_of_objects) const override;
};


class MediumStack{
    public:
        MediumStack(){}
        ~MediumStack();

        Medium* get_medium() const;
        void add_medium(Medium* medium, const int id);
        void pop_medium(const int id);

    private:
        const int MAX_STACK_SIZE = 50;
        int stack_size = 0;
        Medium** medium_stack = new Medium*[MAX_STACK_SIZE];
};

#endif