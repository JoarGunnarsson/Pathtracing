#ifndef MEDIUM_H
#define MEDIUM_H
#include <stdexcept>
#include "utils.h"
#include "objects.h"

class Object;


Hit find_closest_hit(const Ray& ray, Object** objects, const int size);


class Medium{
    public:
        const double refractive_index = 1;
        int id;
        double attenuation_coefficient;
        double scattering_coefficient;
        vec3 absorption_albedo;
        vec3 extinction_albedo;

        Medium(const double _attenuation_coefficient, const double _scattering_coefficient, const vec3& _absorption_albedo){
            attenuation_coefficient = _attenuation_coefficient;
            scattering_coefficient = _scattering_coefficient;
            absorption_albedo = _absorption_albedo;
            extinction_albedo = attenuation_coefficient * absorption_albedo + vec3(scattering_coefficient);
        }

    virtual double sample_distance(){
        if (scattering_coefficient == 0){
            return -1;
        }
        double rand = random_uniform(0, 1);
        return - std::log(1 - rand) / (scattering_coefficient);
    }

    virtual vec3 sample_direction(const vec3& incident_vector){
        return sample_spherical();
    }

    virtual void Integrate(Object** objects, const int number_of_objects, const Ray& incoming_ray, vec3& Lv, vec3& transmittance, vec3& weight, Ray& outgoing_ray){
        Hit hit = find_closest_hit(incoming_ray, objects, number_of_objects);
        if(hit.distance < constants::EPSILON){
            return;
        }

        transmittance = exp_vector(-absorption_albedo * attenuation_coefficient * (hit.intersection_point - incoming_ray.starting_position).length());
    }

};


class BeersLawMedium: public Medium{
    public:
        using Medium::Medium;

    vec3 transmittance_color(const double distance){
        return exp_vector(-absorption_albedo * attenuation_coefficient * distance);
    }

    virtual void Integrate(Object** objects, const int number_of_objects, const Ray& incoming_ray, vec3& Lv, vec3& transmittance, vec3& weight, Ray& outgoing_ray) override {
        Hit hit = find_closest_hit(incoming_ray, objects, number_of_objects);
        if(hit.distance < constants::EPSILON){
            return;
        }

        transmittance = transmittance_color((hit.intersection_point - incoming_ray.starting_position).length());
        weight = 1;
        outgoing_ray = incoming_ray;
    }
        
};



class SingleScatteringHomogenousMedium: public Medium{
    public:
        using Medium::Medium;

    vec3 transmittance_color(const double distance){
        return exp_vector(-extinction_albedo * distance);
    }

    virtual void Integrate(Object** objects, const int number_of_objects, const Ray& incoming_ray, vec3& L, vec3& transmittance, vec3& weight, Ray& outgoing_ray) override {
        Hit hit = find_closest_hit(incoming_ray, objects, number_of_objects);
        if(hit.distance < constants::EPSILON){
            return;
        }

        transmittance = transmittance_color((hit.intersection_point - incoming_ray.starting_position).length());

        double rand = random_uniform(0, 1);
        double scatter_distance = -std::log(1 - rand * (1 - transmittance.mean())) / extinction_albedo.mean();

        vec3 scatter_point = incoming_ray.starting_position + scatter_distance * incoming_ray.direction_vector;

        // Initialize IsotropicBSDF

        L = vec3(0);


        //L += direct_lighting_medium(incoming_ray.starting_position, objects, number_of_objects);

        vec3 tr = transmittance_color(scatter_distance);
        L *= extinction_albedo * scattering_coefficient / extinction_albedo * tr;

        weight = (vec3(1,1,1) - transmittance) / (tr * extinction_albedo); 
        //outgoing_ray = incoming_ray;
    }

};


class MediumStack{
    public:
        const int MAX_STACK_SIZE = 50;
        int stack_size = 0;
        Medium** medium_stack = new Medium*[MAX_STACK_SIZE];

        MediumStack(){}
        ~MediumStack(){
            delete[] medium_stack;
        }

        Medium* get_medium(){
            if(stack_size == 0){
                return nullptr;
            }
            
            return medium_stack[stack_size-1];
        }

        void add_medium(Medium* medium, const int id){
            // Call this when entering a new medium.
            //std::cout << "Add! Size (before addition is made): " << stack_size << ", id: " << id << "\n";
            if (stack_size == MAX_STACK_SIZE){
                throw std::invalid_argument("Cannot add another medium to stack, stack is full!");
            }
            medium -> id = id;
            medium_stack[stack_size] = medium;
            stack_size++;
        }

        void pop_medium(const int id){
            // Call this when exiting a medium.
            //std::cout << "Remove! Size (before addition is made): " << stack_size << ", id: " << id << "\n";
            for (int i = stack_size-1; i >= 0; i--){
                if (medium_stack[i] -> id == id){
                    medium_stack[i] = nullptr;
                    stack_size--;
                    return;
                }
            }

            //throw std::invalid_argument("Could not remove medium from stack because no matching medium was found!");
        }
};

#endif