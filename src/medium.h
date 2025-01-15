#ifndef MEDIUM_H
#define MEDIUM_H
#include <stdexcept>
#include "utils.h"

class Medium{
    public:
        const double refractive_index = 1;
        int id;
        double attenuation_coefficient;
        double scattering_coefficient;
        vec3 absorption_albedo;

        Medium(const double _attenuation_coefficient, const double _scattering_coefficient, const vec3& _absorption_albedo){
            attenuation_coefficient = _attenuation_coefficient;
            scattering_coefficient = _scattering_coefficient;
            absorption_albedo = _absorption_albedo;
        }

    double sample_distance(){
        if (scattering_coefficient == 0){
            return -1;
        }
        double rand = random_uniform(0, 1);
        return - std::log(1 - rand) / (scattering_coefficient);
    }

    vec3 sample_direction(const vec3& incident_vector){
        return sample_spherical();
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
                throw std::invalid_argument("Could get medium since stack is empty.");
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