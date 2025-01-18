#include "medium.h"


Medium::Medium(const double _scattering_coefficient, const vec3& _absorption_albedo) : 
scattering_coefficient(_scattering_coefficient), absorption_albedo(_absorption_albedo) {
    extinction_albedo = absorption_albedo + vec3(scattering_coefficient);
}

double Medium::sample_distance() const{
    return constants::max_ray_distance;
}

vec3 Medium::sample_direction(const vec3& incident_vector) const{
    return sample_spherical();
}

vec3 Medium::transmittance_albedo(const double distance) const{
    return exp_vector(-extinction_albedo * distance);
}

void Medium::Integrate(Object** objects, const int number_of_objects, Ray& incoming_ray, vec3& Lv, vec3& transmittance, vec3& weight, Ray& outgoing_ray) const{
    Hit hit;
    if(!find_closest_hit(hit, incoming_ray, objects, number_of_objects)){
        return;
    }

    transmittance = colors::WHITE;
}

vec3 Medium::sample(Object** objects, const int number_of_objects, const double distance) const{
    return colors::WHITE;
}


BeersLawMedium::BeersLawMedium(const double _scattering_coefficient, const vec3& _absorption_albedo) : Medium(0, _absorption_albedo){}

void BeersLawMedium::Integrate(Object** objects, const int number_of_objects, Ray& incoming_ray, vec3& Lv, vec3& transmittance, vec3& weight, Ray& outgoing_ray) const{
    Hit hit;
    if(!find_closest_hit(hit, incoming_ray, objects, number_of_objects)){
        return;
    }

    transmittance = transmittance_albedo((hit.intersection_point - incoming_ray.starting_position).length());
    weight = 1;
    outgoing_ray = incoming_ray;
}

vec3 BeersLawMedium::sample(Object** object, const int number_of_objects, const double distance) const{
    return transmittance_albedo(distance);
}


void SingleScatteringHomogenousMedium::Integrate(Object** objects, const int number_of_objects, Ray& incoming_ray, 
vec3& L, vec3& transmittance, vec3& weight, Ray& outgoing_ray) const{
    Hit hit;
    if(!find_closest_hit(hit, incoming_ray, objects, number_of_objects)){
        return;
    }

    transmittance = transmittance_albedo((hit.intersection_point - incoming_ray.starting_position).length());
    double rand = random_uniform(0, 1);
    double scatter_distance = -std::log(1 - rand * (1 - transmittance.mean())) / extinction_albedo.mean();
    vec3 scatter_point = incoming_ray.starting_position + scatter_distance * incoming_ray.direction_vector;

    L = vec3(0);

    vec3 sampled_direction;
    L += direct_lighting(incoming_ray.starting_position, objects, number_of_objects, sampled_direction) * 0.25 / M_PI;

    vec3 tr = transmittance_albedo(scatter_distance);
    L *= extinction_albedo * scattering_coefficient / extinction_albedo * tr;

    weight = (vec3(1,1,1) - transmittance) / (tr * extinction_albedo); 
    outgoing_ray = incoming_ray;
}


double ScatteringMediumHomogenous::sample_distance() const{
    return -std::log(1 -random_uniform(0, 1)) / extinction_albedo.mean();
}

vec3 ScatteringMediumHomogenous::sample(Object** objects, const int number_of_objects, const double distance) const{
    return transmittance_albedo(distance);
}


MediumStack::~MediumStack(){
    delete[] medium_stack;
}

Medium* MediumStack::get_medium() const{
    if(stack_size == 0){
        return nullptr;
    }
    
    return medium_stack[stack_size-1];
}

void MediumStack::add_medium(Medium* medium, const int id){
    // Call this when entering a new medium.
    //std::cout << "Add! Size (before addition is made): " << stack_size << ", id: " << id << "\n";
    if (stack_size == MAX_STACK_SIZE){
        throw std::invalid_argument("Cannot add another medium to stack, stack is full!");
    }
    bool found = false;
    for (int i = 0; i < stack_size; i++){
        if (medium_stack[i] -> id == id){
            //std::cout << "Error: 1. Found a medium with same ID!\n"; 
            found = true;
        }
    }
    
    medium -> id = id;
    medium_stack[stack_size] = medium;
    stack_size++;
}

void MediumStack::pop_medium(const int id){
    // Call this when exiting a medium.
    //std::cout << "Remove! Size (before addition is made): " << stack_size << ", id: " << id << "\n";
    for (int i = stack_size-1; i >= 0; i--){
        if (medium_stack[i] -> id == id){
            medium_stack[i] = nullptr;
            stack_size--;
            return;
        }
    }
    //std::cout << "Error: 2. Could not remove medium from stack because no matching medium was found!\n";
    //throw std::invalid_argument("Could not remove medium from stack because no matching medium was found!");
}
