#include "medium.h"


int sample_random_light_2(Object** objects, const int number_of_objects, int& number_of_light_sources){
    int light_source_idx_array[number_of_objects];

    number_of_light_sources = 0;
    
    for (int i = 0; i < number_of_objects; i++){
        if (objects[i] -> is_light_source()){
            light_source_idx_array[number_of_light_sources] = i;
            number_of_light_sources++;
        }
    }
    if (number_of_light_sources == 0){
        return -1;
    }
    
    int random_index = random_int(0, number_of_light_sources);
    int light_index = light_source_idx_array[random_index];
    return light_index;
}


vec3 direct_lighting_2(const vec3& point, Object** objects, const int number_of_objects){
    int number_of_light_sources;
    int light_index = sample_random_light_2(objects, number_of_objects, number_of_light_sources);
    if (light_index == -1){
        return colors::BLACK;
    }

    double inverse_PDF;
    vec3 random_point = objects[light_index] -> random_light_point(point, inverse_PDF);

    vec3 vector_towards_light = random_point - point;
    double distance_to_light = vector_towards_light.length();
    vector_towards_light = normalize_vector(vector_towards_light);
    
    Ray light_ray;
    light_ray.starting_position = point;
    light_ray.direction_vector =  vector_towards_light;
    light_ray.prepare();

    // This does not need to be a full find_closest hit, we can return early. Should allow parameter ignore_index.
    Hit light_hit = find_closest_hit(light_ray, objects, number_of_objects);
    bool in_shadow = false;//light_hit.intersected_object_index != light_index;
    bool same_distance = true;//std::abs(distance_to_light - light_hit.distance) <= constants::EPSILON;

    if ( in_shadow || !same_distance){
        return colors::BLACK;
    }

    vec3 light_emitance = objects[light_index] -> get_light_emittance(light_hit);

    display_vector(light_emitance * inverse_PDF * (double) number_of_light_sources);
    return light_emitance * inverse_PDF * (double) number_of_light_sources;
}


Medium::Medium(const double _scattering_coefficient, const vec3& _absorption_albedo) : 
scattering_coefficient(_scattering_coefficient), absorption_albedo(_absorption_albedo) {
    extinction_albedo = absorption_albedo + vec3(scattering_coefficient);
}

double Medium::sample_distance() const{
    return -1;
}

vec3 Medium::sample_direction(const vec3& incident_vector) const{
    return sample_spherical();
}

void Medium::Integrate(Object** objects, const int number_of_objects, const Ray& incoming_ray, vec3& Lv, vec3& transmittance, vec3& weight, Ray& outgoing_ray){
    Hit hit = find_closest_hit(incoming_ray, objects, number_of_objects);
    if(hit.distance < constants::EPSILON){
        return;
    }

    transmittance = colors::WHITE;
}


vec3 BeersLawMedium::transmittance_color(const double distance) const{
    return exp_vector(-absorption_albedo * distance);
}

void BeersLawMedium::Integrate(Object** objects, const int number_of_objects, const Ray& incoming_ray, vec3& Lv, vec3& transmittance, vec3& weight, Ray& outgoing_ray){
    Hit hit = find_closest_hit(incoming_ray, objects, number_of_objects);
    if(hit.distance < constants::EPSILON){
        return;
    }

    transmittance = transmittance_color((hit.intersection_point - incoming_ray.starting_position).length());
    weight = 1;
    outgoing_ray = incoming_ray;
}


vec3 SingleScatteringHomogenousMedium::transmittance_color(const double distance) const{
    return exp_vector(-extinction_albedo * distance);
}

void SingleScatteringHomogenousMedium::Integrate(Object** objects, const int number_of_objects, const Ray& incoming_ray, vec3& L, vec3& transmittance, vec3& weight, Ray& outgoing_ray){
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

    L += direct_lighting_2(incoming_ray.starting_position, objects, number_of_objects) * 0.25 / M_PI;

    vec3 tr = transmittance_color(scatter_distance);
    L *= extinction_albedo * scattering_coefficient / extinction_albedo * tr;

    weight = (vec3(1,1,1) - transmittance) / (tr * extinction_albedo); 
    outgoing_ray = incoming_ray;
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
            std::cout << "Error: 1. Found a medium with same ID!\n"; 
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
