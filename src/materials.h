#ifndef MATERIALS_H
#define MATERIALS_H
#include "vec3.h"
#include "colors.h"
#include "utils.h"
#include "constants.h"
#include "valuemap.h"
#include "medium.h"


class Object;


Hit find_closest_hit(const Ray& ray, Object** objects, const int size);


struct BrdfData{
    vec3 outgoing_vector;
    vec3 brdf_multiplier;
    int type = DIFFUSE;
};


struct MicrofacetSampleArgs{
    vec3 sampled_half_vector;
    vec3 normal_vector;
    vec3 incident_vector;
    vec3 intersection_point;
    float cosine_factor;
    float eta;
    float u;
    float v;
    float alpha;
    Object** objects;
    int number_of_objects;
    bool outside;
};


struct MicrofacetData{
    bool outside;
    float alpha;
    float eta;
    vec3 normal_into_interface;
    vec3 half_vector;
    float F_r;
};


struct MaterialData{
    ValueMap3D* albedo_map = nullptr;
    float refractive_index = 1;
    float extinction_coefficient = 0;
    float attenuation_coefficient = 0;
    vec3 absorption_albedo = colors::WHITE;
    ValueMap3D* emission_color_map = nullptr;
    ValueMap1D* light_intensity_map = nullptr;
    bool is_dielectric = true;
    ValueMap1D* roughness_map = nullptr; 
    ValueMap1D* percentage_diffuse_map = nullptr;
    bool is_light_source = false;
    double scattering_coefficient = 0;
};


class Material{
    public:
        ValueMap3D* albedo_map;
        float refractive_index;
        float attenuation_coefficient;
        vec3 absorption_albedo;
        ValueMap3D* emission_color_map;
        ValueMap1D* light_intensity_map;
        bool is_dielectric;
        float extinction_coefficient;
        bool is_light_source;
        ValueMap1D* roughness_map;
        ValueMap1D* percentage_diffuse_map;
        Medium* medium;
        Material(){}
        Material(MaterialData data){
            if (!data.albedo_map){
                data.albedo_map = new ValueMap3D(colors::WHITE);
            }

            if (!data.emission_color_map){
                data.emission_color_map = new ValueMap3D(colors::WHITE);
            }

            if (!data.light_intensity_map){
                data.light_intensity_map = new ValueMap1D(0);
            }

            if (!data.roughness_map){
                data.roughness_map = new ValueMap1D(0);
            }

            if (!data.percentage_diffuse_map){
                data.percentage_diffuse_map = new ValueMap1D(1);
            }

            albedo_map = data.albedo_map;
            refractive_index = data.refractive_index;
            attenuation_coefficient = data.attenuation_coefficient;
            absorption_albedo = data.absorption_albedo;
            emission_color_map = data.emission_color_map;
            light_intensity_map = data.light_intensity_map;
            is_dielectric = data.is_dielectric;
            if (is_dielectric){
                extinction_coefficient = 0;
            }
            else{
                extinction_coefficient = data.extinction_coefficient;
            }
            is_light_source = data.is_light_source;

            roughness_map = data.roughness_map;
            percentage_diffuse_map = data.percentage_diffuse_map;

            medium = new Medium(data.attenuation_coefficient, data.scattering_coefficient, absorption_albedo);
        }
    
        ~Material(){
            delete albedo_map;
            delete emission_color_map;
            delete light_intensity_map;
            delete roughness_map;
            delete percentage_diffuse_map;
            delete medium;
        }

    virtual vec3 eval(const Hit& hit, const float u, const float v){
        throw VirtualMethodNotAllowedException("This is a pure virtual method and should not be called.");
        vec3 vec;
        return vec;
    }

    virtual BrdfData sample(const Hit& hit, Object** scene_objects, const int number_of_objects, const float u, const float v){
        throw VirtualMethodNotAllowedException("This is a pure virtual method and should not be called.");
        BrdfData data;
        return data;
    }

    vec3 get_light_emittance(const float u, const float v){
        return emission_color_map -> get(u, v) * light_intensity_map -> get(u, v);
    }
};


class DiffuseMaterial : public Material{
    public:
        using Material::Material;

    vec3 eval(const Hit& hit, const float u, const float v) override{
        return albedo_map -> get(u, v) / M_PI;
    }

    BrdfData sample(const Hit& hit, Object** scene_objects, const int number_of_objects, const float u, const float v) override{
        vec3 adjusted_normal = dot_vectors(hit.incoming_vector, hit.normal_vector) < 0 ? hit.normal_vector : -hit.normal_vector;
        vec3 outgoing_vector = sample_cosine_hemisphere(adjusted_normal);
        vec3 brdf_multiplier = albedo_map -> get(u, v);
        BrdfData data;
        data.outgoing_vector = outgoing_vector;
        data.brdf_multiplier = brdf_multiplier;
        data.type = DIFFUSE;
        return data;
    }
};


class ReflectiveMaterial : public Material{
    public:
        using Material::Material;

    vec3 eval(const Hit& hit, const float u, const float v) override{
        return colors::BLACK;
    }

    BrdfData sample(const Hit& hit, Object** scene_objects, const int number_of_objects, const float u, const float v) override{
        vec3 adjusted_normal = dot_vectors(hit.incoming_vector, hit.normal_vector) < 0 ? hit.normal_vector : -hit.normal_vector;
        vec3 outgoing_vector = reflect_vector(hit.incoming_vector, adjusted_normal);
        BrdfData data;
        data.outgoing_vector = outgoing_vector;
        data.brdf_multiplier = is_dielectric ? colors::WHITE : albedo_map -> get(u, v);
        data.type = REFLECTED;
        return data;
    }
};


class TransparentMaterial : public Material{
    public:
        using Material::Material;

    vec3 eval(const Hit& hit, const float u, const float v) override{
        return colors::BLACK;
    }

    BrdfData sample(const Hit& hit, Object** scene_objects, const int number_of_objects, const float u, const float v) override{
        float incoming_dot_normal = dot_vectors(hit.incoming_vector, hit.normal_vector);
        vec3 normal_into_interface;
        bool inside = incoming_dot_normal > 0.0;
        // Look at the next medium, use that as refractive index!
        float n1;
        float k1;
        float n2;
        float k2;
        if (!inside){
            normal_into_interface = -hit.normal_vector;
            n1 = constants::air_refractive_index;
            k1 = 0;
            n2 = refractive_index;
            k2 = extinction_coefficient;
        }
        else{
            normal_into_interface = hit.normal_vector;
            n1 = refractive_index;
            k1 = extinction_coefficient;
            n2 = constants::air_refractive_index;
            k2 = 0;
        }

        vec3 transmitted_vector = refract_vector(hit.incoming_vector, normal_into_interface, n1 / n2);

        float F_r = 1;
        if (transmitted_vector.length_squared() != 0){
            float cos_incident = dot_vectors(hit.incoming_vector, normal_into_interface);
            F_r = fresnel_multiplier(cos_incident, n1, k1, n2, k2, is_dielectric);
        }

        float random_num = random_uniform(0, 1);
        bool is_reflected = random_num <= F_r;
        
        BrdfData data;
        if (is_reflected){
            data.type = REFLECTED;
            data.brdf_multiplier = is_dielectric ? colors::WHITE : albedo_map -> get(u, v);
            data.outgoing_vector = reflect_vector(hit.incoming_vector, normal_into_interface);
        }
        else{
            data.type = TRANSMITTED;
            /*
            Ray transmission_ray;
            transmission_ray.direction_vector = transmitted_vector;
            transmission_ray.starting_position = hit.intersection_point;
            Hit transmission_hit = find_closest_hit(transmission_ray, scene_objects, number_of_objects);
              
            vec3 attenuation_color;
            float distance = transmission_hit.distance;
            if (distance > 0 && !inside){
                vec3 log_attenuation = absorption_albedo * attenuation_coefficient * (-distance);
                attenuation_color = exp_vector(log_attenuation);
            }
            else{
                attenuation_color = colors::WHITE;
            }
            */
            data.brdf_multiplier = colors::WHITE;//attenuation_color;
            data.outgoing_vector = transmitted_vector;
        }

        return data;
    }
};


class MicrofacetMaterial : public Material{
    public:
        using Material::Material;

    float chi(const float x){
        return x > 0 ? 1 : 0;
    }

    float G1(const vec3& half_vector, const vec3& normal_vector, const vec3& v, const float alpha){
        float cos_theta = dot_vectors(half_vector, v);
        float tan_theta = sqrt((1.0 - cos_theta * cos_theta) / (cos_theta * cos_theta));
        float a = 1.0 / (alpha * tan_theta);

        float multiplicator_approximation = a < 1.6 ? (3.535 * a  + 2.181 * a * a) / (1.0 + 2.276 * a + 2.577 * a * a) : 1.0;
        
        return chi(cos_theta / dot_vectors(v, normal_vector)) * multiplicator_approximation;
    }

    float G(const vec3& half_vector, const vec3& normal_vector, const vec3& incident_vector, const vec3& outgoing_vector, const float alpha){
        return G1(half_vector, normal_vector, incident_vector, alpha) * G1(half_vector, normal_vector, outgoing_vector, alpha);
    }

    MicrofacetData prepare_microfacet_data(const Hit& hit, const float u, const float v){
        float n1;
        float k1;
        float n2;
        float k2;
        vec3 normal_into_interface;
        float incoming_dot_normal = dot_vectors(hit.incoming_vector, hit.normal_vector);
        bool outside = incoming_dot_normal <= 0.0;
        if (outside){
            normal_into_interface = -hit.normal_vector;
            n1 = constants::air_refractive_index;
            k1 = 0;
            n2 = refractive_index;
            k2 = extinction_coefficient;
        }
        else{
            normal_into_interface = hit.normal_vector;
            n1 = refractive_index;
            k1 = extinction_coefficient;
            n2 = constants::air_refractive_index;
            k2 = 0;
        }

        float alpha = std::max(roughness_map -> get(u, v), constants::EPSILON);

        vec3 sampled_half_vector = specular_sample(-normal_into_interface, alpha);
                
        float i_dot_h = dot_vectors(hit.incoming_vector, sampled_half_vector);

        float F_r = fresnel_multiplier(-i_dot_h, n1, k1, n2, k2, is_dielectric);

        MicrofacetData data;
        data.outside = outside;
        data.alpha = alpha;
        data.eta = n1 / n2;
        data.normal_into_interface = normal_into_interface;
        data.half_vector = sampled_half_vector;
        data.F_r = F_r;
        return data;
    }

    vec3 eval(const Hit& hit, const float u, const float v) override{
        MicrofacetData data = prepare_microfacet_data(hit, u, v);
        
        float random_num = random_uniform(0, 1);
        bool reflect_specular = random_num < data.F_r;
        if (reflect_specular){
            return colors::BLACK;
        }
        float transmit = random_uniform(0, 1) > percentage_diffuse_map -> get(u, v);

        if (transmit){
            return colors::BLACK;
        }

        return albedo_map -> get(u, v) / M_PI;
    }

    vec3 specular_sample(const vec3& normal_vector, const float alpha){
        float r1 = random_uniform(0, 1);
        float r2 = random_uniform(0, 1);
        float phi = 2 * M_PI * r2;
        float tan_theta2 = - alpha * alpha * std::log(1 - r1);
        float cos_theta2 = 1.0 / (1.0 + tan_theta2);
        
        float cos_theta = sqrt(cos_theta2);
        float sin_theta = sqrt(1 - cos_theta2);

        vec3 x_hat; 
        vec3 y_hat;
        set_perpendicular_vectors(normal_vector, x_hat, y_hat);
        return x_hat * sin_theta * cos(phi) + y_hat * sin_theta * sin(phi) + normal_vector * cos_theta;
    }

    BrdfData sample_diffuse(const MicrofacetSampleArgs& args){
        BrdfData data;
        data.brdf_multiplier = albedo_map -> get(args.u, args.v);
        data.outgoing_vector = sample_cosine_hemisphere(args.normal_vector);
        data.type = DIFFUSE;
        return data;
    }

    BrdfData sample_reflection(const MicrofacetSampleArgs& args){
        BrdfData data;
        vec3 reflection_color = is_dielectric ? colors::WHITE : albedo_map -> get(args.u, args.v);
        data.outgoing_vector = reflect_vector(-args.incident_vector, args.sampled_half_vector);
        data.brdf_multiplier = reflection_color * G(args.sampled_half_vector, args.normal_vector, args.incident_vector, data.outgoing_vector, args.alpha) * args.cosine_factor;
        data.type = REFLECTED;
        return data;
    }

    vec3 compute_attenuated_color(const MicrofacetSampleArgs& args, const vec3& outgoing_vector){
        Ray transmission_ray;
        transmission_ray.direction_vector = outgoing_vector;
        transmission_ray.starting_position = args.intersection_point;
        Hit transmission_hit = find_closest_hit(transmission_ray, args.objects, args.number_of_objects);
        vec3 attenuation_color;
        float distance = transmission_hit.distance;
        if (distance > 0 && args.outside){
            vec3 log_attenuation = absorption_albedo * attenuation_coefficient * (-distance);
            attenuation_color = exp_vector(log_attenuation);
        }
        else{
            attenuation_color = colors::WHITE;
        }

        return attenuation_color;
    }

    BrdfData sample_transmission(const MicrofacetSampleArgs& args){
        vec3 refracted_vector = refract_vector(-args.incident_vector, -args.sampled_half_vector, args.eta);

        if (refracted_vector.length() == 0){
            return sample_reflection(args);
        }

        vec3 attenuated_color = colors::WHITE;//compute_attenuated_color(args, refracted_vector);

        BrdfData data;
        data.outgoing_vector = refracted_vector;

        data.brdf_multiplier = attenuated_color * G(args.sampled_half_vector,args.normal_vector, args.incident_vector, data.outgoing_vector, args.alpha) * args.cosine_factor;
        data.type = TRANSMITTED;
        return data;
    }
    
    
    BrdfData sample(const Hit& hit, Object** objects, const int number_of_objects, const float u, const float v) override{
        MicrofacetData data = prepare_microfacet_data(hit, u, v);
                
        float i_dot_h = dot_vectors(hit.incoming_vector, data.half_vector);
        float i_dot_n = dot_vectors(hit.incoming_vector, data.normal_into_interface);
        float n_dot_h = dot_vectors(data.half_vector, data.normal_into_interface);

        float cosine_factor = std::abs(i_dot_h / (i_dot_n * n_dot_h));

        float random_num = random_uniform(0, 1);
        bool reflect_specular = random_num < data.F_r;
        
        MicrofacetSampleArgs args;
        args.sampled_half_vector = data.half_vector;
        args.normal_vector = -data.normal_into_interface;
        args.incident_vector = -hit.incoming_vector;
        args.cosine_factor = cosine_factor;
        args.u = u;
        args.v = v;
        args.alpha = data.alpha;
        if (reflect_specular){
            return sample_reflection(args);
        }
        else{
            args.intersection_point = hit.intersection_point;
            args.eta = data.eta;
            args.objects = objects;
            args.number_of_objects = number_of_objects;
            args.outside = data.outside;
            float random_num2 = random_uniform(0, 1);
            bool diffuse = random_num2 < percentage_diffuse_map -> get(u, v);
            return diffuse ? sample_diffuse(args) : sample_transmission(args);
        }
    }
};


class MaterialManager{
    const int MAX_MATERIALS = 100;
    Material** material_array = new Material*[MAX_MATERIALS];
    int current_idx = 0;

    public:
        MaterialManager(){}
        ~MaterialManager(){
            for (int i = 0; i < current_idx; i++){
                delete material_array[i];
            }
        }

        void add_material(Material* material){
            material_array[current_idx++] = material;
        };
};

#endif