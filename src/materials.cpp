#include "materials.h"


Material::Material(MaterialData data){
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

    medium = data.medium;
}

Material::~Material(){
    delete albedo_map;
    delete emission_color_map;
    delete light_intensity_map;
    delete roughness_map;
    delete percentage_diffuse_map;
    delete medium;
}

bool Material::allow_direct_light() const { return false; }
bool Material::compute_direct_light() const { return false; }
vec3 Material::eval(const Hit& hit, const double u, const double v) const{ return vec3(); }
BrdfData Material::sample(const Hit& hit, const double u, const double v) const{ return BrdfData(); }
double Material::brdf_pdf(const vec3& outgoing_vector, const vec3& incident_vector, const vec3& normal_vector, const double u, const double v) const{ return 0; }

vec3 Material::get_light_emittance(const double u, const double v) const{
    return emission_color_map -> get(u, v) * light_intensity_map -> get(u, v);
}


bool DiffuseMaterial::compute_direct_light() const{ 
    return true; 
}

vec3 DiffuseMaterial::eval(const Hit& hit, const double u, const double v) const{
    return albedo_map -> get(u, v) / M_PI;
}

BrdfData DiffuseMaterial::sample(const Hit& hit, const double u, const double v) const{
    vec3 adjusted_normal = dot_vectors(hit.incident_vector, hit.normal_vector) < 0 ? hit.normal_vector : -hit.normal_vector;
    vec3 outgoing_vector = sample_cosine_hemisphere(adjusted_normal);
    BrdfData data;
    data.outgoing_vector = outgoing_vector;
    data.brdf_over_pdf = albedo_map -> get(u, v);
    data.pdf = brdf_pdf(outgoing_vector, hit.incident_vector, hit.normal_vector, u, v);
    data.type = DIFFUSE;
    return data;
}

double DiffuseMaterial::brdf_pdf(const vec3& outgoing_vector, const vec3& incident_vector, const vec3& normal_vector, const double u, const double v) const{ 
    return std::abs(dot_vectors(outgoing_vector, normal_vector)) / M_PI; 
}


vec3 ReflectiveMaterial::eval(const Hit& hit, const double u, const double v) const{
    return colors::BLACK;
}

BrdfData ReflectiveMaterial::sample(const Hit& hit, const double u, const double v) const{
    vec3 adjusted_normal = dot_vectors(hit.incident_vector, hit.normal_vector) < 0 ? hit.normal_vector : -hit.normal_vector;
    vec3 outgoing_vector = reflect_vector(hit.incident_vector, adjusted_normal);
    BrdfData data;
    data.outgoing_vector = outgoing_vector;
    data.brdf_over_pdf = is_dielectric ? colors::WHITE : albedo_map -> get(u, v);
    data.type = REFLECTED;
    return data;
}


double ReflectiveMaterial::brdf_pdf(const vec3& outgoing_vector, const vec3& incident_vector, const vec3& normal_vector, const double u, const double v) const{ 
    return 0;
}


bool TransparentMaterial::allow_direct_light() const { return refractive_index == 1; }
vec3 TransparentMaterial::eval(const Hit& hit, const double u, const double v) const{
    return colors::BLACK;
}

BrdfData TransparentMaterial::sample(const Hit& hit, const double u, const double v) const{
    double incident_dot_normal = dot_vectors(hit.incident_vector, hit.normal_vector);
    vec3 normal_into_interface;
    bool outside = incident_dot_normal < 0.0;
    // TODO: Look at the next medium, use that as refractive index!
    double n1;
    double k1;
    double n2;
    double k2;
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

    vec3 transmitted_vector = refract_vector(hit.incident_vector, normal_into_interface, n1 / n2);

    double F_r = 1;
    if (transmitted_vector.length_squared() != 0){
        double cos_incident = dot_vectors(hit.incident_vector, normal_into_interface);
        F_r = fresnel_multiplier(cos_incident, n1, k1, n2, k2, is_dielectric);
    }

    double random_num = random_uniform(0, 1);
    bool is_reflected = random_num <= F_r;
    
    BrdfData data;
    if (is_reflected){
        data.type = REFLECTED;
        data.outgoing_vector = reflect_vector(hit.incident_vector, normal_into_interface);
        data.brdf_over_pdf = is_dielectric ? colors::WHITE : albedo_map -> get(u, v);
    }
    else{
        data.type = TRANSMITTED;
        data.outgoing_vector = transmitted_vector;
        data.brdf_over_pdf = colors::WHITE;
    }

    return data;
}

double TransparentMaterial::brdf_pdf(const vec3& outgoing_vector, const vec3& incident_vector, const vec3& normal_vector, const double u, const double v) const{ 
    return 0;
}


bool MicrofacetMaterial::compute_direct_light() const{ 
    return true;
}

double MicrofacetMaterial::chi(const double x) const{
    return x > 0 ? 1 : 0;
}

double MicrofacetMaterial::D(const vec3& half_vector, const vec3& normal_vector, const double alpha) const{
    double cos_theta = dot_vectors(half_vector, normal_vector);
    double cos_theta_2 = cos_theta * cos_theta;
    double cos_theta_4 = cos_theta_2 * cos_theta_2;
    double tan_theta = sqrt((1.0 - cos_theta_2) / (cos_theta_2));
    double chi_factor = chi(cos_theta);
    double alpha_b_2 = alpha * alpha;
    double fraction = 1.0 / (M_PI * alpha_b_2 * cos_theta_4);
    double exp_factor = std::exp(- (tan_theta * tan_theta) / alpha_b_2);
    return chi_factor * fraction * exp_factor;
}


double MicrofacetMaterial::G1(const vec3& half_vector, const vec3& normal_vector, const vec3& v, const double alpha) const{
    double cos_theta = dot_vectors(half_vector, v);
    double tan_theta = sqrt((1.0 - cos_theta * cos_theta) / (cos_theta * cos_theta));
    double a = 1.0 / (alpha * tan_theta);

    double approximation = a < 1.6 ? (3.535 * a  + 2.181 * a * a) / (1.0 + 2.276 * a + 2.577 * a * a) : 1.0;
    
    return chi(cos_theta / dot_vectors(v, normal_vector)) * approximation;
}

double MicrofacetMaterial::G(const vec3& half_vector, const vec3& normal_vector, const vec3& incident_vector, const vec3& outgoing_vector, const double alpha) const{
    return G1(half_vector, normal_vector, incident_vector, alpha) * G1(half_vector, normal_vector, outgoing_vector, alpha);
}

MicrofacetData MicrofacetMaterial::prepare_microfacet_data(const Hit& hit, const double u, const double v) const{
    double n1;
    double k1;
    double n2;
    double k2;
    vec3 normal_into_interface;
    double incoming_dot_normal = dot_vectors(hit.incident_vector, hit.normal_vector);
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

    double alpha = std::max(roughness_map -> get(u, v), constants::EPSILON);

    vec3 sampled_half_vector = specular_sample(-normal_into_interface, alpha);
            
    double i_dot_h = dot_vectors(hit.incident_vector, sampled_half_vector);

    double F_r = fresnel_multiplier(-i_dot_h, n1, k1, n2, k2, is_dielectric);

    MicrofacetData data;
    data.outside = outside;
    data.alpha = alpha;
    data.eta = n1 / n2;
    data.normal_into_interface = normal_into_interface;
    data.half_vector = sampled_half_vector;
    data.F_r = F_r;
    return data;
}

vec3 MicrofacetMaterial::eval(const Hit& hit, const double u, const double v) const{
    // TODO: This needs outgoing vector if I want to return non-zero values for specular reflection (might result in lower variance, but only if alpha != 0!).
    MicrofacetData data = prepare_microfacet_data(hit, u, v);
    double random_num = random_uniform(0, 1);
    bool reflect_specular = random_num < data.F_r;
    if (reflect_specular){
        return colors::BLACK; // pdf * reflected_color * G(args.sampled_half_vector, args.normal_vector, args.incident_vector, data.outgoing_vector, args.alpha) * args.cosine_factor;
    }
    double transmit = random_uniform(0, 1) > percentage_diffuse_map -> get(u, v);

    if (transmit){
        return colors::BLACK; // pdf * G(args.sampled_half_vector, args.normal_vector, args.incident_vector, data.outgoing_vector, args.alpha) * args.cosine_factor;
    }

    return albedo_map -> get(u, v) / M_PI;
}

vec3 MicrofacetMaterial::specular_sample(const vec3& normal_vector, const double alpha) const{
    double r1 = random_uniform(0, 1);
    double r2 = random_uniform(0, 1);
    double phi = 2 * M_PI * r2;
    double tan_theta2 = - alpha * alpha * std::log(1 - r1);
    double cos_theta2 = 1.0 / (1.0 + tan_theta2);
    
    double cos_theta = sqrt(cos_theta2);
    double sin_theta = sqrt(1 - cos_theta2);

    vec3 x_hat; 
    vec3 y_hat;
    set_perpendicular_vectors(normal_vector, x_hat, y_hat);
    return x_hat * sin_theta * cos(phi) + y_hat * sin_theta * sin(phi) + normal_vector * cos_theta;
}

BrdfData MicrofacetMaterial::sample_diffuse(const MicrofacetSampleArgs& args) const{
    BrdfData data;
    data.outgoing_vector = sample_cosine_hemisphere(args.normal_vector);
    data.brdf_over_pdf = albedo_map -> get(args.u, args.v);
    data.pdf = dot_vectors(args.normal_vector, data.outgoing_vector) / M_PI;
    data.type = DIFFUSE;
    return data;
}

BrdfData MicrofacetMaterial::sample_reflection(const MicrofacetSampleArgs& args) const{
    BrdfData data;
    vec3 reflection_color = is_dielectric ? colors::WHITE : albedo_map -> get(args.u, args.v);
    data.outgoing_vector = reflect_vector(-args.incident_vector, args.sampled_half_vector);
    data.pdf = brdf_pdf(data.outgoing_vector, args.incident_vector, args.normal_vector, args.u, args.v);
    data.brdf_over_pdf = reflection_color * G(args.sampled_half_vector, args.normal_vector, args.incident_vector, data.outgoing_vector, args.alpha) * args.cosine_factor;
    data.type = REFLECTED;
    return data;
}

BrdfData MicrofacetMaterial::sample_transmission(const MicrofacetSampleArgs& args) const{
    vec3 refracted_vector = refract_vector(-args.incident_vector, -args.sampled_half_vector, args.eta);

    if (refracted_vector.length() == 0){
        return sample_reflection(args);
    }

    BrdfData data;
    data.outgoing_vector = refracted_vector;

    data.brdf_over_pdf = G(args.sampled_half_vector, args.normal_vector, args.incident_vector, data.outgoing_vector, args.alpha) * args.cosine_factor;
    data.type = TRANSMITTED;
    data.pdf = brdf_pdf(data.outgoing_vector, args.incident_vector, args.normal_vector, args.u, args.v);
    return data;
}


BrdfData MicrofacetMaterial::sample(const Hit& hit, const double u, const double v) const{
    MicrofacetData data = prepare_microfacet_data(hit, u, v);
            
    double i_dot_h = dot_vectors(hit.incident_vector, data.half_vector);
    double i_dot_n = dot_vectors(hit.incident_vector, data.normal_into_interface);
    double n_dot_h = dot_vectors(data.half_vector, data.normal_into_interface);

    double cosine_factor = std::abs(i_dot_h / (i_dot_n * n_dot_h));

    double random_num = random_uniform(0, 1);
    bool reflect_specular = random_num <= data.F_r;
    
    MicrofacetSampleArgs args;
    args.sampled_half_vector = data.half_vector;
    args.normal_vector = -data.normal_into_interface;
    args.incident_vector = -hit.incident_vector;
    args.cosine_factor = cosine_factor;
    args.u = u;
    args.v = v;
    args.alpha = data.alpha;
    if (reflect_specular){
        return sample_reflection(args);
    }
    else{
        args.eta = data.eta;
        args.outside = data.outside;
        double random_num2 = random_uniform(0, 1);
        bool diffuse = random_num2 < percentage_diffuse_map -> get(u, v);
        return diffuse ? sample_diffuse(args) : sample_transmission(args);
    }
}

double MicrofacetMaterial::brdf_pdf(const vec3& outgoing_vector, const vec3& incident_vector, const vec3& normal_vector, const double u, const double v) const{  
    vec3 half_vector;
    bool reflected = dot_vectors(incident_vector, normal_vector) * dot_vectors(outgoing_vector, normal_vector) < 0;
    if (reflected){
        half_vector = normalize_vector(outgoing_vector + incident_vector);
    }
    else{
        half_vector = normalize_vector(incident_vector - outgoing_vector); // TODO: Add refractive indices here. ###
    }
    return D(half_vector, normal_vector, std::max(roughness_map -> get(u, v), constants::EPSILON)) * dot_vectors(half_vector, normal_vector);
}


MaterialManager::~MaterialManager(){
    for (int i = 0; i < current_idx; i++){
        delete material_array[i];
    }
}

void MaterialManager::add_material(Material* material){
    material_array[current_idx++] = material;
};
