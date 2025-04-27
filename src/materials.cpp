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

    medium = data.medium;
}

Material::~Material(){
    delete albedo_map;
    delete emission_color_map;
    delete light_intensity_map;
    delete roughness_map;
    delete medium;
}

bool Material::allow_direct_light() const { return false; }
bool Material::compute_direct_light() const { return false; }
vec3 Material::eval(const Hit& hit, const vec3& outgoing_vector, const double u, const double v) const{ return vec3(0); }
BrdfData Material::sample(const Hit& hit, const double u, const double v) const{ return BrdfData(); }
double Material::brdf_pdf(const vec3& outgoing_vector, const vec3& incident_vector, const vec3& normal_vector, const double u, const double v) const{ return 0; }

vec3 Material::get_light_emittance(const double u, const double v) const{
    return emission_color_map -> get(u, v) * light_intensity_map -> get(u, v);
}


bool DiffuseMaterial::compute_direct_light() const{ 
    return true; 
}

vec3 DiffuseMaterial::eval(const Hit& hit, const vec3& outgoing_vector, const double u, const double v) const{
    return albedo_map -> get(u, v) / M_PI;
}

BrdfData DiffuseMaterial::sample(const Hit& hit, const double u, const double v) const{
    vec3 outgoing_vector = sample_cosine_hemisphere(hit.normal_vector);
    BrdfData data;
    data.outgoing_vector = outgoing_vector;
    data.brdf_over_pdf = albedo_map -> get(u, v);
    data.pdf = brdf_pdf(outgoing_vector, hit.incident_vector, hit.normal_vector, u, v);
    data.type = DIFFUSE;
    return data;
}

double DiffuseMaterial::brdf_pdf(const vec3& outgoing_vector, const vec3& incident_vector, const vec3& normal_vector, const double u, const double v) const{ 
    return std::max(dot_vectors(outgoing_vector, normal_vector), 0.0) / M_PI; 
}


vec3 ReflectiveMaterial::eval(const Hit& hit, const vec3& outgoing_vector, const double u, const double v) const{
    return colors::BLACK;
}

BrdfData ReflectiveMaterial::sample(const Hit& hit, const double u, const double v) const{
    vec3 outgoing_vector = reflect_vector(hit.incident_vector, hit.normal_vector);
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
vec3 TransparentMaterial::eval(const Hit& hit, const vec3& outgoing_vector, const double u, const double v) const{
    return colors::BLACK;
}

BrdfData TransparentMaterial::sample(const Hit& hit, const double u, const double v) const{
    // TODO: Look at the current medium, use that as refractive index! Update for extinction coefficient!
    double n1, k1, n2, k2;
    if (hit.outside){
        n1 = constants::air_refractive_index;
        k1 = 0;
        n2 = refractive_index;
        k2 = extinction_coefficient;
    }
    else{
        n1 = refractive_index;
        k1 = extinction_coefficient;
        n2 = constants::air_refractive_index;
        k2 = 0;
    }

    vec3 transmitted_vector = refract_vector(hit.incident_vector, -hit.normal_vector, n1 / n2);

    double F_r = 1;
    if (transmitted_vector.length_squared() != 0){
        double cos_incident = -dot_vectors(hit.incident_vector, hit.normal_vector);
        F_r = fresnel_multiplier(cos_incident, n1, k1, n2, k2, is_dielectric);
    }

    double random_num = random_uniform(0, 1);
    bool is_reflected = random_num <= F_r;
    
    BrdfData data;
    if (is_reflected){
        data.type = REFLECTED;
        data.outgoing_vector = reflect_vector(hit.incident_vector, hit.normal_vector);
        data.brdf_over_pdf = is_dielectric ? colors::WHITE : albedo_map -> get(u, v);
    }
    else{
        data.type = TRANSMITTED;
        data.outgoing_vector = transmitted_vector;
        data.brdf_over_pdf = colors::WHITE;
    }

    return data;
}


bool MicrofacetMaterial::compute_direct_light() const{ 
    return true;
}

double MicrofacetMaterial::chi(const double x) const{
    return x > 0 ? 1 : 0;
}

double MicrofacetMaterial::get_alpha(const double u, const double v) const{
    // roughness for 0.01 and above seems to work without nans - perhaps set that as the limit?
    return std::max(roughness_map -> get(u, v), constants::EPSILON);
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
    // Uses convention of v facing away from intersection.
    double cos_theta = dot_vectors(half_vector, v);
    double tan_theta = sqrt((1.0 - cos_theta * cos_theta) / (cos_theta * cos_theta));
    double a = 1.0 / (alpha * tan_theta);

    double approximation = a < 1.6 ? (3.535 * a  + 2.181 * a * a) / (1.0 + 2.276 * a + 2.577 * a * a) : 1.0;
    
    return chi(cos_theta / dot_vectors(v, normal_vector)) * approximation;
}

double MicrofacetMaterial::G(const vec3& half_vector, const vec3& normal_vector, const vec3& incident_vector, const vec3& outgoing_vector, const double alpha) const{
    return G1(half_vector, normal_vector, -incident_vector, alpha) * G1(half_vector, normal_vector, outgoing_vector, alpha);
}

vec3 MicrofacetMaterial::sample_half_vector(const vec3& normal_vector, const double alpha) const{
    // Samples half_angle_vector
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

double MicrofacetMaterial::diffuse_pdf(const vec3& outgoing_vector, const vec3& normal_vector) const{
    return std::max(dot_vectors(normal_vector, outgoing_vector) / M_PI, 0.0);
}

double MicrofacetMaterial::specular_pdf(const vec3& outgoing_vector, const vec3& incident_vector, const vec3& normal_vector, const double u, const double v) const{
    vec3 half_vector = normalize_vector(outgoing_vector - incident_vector); 

    double half_vector_pdf = D(half_vector, normal_vector, get_alpha(u, v)) * dot_vectors(half_vector, normal_vector);
    return std::max(half_vector_pdf / (4.0 * dot_vectors(outgoing_vector, half_vector)), 0.0); 
}


vec3 GlossyMaterial::eval(const Hit& hit, const vec3& outgoing_vector, const double u, const double v) const{
    // compute fresnel glossy is kind of redundan. New method that returns refractive indices etc.
    vec3 half_vector = normalize_vector(outgoing_vector - hit.incident_vector);

    double n1, k1, n2, k2;
    if (hit.outside){
        n1 = constants::air_refractive_index;
        k1 = 0;
        n2 = refractive_index;
        k2 = extinction_coefficient;
    }
    else{
        n1 = refractive_index;
        k1 = extinction_coefficient;
        n2 = constants::air_refractive_index;
        k2 = 0;
    }

    double i_dot_h = -dot_vectors(hit.incident_vector, half_vector);

    double F_r = schlick_fresnel(i_dot_h, n1, k1, n2, k2);

    double R_0_sqrt = (n1 - n2) / (n1 + n2);
    double R_0 = R_0_sqrt * R_0_sqrt;
    double fac1 = std::min((1.0-dot_vectors(hit.normal_vector, -hit.incident_vector)/2.0), 1.0);
    double fac2 = std::min((1.0-dot_vectors(hit.normal_vector, outgoing_vector)/2.0), 1.0);
    vec3 diffuse = albedo_map -> get(u, v) * 28.0 / (23.0 * M_PI) * (1 - R_0) * (1-fac1*fac1*fac1*fac1*fac1) * (1-fac2*fac2*fac2*fac2*fac2);

    double alpha = get_alpha(u, v);
    vec3 reflection_color = is_dielectric ? colors::WHITE : albedo_map -> get(u, v);
    double d_factor = D(half_vector, hit.normal_vector, alpha) * dot_vectors(half_vector, hit.normal_vector);
    double g_factor = G(half_vector, hit.normal_vector, hit.incident_vector, outgoing_vector, alpha);
    double denom_factor = -1.0 / (4.0 * dot_vectors(hit.incident_vector, hit.normal_vector) * dot_vectors(hit.normal_vector, outgoing_vector));
    vec3 specular = reflection_color * F_r * d_factor * g_factor * denom_factor;
    return diffuse + specular;
}

vec3 GlossyMaterial::sample_outgoing(const vec3& incident_vector, const vec3& normal_vector, const double u, const double v) const{
    double rand_1 = random_uniform(0, 1);
    if (rand_1 <= 0.5){
        return sample_cosine_hemisphere(normal_vector);
    }
    else{
        vec3 half_vector = sample_half_vector(normal_vector, get_alpha(u, v));
        return reflect_vector(incident_vector, half_vector);
    }
}

BrdfData GlossyMaterial::sample(const Hit& hit, const double u, const double v) const{
    BrdfData brdf_data;
    brdf_data.outgoing_vector = sample_outgoing(hit.incident_vector, hit.normal_vector, u, v);
    brdf_data.pdf = brdf_pdf(brdf_data.outgoing_vector, hit.incident_vector, hit.normal_vector, u, v);
    brdf_data.brdf_over_pdf = brdf_data.pdf == 0 ? 0 : eval(hit, brdf_data.outgoing_vector, u, v) * dot_vectors(brdf_data.outgoing_vector, hit.normal_vector) / brdf_data.pdf;
    brdf_data.type = DIFFUSE;
    // TODO: The above cos term, should it be normal vector or half vector? I think maybe it was half vector before but I changed it? But test anyways.
    return brdf_data;
}

double GlossyMaterial::brdf_pdf(const vec3& outgoing_vector, const vec3& incident_vector, const vec3& normal_vector, const double u, const double v) const { 
    return 0.5 * (diffuse_pdf(outgoing_vector, normal_vector) + specular_pdf(outgoing_vector, incident_vector, normal_vector, u, v));
}


vec3 MetallicMicrofacet::eval(const Hit& hit, const vec3& outgoing_vector, const double u, const double v) const{
    // compute fresnel glossy is kind of redundan. New method that returns refractive indices etc.
    vec3 half_vector = normalize_vector(outgoing_vector - hit.incident_vector);

    double n1, k1, n2, k2;
    double incoming_dot_normal = dot_vectors(hit.incident_vector, hit.normal_vector); // TODO: This is wrong, use hit variable instead.
    bool outside = incoming_dot_normal <= 0.0;
    if (outside){
        n1 = constants::air_refractive_index;
        k1 = 0;
        n2 = refractive_index;
        k2 = extinction_coefficient;
    }
    else{
        n1 = refractive_index;
        k1 = extinction_coefficient;
        n2 = constants::air_refractive_index;
        k2 = 0;
    }

    double i_dot_h = -dot_vectors(hit.incident_vector, half_vector);

    double F_r = fresnel_multiplier(i_dot_h, n1, k1, n2, k2, is_dielectric);

    double alpha = get_alpha(u, v);
    vec3 reflection_color = is_dielectric ? colors::WHITE : albedo_map -> get(u, v);
    double d_factor = D(half_vector, hit.normal_vector, alpha) * dot_vectors(half_vector, hit.normal_vector);
    double g_factor = G(half_vector, hit.normal_vector, hit.incident_vector, outgoing_vector, alpha);
    double denom_factor = -1.0 / (4.0 * dot_vectors(hit.incident_vector, hit.normal_vector) * dot_vectors(hit.normal_vector, outgoing_vector));
    vec3 specular = reflection_color * F_r * d_factor * g_factor * denom_factor;
    return specular;
}

vec3 MetallicMicrofacet::sample_outgoing(const vec3& incident_vector, const vec3& normal_vector, const double u, const double v) const{
    vec3 half_vector = sample_half_vector(normal_vector, get_alpha(u, v));
    return reflect_vector(incident_vector, half_vector);
}

BrdfData MetallicMicrofacet::sample(const Hit& hit, const double u, const double v) const{
    BrdfData brdf_data;
    brdf_data.outgoing_vector = sample_outgoing(hit.incident_vector, hit.normal_vector, u, v);
    brdf_data.pdf = brdf_pdf(brdf_data.outgoing_vector, hit.incident_vector, hit.normal_vector, u, v);
    brdf_data.brdf_over_pdf = brdf_data.pdf == 0 ? 0 : eval(hit, brdf_data.outgoing_vector, u, v) * dot_vectors(brdf_data.outgoing_vector, hit.normal_vector) / brdf_data.pdf;
    brdf_data.type = DIFFUSE;
    // TODO: The above cos term, should it be normal vector or half vector? I think maybe it was half vector before but I changed it? But test anyways.
    return brdf_data;
}

double MetallicMicrofacet::brdf_pdf(const vec3& outgoing_vector, const vec3& incident_vector, const vec3& normal_vector, const double u, const double v) const {
    return specular_pdf(outgoing_vector, incident_vector, normal_vector, u, v);
}


bool TransparentMicrofacetMaterial::compute_direct_light() const{ 
    return false;
}

vec3 TransparentMicrofacetMaterial::eval(const Hit& hit, const vec3& outgoing_vector, const double u, const double v) const{
    return vec3(0.0);
};

vec3 TransparentMicrofacetMaterial::sample_outgoing(vec3& half_vector, const vec3& incident_vector, const vec3& normal_vector, const bool outside, const double u, const double v) const{
    double n1, k1, n2, k2;
    if (outside){
        n1 = constants::air_refractive_index;
        k1 = 0;
        n2 = refractive_index;
        k2 = extinction_coefficient;
    }
    else{
        n1 = refractive_index;
        k1 = extinction_coefficient;
        n2 = constants::air_refractive_index;
        k2 = 0;
    }

    half_vector = sample_half_vector(normal_vector, get_alpha(u, v));

    double i_dot_h = -dot_vectors(incident_vector, half_vector);
    double F_r = fresnel_multiplier(i_dot_h, n1, k1, n2, k2, is_dielectric);

    vec3 refracted_vector = refract_vector(incident_vector, -half_vector, n1 / n2);
    double rand_num = random_uniform(0, 1);
    if (rand_num <= F_r || refracted_vector.length_squared() == 0){
        return reflect_vector(incident_vector, half_vector);
    }

    return refracted_vector;
}

BrdfData TransparentMicrofacetMaterial::sample(const Hit& hit, const double u, const double v) const{
    BrdfData brdf_data;

    vec3 half_vector;
    brdf_data.outgoing_vector = sample_outgoing(half_vector, hit.incident_vector, hit.normal_vector, hit.outside, u, v);
    
    double cosine_factor = dot_vectors(hit.incident_vector, half_vector) / (dot_vectors(hit.incident_vector, hit.normal_vector) * dot_vectors(half_vector, hit.normal_vector));
    brdf_data.brdf_over_pdf = G(half_vector, hit.normal_vector, hit.incident_vector, brdf_data.outgoing_vector, get_alpha(u, v)) * cosine_factor;
    brdf_data.type = REFLECTED; // Change to just specular?
    return brdf_data;
}

double TransparentMicrofacetMaterial::brdf_pdf(const vec3& outgoing_vector, const vec3& incident_vector, const vec3& normal_vector, const double u, const double v) const{
    return 0.0;
}

MaterialManager::~MaterialManager(){
    for (int i = 0; i < current_idx; i++){
        delete material_array[i];
    }
}

void MaterialManager::add_material(Material* material){
    material_array[current_idx++] = material;
};
