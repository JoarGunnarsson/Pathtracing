#ifndef MATERIALS_H
#define MATERIALS_H

#include "vec3.h"
#include "colors.h"
#include "utils.h"
#include "constants.h"
#include "valuemap.h"
#include "medium.h"


class Medium;


struct BrdfData{
    vec3 outgoing_vector;
    vec3 brdf_multiplier;
    int type = DIFFUSE;
};


struct MicrofacetSampleArgs{
    vec3 sampled_half_vector;
    vec3 normal_vector;
    vec3 incident_vector;
    float cosine_factor;
    float eta;
    float u;
    float v;
    float alpha;
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
    ValueMap3D* emission_color_map = nullptr;
    ValueMap1D* light_intensity_map = nullptr;
    bool is_dielectric = true;
    ValueMap1D* roughness_map = nullptr; 
    ValueMap1D* percentage_diffuse_map = nullptr;
    bool is_light_source = false;
    Medium* medium = nullptr;
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
        Material(MaterialData data);
        ~Material();
        
    virtual vec3 eval(const Hit& hit, const float u, const float v);
    virtual BrdfData sample(const Hit& hit, const float u, const float v);
    vec3 get_light_emittance(const float u, const float v);
};


class DiffuseMaterial : public Material{
    public:
        using Material::Material;

    vec3 eval(const Hit& hit, const float u, const float v) override;
    BrdfData sample(const Hit& hit, const float u, const float v) override;
};


class ReflectiveMaterial : public Material{
    public:
        using Material::Material;

    vec3 eval(const Hit& hit, const float u, const float v) override;
    BrdfData sample(const Hit& hit, const float u, const float v) override;
};


class TransparentMaterial : public Material{
    public:
        using Material::Material;

    vec3 eval(const Hit& hit, const float u, const float v) override;
    BrdfData sample(const Hit& hit, const float u, const float v) override;
};


class MicrofacetMaterial : public Material{
    public:
        using Material::Material;

    float chi(const float x);
    float G1(const vec3& half_vector, const vec3& normal_vector, const vec3& v, const float alpha);
    float G(const vec3& half_vector, const vec3& normal_vector, const vec3& incident_vector, const vec3& outgoing_vector, const float alpha);
    MicrofacetData prepare_microfacet_data(const Hit& hit, const float u, const float v);
    vec3 eval(const Hit& hit, const float u, const float v) override;
    vec3 specular_sample(const vec3& normal_vector, const float alpha);
    BrdfData sample_diffuse(const MicrofacetSampleArgs& args);
    BrdfData sample_reflection(const MicrofacetSampleArgs& args);
    BrdfData sample_transmission(const MicrofacetSampleArgs& args);
    BrdfData sample(const Hit& hit, const float u, const float v) override;
};


class MaterialManager{
    const int MAX_MATERIALS = 100;
    Material** material_array = new Material*[MAX_MATERIALS];
    int current_idx = 0;

    public:
        MaterialManager(){}
        ~MaterialManager();

        void add_material(Material* material);
};

#endif
