#ifndef MATERIALS_H
#define MATERIALS_H

#include "vec3.h"
#include "colors.h"
#include "utils.h"
#include "constants.h"
#include "valuemap.h"
#include "medium.h"

class Medium;

struct BrdfData {
    vec3 outgoing_vector;
    vec3 brdf_over_pdf;
    int type = DIFFUSE;
    double pdf;
};

struct MicrofacetSampleArgs {
    vec3 sampled_half_vector;
    vec3 normal_vector;
    vec3 incident_vector;
    double cosine_factor;
    double eta;
    double u;
    double v;
    double alpha;
    bool outside;
};

struct MicrofacetData {
    bool outside;
    double alpha;
    double eta;
    vec3 normal_into_interface;
    vec3 half_vector;
    double F_r;
};

struct MaterialData {
    ValueMap3D* albedo_map = nullptr;
    double refractive_index = 1;
    double extinction_coefficient = 0;
    ValueMap3D* emission_color_map = nullptr;
    ValueMap1D* light_intensity_map = nullptr;
    bool is_dielectric = true;
    ValueMap1D* roughness_map = nullptr;
    bool is_light_source = false;
    Medium* medium = nullptr;
};

class Material {
  public:
    ValueMap3D* albedo_map;
    double refractive_index;
    double attenuation_coefficient;
    vec3 absorption_albedo;
    ValueMap3D* emission_color_map;
    ValueMap1D* light_intensity_map;
    bool is_dielectric;
    double extinction_coefficient;
    bool is_light_source;
    ValueMap1D* roughness_map;
    Medium* medium;

    Material() {}
    Material(MaterialData data);
    ~Material();

    virtual bool allow_direct_light() const;
    virtual bool compute_direct_light() const;
    virtual vec3 eval(const Hit& hit, const vec3& outgoing_vector, const double u, const double v) const;
    virtual BrdfData sample(const Hit& hit, const double u, const double v) const;
    virtual double brdf_pdf(const vec3& outgoing_vector, const vec3& incident_vector, const vec3& normal_vector,
                            const double u, const double v) const;
    vec3 get_light_emittance(const double u, const double v) const;
};

class DiffuseMaterial : public Material {
  public:
    using Material::Material;

    bool compute_direct_light() const override;
    vec3 eval(const Hit& hit, const vec3& outgoing_vector, const double u, const double v) const override;
    BrdfData sample(const Hit& hit, const double u, const double v) const override;
    double brdf_pdf(const vec3& outgoing_vector, const vec3& incident_vector, const vec3& normal_vector, const double u,
                    const double v) const override;
};

class ReflectiveMaterial : public Material {
  public:
    using Material::Material;

    vec3 eval(const Hit& hit, const vec3& outgoing_vector, const double u, const double v) const override;
    BrdfData sample(const Hit& hit, const double u, const double v) const override;
    double brdf_pdf(const vec3& outgoing_vector, const vec3& incident_vector, const vec3& normal_vector, const double u,
                    const double v) const override;
};

class TransparentMaterial : public Material {
  public:
    using Material::Material;

    bool allow_direct_light() const override;
    vec3 eval(const Hit& hit, const vec3& outgoing_vector, const double u, const double v) const override;
    BrdfData sample(const Hit& hit, const double u, const double v) const override;
};

class MicrofacetMaterial : public Material {
  public:
    using Material::Material;

    bool compute_direct_light() const override;
    double chi(const double x) const;
    double get_alpha(const double u, const double v) const;
    double D(const vec3& half_vector, const vec3& normal_vector, const double alpha) const;
    double G1(const vec3& half_vector, const vec3& normal_vector, const vec3& v, const double alpha) const;
    double G(const vec3& half_vector, const vec3& normal_vector, const vec3& incident_vector,
             const vec3& outgoing_vector, const double alpha) const;
    vec3 sample_half_vector(const vec3& normal_vector, const double alpha) const;
    BrdfData sample_transmission(const MicrofacetSampleArgs& args) const;
    double diffuse_pdf(const vec3& outgoing_vector, const vec3& normal_vector) const;
    double specular_pdf(const vec3& outgoing_vector, const vec3& incident_vector, const vec3& normal_vector,
                        const double u, const double v) const;
};

class GlossyMaterial : public MicrofacetMaterial {
  public:
    using MicrofacetMaterial::MicrofacetMaterial;
    vec3 eval(const Hit& hit, const vec3& outgoing_vector, const double u, const double v) const override;
    vec3 sample_outgoing(const vec3& incident_vector, const vec3& normal_vector, const double u, const double v) const;
    BrdfData sample(const Hit& hit, const double u, const double v) const override;
    double brdf_pdf(const vec3& outgoing_vector, const vec3& incident_vector, const vec3& normal_vector, const double u,
                    const double v) const override;
};

class MetallicMicrofacet : public MicrofacetMaterial {
  public:
    using MicrofacetMaterial::MicrofacetMaterial;

    vec3 eval(const Hit& hit, const vec3& outgoing_vector, const double u, const double v) const override;
    vec3 sample_outgoing(const vec3& incident_vector, const vec3& normal_vector, const double u, const double v) const;
    BrdfData sample(const Hit& hit, const double u, const double v) const override;
    double brdf_pdf(const vec3& outgoing_vector, const vec3& incident_vector, const vec3& normal_vector, const double u,
                    const double v) const override;
};

class TransparentMicrofacetMaterial : public MicrofacetMaterial {
  public:
    using MicrofacetMaterial::MicrofacetMaterial;

    bool compute_direct_light() const override;
    vec3 eval(const Hit& hit, const vec3& outgoing_vector, const double u, const double v) const override;
    vec3 sample_outgoing(vec3& half_vector, const vec3& incident_vector, const vec3& normal_vector, const bool outside,
                         const double u, const double v) const;
    BrdfData sample(const Hit& hit, const double u, const double v) const override;
    double brdf_pdf(const vec3& outgoing_vector, const vec3& incident_vector, const vec3& normal_vector, const double u,
                    const double v) const override;
};

class MaterialManager {
    const int MAX_MATERIALS = 100;
    Material** material_array = new Material*[MAX_MATERIALS];
    int current_idx = 0;

  public:
    MaterialManager() {}
    ~MaterialManager();

    void add_material(Material* material);
};

#endif
