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
    ReflectionType type = ReflectionType::DIFFUSE;
    double pdf;
    vec3 outgoing_vector;
    vec3 brdf_over_pdf;
};

struct MicrofacetSampleArgs {
    double cosine_factor;
    double eta;
    double u;
    double v;
    double alpha;
    bool outside;
    vec3 sampled_half_vector;
    vec3 normal_vector;
    vec3 incident_vector;
};

struct MicrofacetData {
    bool outside;
    double alpha;
    double eta;
    double F_r;
    vec3 normal_into_interface;
    vec3 half_vector;
};

struct MaterialData {
    bool is_light_source = false;
    bool is_dielectric = true;
    double surface_refractive_index = 1.0;
    double extinction_coefficient = 2.0;
    ValueMap1D const* light_intensity_map = nullptr;
    ValueMap1D const* roughness_map = nullptr;
    ValueMap1D const* transparency_map = nullptr;
    ValueMap3D const* albedo_map = nullptr;
    ValueMap3D const* emission_color_map = nullptr;
    Medium const* internal_medium = nullptr;
    Medium const* external_medium = nullptr;
};

class Material {
  public:
    bool is_light_source;
    bool is_dielectric;
    double surface_refractive_index;
    double extinction_coefficient;
    double attenuation_coefficient;
    vec3 absorption_albedo;
    ValueMap1D const* light_intensity_map;
    ValueMap1D const* roughness_map;
    ValueMap1D const* transparency_map;
    ValueMap3D const* albedo_map;
    ValueMap3D const* emission_color_map;
    Medium const* internal_medium;
    Medium const* external_medium;

    Material() {}
    Material(MaterialData data);
    virtual ~Material() {}

    virtual bool allow_direct_light(const double u, const double v) const;
    virtual vec3 eval(const Hit& hit, const vec3& outgoing_vector, const double u, const double v) const;
    virtual BrdfData sample(const Hit& hit, const double u, const double v) const;
    virtual double brdf_pdf(const vec3& outgoing_vector, const vec3& incident_vector, const vec3& normal_vector,
                            const double u, const double v) const;
    vec3 get_light_emittance(const double u, const double v) const;
    bool sample_transparency_map(const double u, const double v) const;
};

class DiffuseMaterial : public Material {
  public:
    using Material::Material;

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

    bool allow_direct_light(const double u, const double v) const override;
    vec3 eval(const Hit& hit, const vec3& outgoing_vector, const double u, const double v) const override;
    BrdfData sample(const Hit& hit, const double u, const double v) const override;
};

class MicrofacetMaterial : public Material {
  public:
    using Material::Material;

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

class MetallicMicrofacetMaterial : public MicrofacetMaterial {
  public:
    MetallicMicrofacetMaterial(MaterialData data);

    vec3 eval(const Hit& hit, const vec3& outgoing_vector, const double u, const double v) const override;
    vec3 sample_outgoing(const vec3& incident_vector, const vec3& normal_vector, const double u, const double v) const;
    BrdfData sample(const Hit& hit, const double u, const double v) const override;
    double brdf_pdf(const vec3& outgoing_vector, const vec3& incident_vector, const vec3& normal_vector, const double u,
                    const double v) const override;
};

class ReflectiveMicrofacetMaterial : public MetallicMicrofacetMaterial {
  public:
    using MetallicMicrofacetMaterial::MetallicMicrofacetMaterial;

    vec3 eval(const Hit& hit, const vec3& outgoing_vector, const double u, const double v) const override;
};

class TransparentMicrofacetMaterial : public MicrofacetMaterial {
  public:
    using MicrofacetMaterial::MicrofacetMaterial;

    vec3 eval(const Hit& hit, const vec3& outgoing_vector, const double u, const double v) const override;
    vec3 sample_outgoing(vec3& half_vector, const vec3& incident_vector, const vec3& normal_vector, const bool outside,
                         const double u, const double v) const;
    BrdfData sample(const Hit& hit, const double u, const double v) const override;
    double brdf_pdf(const vec3& outgoing_vector, const vec3& incident_vector, const vec3& normal_vector, const double u,
                    const double v) const override;
};

#endif
