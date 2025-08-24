#ifndef OBJECTS_H
#define OBJECTS_H

#include "vec3.h"
#include "utils.h"
#include "constants.h"
#include "colors.h"
#include "materials.h"

class Material;

struct BrdfData;

class MediumStack;

class Object{
    public:
        Material* material;
        double area;
        int primitive_ID; // Used when object belongs to an ObjectUnion.
        Object(){}
        Object(Material* _material);

        virtual vec3 max_axis_point() const;
        virtual vec3 min_axis_point() const;
        virtual vec3 compute_centroid() const;
        virtual vec3 get_UV(const vec3& point) const;
        virtual Material* get_material(const int primitive_ID) const;
        virtual bool is_light_source() const;
        virtual vec3 eval(const Hit& hit, const vec3& outgoing_vector) const;
        vec3 sample_direct(const Hit& hit, Object** objects, const int number_of_objects, const MediumStack& current_medium_stack) const;
        virtual BrdfData sample(const Hit& hit) const;
        virtual double brdf_pdf(const vec3& outgoing_vector, const Hit& hit) const;
        virtual vec3 get_light_emittance(const Hit& hit) const;
        virtual bool find_closest_object_hit(Hit& hit, Ray& ray) const;
        virtual vec3 get_normal_vector(const vec3& surface_point, const int primitive_ID) const;
        virtual vec3 generate_random_surface_point() const;
        double area_to_angle_PDF_factor(const vec3& surface_point, const vec3& intersection_point, const int primitive_ID) const;
        virtual double light_pdf(const vec3& surface_point, const vec3& intersection_point, const int primitive_id) const;
        virtual vec3 random_light_point(const vec3& intersection_point, double& inverse_PDF) const;
};


class Sphere: public Object{
    public:
        Sphere(){}
        Sphere(const vec3& _position, const double _radius);
        Sphere(const vec3& _position, const double _radius, Material*_material);

        vec3 get_UV(const vec3& point) const override;
        bool find_closest_object_hit(Hit& hit, Ray& ray) const override;
        vec3 get_normal_vector(const vec3& surface_point, const int primitive_ID) const override;
        vec3 generate_random_surface_point() const override;
        double light_pdf(const vec3& surface_point, const vec3& intersection_point, const int primitive_id) const override;
        vec3 random_light_point(const vec3& intersection_point, double& inverse_PDF) const override;

    private:
        vec3 position;
        double radius;
};


class Plane: public Object{
    public:
        Plane(){}
        Plane(const vec3& _position, const vec3& _v1, const vec3& _v2, Material*_material);

        vec3 get_UV(const vec3& point) const override;
        bool compute_distance_in_centered_system(const vec3& starting_point, const Ray& ray, double& distance) const;
        bool find_closest_object_hit(Hit& hit, Ray& ray) const override;
        vec3 get_normal_vector(const vec3& surface_point, const int primitive_ID) const override;
        double light_pdf(const vec3& surface_point, const vec3& intersection_point, const int primitive_id) const override;

    protected:
        vec3 position;
        vec3 v1;
        vec3 v2;
        vec3 normal_vector;
};


class Rectangle: public Plane{
    public:
        Rectangle(){}
        Rectangle(const vec3& _position, const vec3& _v1, const vec3& _v2, const double _L1, const double _L2, Material*_material);

        vec3 get_UV(const vec3& point) const override;
        bool find_closest_object_hit(Hit& hit, Ray& ray) const override;
        double light_pdf(const vec3& surface_point, const vec3& intersection_point, const int primitive_id) const override;
        vec3 generate_random_surface_point() const override;

    private:
        double L1;
        double L2;
};


class Triangle: public Object{
    public:
        Triangle(){}
        Triangle(const vec3& _p1, const vec3& _p2, const vec3& _p3, Material*_material);

        vec3 max_axis_point() const override;
        vec3 min_axis_point() const override;
        vec3 compute_centroid() const override;
        void set_vertex_UV(const vec3& _uv1, const vec3& _uv2, const vec3& _uv3);
        void set_vertex_normals(const vec3& _n1, const vec3& _n2, const vec3& _n3);
        vec3 get_normal_vector_smoothed(const vec3& surface_point, const int primitive_ID) const;
        vec3 get_normal_vector(const vec3& surface_point, const int primitive_ID) const override;
        vec3 compute_barycentric(const vec3& point) const;
        vec3 get_UV(const vec3& point) const override;
        bool find_closest_object_hit(Hit& hit, Ray& ray) const override;
        vec3 generate_random_surface_point() const override;

    private:
        vec3 position;
        vec3 normal_vector;
        vec3 p1;
        vec3 p2;
        vec3 p3;
        vec3 v1;
        vec3 v2;
        double x1;
        double y1;
        double x2;
        double y2;
        double x3;
        double y3;
        double det_T;

        vec3 uv1;
        vec3 uv2;
        vec3 uv3;
        vec3 n1;
        vec3 n2;
        vec3 n3;

        bool smooth_shaded = false;
};


bool find_closest_hit(Hit& closest_hit, Ray& ray, Object** objects, const int number_of_objects);
int sample_random_light(Object** objects, const int number_of_objects, int& number_of_light_sources);

vec3 direct_lighting(const vec3& point, Object** objects, const int number_of_objects, vec3& sampled_direction, const MediumStack& current_medium_stack);
double mis_weight(const int n_a, const double pdf_a, const int n_b, const double pdf_b);
vec3 sample_light(const Hit& hit, Object** objects, const int number_of_objects, const MediumStack& current_medium_stack, const bool is_scatter);


#endif
