#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "shape.h"

class Triangle : public Shape
{
public:
    Vector a, b, c;

    Triangle(const Vector &a_, const Vector &b_, const Vector &c_,
             const Color &emission_, const Color &color_, Refl_t refl_) :
            a(a_), b(b_), c(c_)
    {
        Vector temp_v0, temp_v1;
        temp_v0 = b - a;
        temp_v1 = c - a;
        normal = temp_v0.Cross(temp_v1);
        normal = normal.Normalized();
        center = (a+b+c) / 3.0;

		emission = emission_;
		color = color_;
		refl = refl_;

        //farthest point away from center has to be vertice
        double length_to_a = (center - a).Length();
        double length_to_b = (center - b).Length();
        double length_to_c = (center - c).Length();
        radius = std::max(std::max(length_to_a, length_to_b), length_to_c);

    }

    /* triangle-ray intersection */
    double intersect(const Ray &ray) const
    {
        Vector edge1 = b - a;
        Vector edge2 = c - a;
        Vector dir = ray.dir;
        Vector orig = ray.org;
        Vector pvec = dir.Cross(edge2);
        double det = edge1.Dot(pvec);

        if (det == 0)
            return 0.0;

        double invDet = 1.0 / det;
        Vector tvec = orig - a;
        double u = tvec.Dot(pvec) * invDet;

        if (u < 0 || u > 1)
            return 0.0;

        Vector qvec = tvec.Cross(edge1);
        double v = dir.Dot(qvec) * invDet;

        if (v < 0 || u + v > 1)
            return 0.0;

        double t = edge2.Dot(qvec) * invDet;

        if (t <= 0.00000001)
            return 0.0;

        return t;
    }

    const Vector & get_position() const
    {
        return center;
    }

    Vector get_normal(const Vector &hitpoint) const
    {
        return normal;
    }

private:
    Vector center;
    Vector normal;

};

#endif // TRIANGLE_H
