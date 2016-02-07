#ifndef SPHERE_H
#define SPHERE_H

#include "shape.h"

class Sphere : public Shape
{
public:

    Sphere(double radius_, Vector position_, Vector emission_,
           Vector color_, Refl_t refl_)
    {
        radius = radius_;
		color = color_;
		refl = refl_;
		emission = emission_;
		position = position_;
    }

    double intersect(const Ray &ray) const
    {
        /* Check for ray-sphere intersection by solving for t:
            t^2*d.d + 2*t*(o-p).d + (o-p).(o-p) - R^2 = 0 */
        Vector op = position - ray.org;
        double eps = 1e-4;
        double b = op.Dot(ray.dir);
        double radicant = b * b - op.Dot(op) + radius * radius;
        if (radicant < 0.0)
            return 0.0;      /* No intersection */
        else
            radicant = sqrt(radicant);

        double t;
        t = b - radicant;    /* Check smaller root first */
        if (t > eps)
            return t;

        t = b + radicant;
        if (t > eps)          /* Check second root */
            return t;

        return 0.0;          /* No intersection in ray direction */
    }

    const Vector & get_position() const
    {
        return position;
    }

    Vector get_normal(const Vector &hitpoint) const
    {
        return (hitpoint - position).Normalized();  /* Normal at intersection */

    }
private:
    Vector position;
};

#endif // SPHERE_H
