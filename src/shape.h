#ifndef SHAPE_H
#define SHAPE_H

#include "ray.h"

typedef Vector Color;

enum Refl_t
{
    DIFF, SPEC, REFR, GLOSSY, TRANSL 
};

class Shape
{
public:
    virtual double intersect(const Ray &ray) const  = 0;
    virtual Vector get_normal(const Vector &hitpoint) const = 0;
    virtual const Vector & get_position() const  = 0;

    double radius;
    Color emission, color;
    Refl_t refl;
};

#endif // SHAPE_H
