#ifndef RAY_H
#define RAY_H

#include "vector.h"


class Ray
{
public:
    Vector org, dir;

    /* Origin and direction */
    Ray(const Vector org_, const Vector &dir_) : org(org_), dir(dir_)
    { }
	Ray()
	{
		org = Vector(0,0,0);
		dir = Vector(0,0,0);
	}
};

#endif // RAY_H
