#ifndef VECTOR_H
#define VECTOR_H

#include <cmath>

class Vector
{
public:
    double x, y, z;

    /* Position XYZ or color RGB */

    Vector(const Vector &b) : x(b.x), y(b.y), z(b.z)
    { }

    Vector(double x_ = 0, double y_ = 0, double z_ = 0) : x(x_), y(y_), z(z_)
    { }

    Vector operator+(const Vector &b) const
    {
        return Vector(x + b.x, y + b.y, z + b.z);
    }

    Vector operator-(const Vector &b) const
    {
        return Vector(x - b.x, y - b.y, z - b.z);
    }

    Vector operator/(double c) const
    {
        return Vector(x / c, y / c, z / c);
    }

    Vector operator*(double c) const
    {
        return Vector(x * c, y * c, z * c);
    }

    friend Vector operator*(double c, const Vector &b)
    {
        return b * c;
    }

    Vector MultComponents(const Vector &b) const
    {
        return Vector(x * b.x, y * b.y, z * b.z);
    }

    Vector Normalized() const
    {
        return Vector(x, y, z) / sqrt(x * x + y * y + z * z);
    }

    const double Dot(const Vector &b) const
    {
        return x * b.x + y * b.y + z * b.z;
    }

    const Vector Cross(const Vector &b) const
    {
        return Vector((y * b.z) - (z * b.y),
                      (z * b.x) - (x * b.z),
                      (x * b.y) - (y * b.x));
    }

    const double LengthSquared() const
    {
        return x * x + y * y + z * z;
    }

    const double Length() const
    {
        return sqrt(LengthSquared());
    }

    const double Max()
    {
        return fmax(x, fmax(x, y));
    }

    Vector &clamp()
    {
        x = x < 0 ? 0.0 : x > 1.0 ? 1.0 : x;
        y = y < 0 ? 0.0 : y > 1.0 ? 1.0 : y;
        z = z < 0 ? 0.0 : z > 1.0 ? 1.0 : z;
        return *this;
    }
};
#endif // VECTOR_H
