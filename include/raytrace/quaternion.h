#ifndef RAYTRACE_QUATERNION_H_GUARD
#define RAYTRACE_QUATERNION_H_GUARD

#include <iostream>

typedef double real;

struct quaternion
{
    real r;
    real x;
    real y;
    real z;
};

typedef quaternion color;

std::ostream& operator<<(std::ostream& os, const quaternion& v);

real norm(quaternion v);

quaternion conjugate(quaternion v);

quaternion operator -(const quaternion& v);

quaternion operator *(const quaternion& a, const real b);

quaternion operator *(const real b, const quaternion& a);

quaternion operator /(const quaternion& a, const real b);

quaternion operator /(const real a, const quaternion& b);

quaternion operator +(const quaternion& a, const quaternion& b);

quaternion operator +(const real& a, const quaternion& b);

quaternion operator -(const quaternion& a, const quaternion& b);

quaternion operator *(const quaternion& a, const quaternion& b);

quaternion operator /(const quaternion& a, const quaternion& b);

quaternion infinity ();

real dot(const quaternion& a, const quaternion& b);

quaternion project(const quaternion& a, const quaternion& b);

quaternion cross_align(const quaternion& a, const quaternion& b);

color clip_color(const color& a);

#endif
