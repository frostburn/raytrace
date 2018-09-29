#include <cmath>

#include "raytrace/quaternion.h"
#include "raytrace/power3d.h"

quaternion pow3d(const quaternion& v, const int& n) {
    real x = v.x;
    real y = v.y;
    real z = v.z;

    real r = sqrt(x*x + y*y);
    if (r != 0) {
        x /= r;
        y /= r;
    }

    real temp;
    real res_x = x;
    real res_y = y;
    real res_z = z;
    real res_r = r;
    for (int i = 1; i < n; ++i) {
        temp = res_x;
        res_x = res_x*x - res_y*y;
        res_y = temp*y + res_y*x;

        temp = res_r;
        res_r = res_r*r - res_z*z;
        res_z = temp*z + res_z*r;
    }
    return (quaternion){0, res_x*res_r, res_y*res_r, res_z};
}

// std::tuple<quaternion, quaternion, quaternion, quaternion> pow3d_grad(const quaternion& v, const int& n) {
quaternion* pow3d_grad(const quaternion& v, const int& n) {
    real x = v.x;
    real y = v.y;
    real z = v.z;

    // Jacobian matrix with the real part used for temporary storage
    quaternion dx = {0, 1, 0, 0};
    quaternion dy = {0, 0, 1, 0};
    quaternion dz = {0, 0, 0, 1};

    quaternion ddx = {0, 1, 0, 0};
    quaternion ddy = {0, 0, 1, 0};

    real r = sqrt(x*x + y*y);
    if (r != 0) {
        ddx.r = x / r;
        ddx.x = (ddx.x*r - x*ddx.r) / (r*r);
        ddx.y = (ddx.y*r - y*ddx.r) / (r*r);

        ddy.r = y / r;
        ddy.x = (ddy.x*r - x*ddy.r) / (r*r);
        ddy.y = (ddy.y*r - y*ddy.r) / (r*r);

        x /= r;
        y /= r;
    }
    dx.r = ddx.r;
    dx.x = ddx.x;
    dx.y = ddx.y;

    dy.r = ddy.r;
    dy.x = ddy.x;
    dy.y = ddy.y;

    real temp;
    real res_x = x;
    real res_y = y;
    real res_z = z;
    real res_r = r;
    for (int i = 1; i < n; ++i) {
        temp = dx.x;
        dx.x = dx.x*x + res_x*ddx.x - dx.y*y - res_y*ddx.y;
        dx.y = temp*y + res_x*ddx.y + dx.y*x + res_y*ddx.x;

        temp = dy.x;
        dy.x = dy.x*x + res_x*ddy.x - dy.y*y - res_y*ddy.y;
        dy.y = temp*y + res_x*ddy.y + dy.y*x + res_y*ddy.x;

        temp = res_x;
        res_x = res_x*x - res_y*y;
        res_y = temp*y + res_y*x;

        temp = dx.r;
        dx.r = dx.r*r + res_r*ddx.r - dx.z*z;
        dx.z = temp*z + dx.z*r + res_z*ddx.r;

        temp = dy.r;
        dy.r = dy.r*r + res_r*ddy.r - dy.z*z;
        dy.z = temp*z + dy.z*r + res_z*ddy.r;

        temp = dz.r;
        dz.r = dz.r*r - dz.z*z - res_z;
        dz.z = temp*z + res_r + dz.z*r;

        temp = res_r;
        res_r = res_r*r - res_z*z;
        res_z = temp*z + res_z*r;
    }
    dx.x = dx.r*res_x + dx.x*res_r;
    dx.y = dx.r*res_y + dx.y*res_r;
    dx.r = 0;

    dy.x = dy.r*res_x + dy.x*res_r;
    dy.y = dy.r*res_y + dy.y*res_r;
    dy.r = 0;

    dz.x = res_x * dz.r;
    dz.y = res_y * dz.r;
    dz.r = 0;

    quaternion *res = new quaternion[4];
    res[0] = (quaternion){0, res_x*res_r, res_y*res_r, res_z};
    res[1] = dx;
    res[2] = dy;
    res[3] = dz;
    return res;
    // return {
    //     res,
    //     dx, dy, dz
    // };
}
