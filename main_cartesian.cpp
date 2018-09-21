#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <vector>
#include <memory>
#include <random>
#include <cassert>

#include "raytrace/quaternion.h"
#include "raytrace/ray_traceable.h"
#include "raytrace/gradient_traceable.h"
#include "raytrace/sphere.h"
#include "raytrace/clifford_torus.h"
#include "raytrace/trace.h"

using namespace std;

int main ()
{
    quaternion camera_pos = {0, 0, 1, -2.2};
    quaternion look_at = {0, 0, 0, 0};
    quaternion up = {0, 0, 1, 0};
    real view_width = 1;

    // Hacky cross products
    quaternion vleft = (look_at - camera_pos) * up;
    vleft.r = 0;
    up = vleft * (look_at - camera_pos);
    up.r = 0;

    vleft = vleft / norm(vleft) * view_width;
    up = up / norm(up) * view_width;

    shared_ptr<Sphere> sphere = make_shared<Sphere>();
    sphere->location = {0, 0, 0, 0};
    // sphere->scale = sphere->scale * 1.4;
    sphere->scale.y = 0.5;
    sphere->reflective = false;

    shared_ptr<Sphere> another_sphere = make_shared<Sphere>();
    another_sphere->location = {0, 0.2, 0.4, 0.5};
    another_sphere->scale = another_sphere->scale * 0.2;
    another_sphere->reflective = false;

    shared_ptr<CliffordTorus> cylinder = make_shared<CliffordTorus>();
    cylinder->location = {0, -0.3, 0.5, -0.2};
    cylinder->scale = cylinder->scale * 0.3;
    cylinder->left_transform = (quaternion){1, 0.3, 0.3, 0.1};
    cylinder->right_transform = 1.0 / cylinder->left_transform;
    cylinder->reflective = false;

    vector<shared_ptr<RayTraceable>> objects;
    objects.push_back(sphere);
    // objects.push_back(another_sphere);
    // objects.push_back(cylinder);

    int width = 50;
    int height = 50;

    default_random_engine generator;
    generator.seed(random_device()());
    normal_distribution<real> distribution(0.0, 0.2 / (real) width);

    cout << "P2" << endl;
    cout << width << " " << height << endl;
    cout << 255 << endl;
    for (int j = 0; j < height; ++j) {
        cerr << (j * 100) / height << "%" << endl;
        real y_ = 1 - 2 * (j / (real) height);
        for (int i = 0; i < width; ++i) {
            real x = 1 - 2 * (i / (real) width) + distribution(generator);
            real y = y_ + distribution(generator);
            quaternion target = look_at + x * vleft + y * up;
            color pixel = raytrace(camera_pos, target, 10, objects);
            cout << setw(4) << left << (int) (pixel.x * 255);
        }
        cout << endl;
    }
    return EXIT_SUCCESS;
}
