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
#include "raytrace/hyper_torus.h"
#include "raytrace/trace.h"

using namespace std;

int main ()
{
    quaternion camera_pos = {0, 1, 0.5, -3.2};
    quaternion look_at = {0, 0, 0.3, 0};
    quaternion up = {0, 0, 1, 0};
    real view_width = 1.5;

    // Hacky cross products
    quaternion vleft = (look_at - camera_pos) * up;
    vleft.r = 0;
    up = vleft * (look_at - camera_pos);
    up.r = 0;

    vleft = vleft / norm(vleft) * view_width;
    up = up / norm(up) * view_width;

    shared_ptr<Sphere> sphere = make_shared<Sphere>();
    sphere->location = {0, 0, -0.2, 0};
    sphere->scale = sphere->scale * 1.4;
    sphere->scale.y = 0.3;
    sphere->reflection = 0.9;

    shared_ptr<Sphere> another_sphere = make_shared<Sphere>();
    another_sphere->location = {0, 0.2, 0.4, 0.5};
    another_sphere->scale = another_sphere->scale * 0.2;

    shared_ptr<CliffordTorus> cylinder = make_shared<CliffordTorus>();
    cylinder->pigment_a = {0, 0.5, 0.4, 0.3};
    cylinder->pigment_b = {0, 0.3, 0.5, 0.2};
    cylinder->pigment_c = {0, 0.1, 0.4, 0.5};
    cylinder->location = {0, -0.3, 0.5, -0.2};
    cylinder->scale = cylinder->scale * 0.3;
    cylinder->left_transform = (quaternion){1, 0.3, 0.3, 0.1};
    cylinder->right_transform = 1.0 / cylinder->left_transform;

    shared_ptr<HyperTorus> torus = make_shared<HyperTorus>();
    torus->minor_radius = 0.5;
    torus->pigment_a = {0, 0.8, 0.4, 0.2};
    torus->pigment_b = {0, 0.1, 0.1, 0.2};
    torus->location = {0, 0.9, 0.5, -0.8};
    torus->scale = {1, 0.4, 0.3, 0.2};
    torus->left_transform = (quaternion){1, 0.3, 0.8, 0.1};
    torus->right_transform = 1.0 / torus->left_transform;

    vector<shared_ptr<RayTraceable>> objects;
    objects.push_back(torus);
    objects.push_back(sphere);
    objects.push_back(another_sphere);
    objects.push_back(cylinder);

    int width = 500;
    int height = 500;

    default_random_engine generator;
    generator.seed(random_device()());
    normal_distribution<real> distribution(0.0, 0.2 / (real) width);

    cout << "P3" << endl;
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
            pixel = clip_color(pixel);
            cout << setw(4) << left << (int) (pixel.x * 255);
            cout << setw(4) << left << (int) (pixel.y * 255);
            cout << setw(4) << left << (int) (pixel.z * 255);
            cout << " ";
        }
        cout << endl;
    }
    return EXIT_SUCCESS;
}
