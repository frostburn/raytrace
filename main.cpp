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

int main()
{
    quaternion camera_pos = {1, 0.1, 0.05, 0};
    camera_pos = camera_pos / norm(camera_pos);

    // shared_ptr<Sphere> sphere = make_shared<Sphere>();
    // sphere->location = {1, 0, 0, 1};
    // sphere->location = sphere->location / norm(sphere->location);
    // sphere->scale = sphere->scale * 0.2;
    // sphere->reflective = true;

    // shared_ptr<Sphere> sphere2 = make_shared<Sphere>();
    // sphere2->location = {1, -0.03, -0.5, 0.7};
    // sphere2->location = sphere2->location / norm(sphere2->location);
    // sphere2->scale = sphere2->scale * 0.21;
    // sphere2->reflective = false;

    // shared_ptr<CliffordTorus> box = make_shared<CliffordTorus>();
    // box->location = {1, 0.82, -0.2, 0.5};
    // box->location = box->location / norm(box->location);
    // box->scale = (quaternion) {0.1, 0.05, 0.03, 0.03};
    // box->reflective = false;

    vector<shared_ptr<RayTraceable>> objects;
    // objects.push_back(sphere);
    // objects.push_back(sphere2);
    // objects.push_back(box);

    int grid_height = 4 - 1;
    int grid_width = 4 - 1;
    int grid_depth = 3;
    for (int j = 0; j <= grid_height; ++j) {
        real y = 1 - 2 * (j / (real) grid_height);
        for (int i = 0; i <= grid_width; ++i) {
            real x = 1 - 2 * (i / (real) grid_width);
            for (int k = 0; k < grid_depth; ++k) {
                real z = 1 + k;
                shared_ptr<Sphere> point = make_shared<Sphere>();
                point->location = {1, 0.2 * x, 0.2 * y, 0.2 * z};
                point->location = point->location / norm(point->location);
                point->scale = point->scale * 0.04;
                point->reflective = false;
                objects.push_back(point);
                quaternion q = {1, 0, 0, 0.2};
                q = q / norm(q);
                point->location = q * point->location;
            }
        }
    }
    int width = 500;
    int height = 500;

    default_random_engine generator;
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
            quaternion target = {1, 0.3 * x, 0.3 * y, 0.2};
            target = target / norm(target);
            color pixel = raytrace_S3(camera_pos, target, 6, objects);
            cout << setw(4) << left << (int) (pixel.x * 255);
        }
        cout << endl;
    }
    return EXIT_SUCCESS;
}
