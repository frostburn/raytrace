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
    quaternion camera_pos = {1, 0.1, 0.05, -0.3};
    camera_pos = camera_pos / norm(camera_pos);

    vector<shared_ptr<RayTraceable>> objects;

    int grid_height = 3 - 1;
    int grid_width = 3 - 1;
    int grid_depth = 3 - 1;
    for (int j = 0; j <= grid_height; ++j) {
        real y = 1 - 2 * (j / (real) grid_height);
        for (int i = 0; i <= grid_width; ++i) {
            real x = 1 - 2 * (i / (real) grid_width);
            for (int k = 0; k <= grid_depth; ++k) {
                real z = 1 - 2 * (k / (real) grid_depth);
                shared_ptr<Sphere> point = make_shared<Sphere>();
                point->location = {1, 0.2 * x, 0.2 * y, 0.2 * z};
                point->location = point->location / norm(point->location);
                point->pigment = {0, 2.1 + x, 2.1 + y, 2.1 + z};
                point->pigment = point->pigment / 3.2;
                point->scale = point->scale * 0.04;
                objects.push_back(point);
                quaternion q = {1, 0, 0, 0.5};
                q = q / norm(q);
                point->location = q * point->location;
            }
        }
    }

    quaternion unit_scale = (quaternion){1, 1, 1, 1} * 0.05;
    shared_ptr<Sphere> unit = make_shared<Sphere>();
    unit->location = {1, 0, 0, 0};
    unit->scale = unit_scale;
    unit->reflection = 0.3;
    objects.push_back(unit);

    unit = make_shared<Sphere>();
    unit->location = {0, 1, 0, 0};
    unit->pigment = {0, 1, 0, 0};
    unit->scale = unit_scale;
    unit->reflection = 0.3;
    objects.push_back(unit);

    unit = make_shared<Sphere>();
    unit->location = {0, 0, 1, 0};
    unit->pigment = {0, 0, 1, 0};
    unit->scale = unit_scale;
    unit->reflection = 0.3;
    objects.push_back(unit);

    unit = make_shared<Sphere>();
    unit->location = {0, 0, 0, 1};
    unit->pigment = {0, 0, 0, 1};
    unit->scale = unit_scale;
    unit->reflection = 0.3;
    objects.push_back(unit);

    unit = make_shared<Sphere>();
    unit->location = {-1, 0, 0, 0};
    unit->pigment = {0, 0.2, 0.2, 0.2};
    unit->scale = unit_scale;
    unit->reflection = 0.5;
    objects.push_back(unit);


    int width = 500;
    int height = 500;

    default_random_engine generator;
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
            quaternion target = {1, 0.3 * x, 0.3 * y, 0.2};
            target = target / norm(target);
            color pixel = raytrace_S3(camera_pos, target, 6, objects);
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
