#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <vector>
#include <memory>
#include <random>
#include <cassert>
#include <string>

#include "raytrace/quaternion.h"
#include "raytrace/ray_traceable.h"
#include "raytrace/gradient_traceable.h"
#include "raytrace/sphere.h"
#include "raytrace/clifford_torus.h"
#include "raytrace/hyper_torus.h"
#include "raytrace/trace.h"
#include "raytrace/input_parser.h"

using namespace std;

int main(int argc, char *argv[])
{
    InputParser input(argc, argv);
    if(input.cmdOptionExists("-h") || input.cmdOptionExists("--help")){
        cerr << "Usage: " << argv[0] << " --time 3.14 --width 200 --height 200" << endl;
        return EXIT_SUCCESS;
    }

    real t = 0;
    int width = 100;
    int height = 100;

    const std::string &time = input.getCmdOption("--time");
    if (!time.empty()){
        t = std::stod(time);
    }
    const std::string &width_ = input.getCmdOption("--width");
    if (!width_.empty()){
        width = std::stoi(width_);
    }
    const std::string &height_ = input.getCmdOption("--height");
    if (!height_.empty()){
        height = std::stoi(height_);
    }
    real aspect_ratio = (real)width / (real)height;

    quaternion camera_transform = {1, 0.2, 0.05, -0.5};
    camera_transform = camera_transform / norm(camera_transform);
    quaternion target_transform = {1, 0.05, -0.15, 0.05};
    target_transform = target_transform / norm(target_transform);

    // TODO: Sensible camera definitions
    // quaternion camera_pos = {1, 0.1, 0.05, -0.3};
    // camera_pos = camera_pos / norm(camera_pos);

    // quaternion look_at = {1, 0, 0, 0};
    // quaternion look_direction = look_at - camera_pos;
    // look_direction = look_direction / norm(look_direction);

    real view_width = 0.2;
    real view_depth = 0.5;

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
                point->scale = point->scale * 0.02;
                objects.push_back(point);
                quaternion q = {0, 0, 0, 1};
                q = exp(q * t);
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

    shared_ptr<HyperTorus2> axis = make_shared<HyperTorus2>();
    axis->location = {0, 0, 0, 0};
    axis->minor_radius = 0.01;
    quaternion q = {0, 1, 0, 0};
    axis->left_transform = exp(q * t * 0.25);
    axis->right_transform = exp(q * t * 0.25);
    objects.push_back(axis);

    axis = make_shared<HyperTorus2>();
    axis->location = {0, 0, 0, 0};
    axis->minor_radius = 0.01;
    quaternion p = {sqrt(0.5), sqrt(0.5), 0, 0};
    q = {0, 0, 1, 0};
    axis->left_transform = p * exp(q * t * 0.25);
    axis->right_transform = exp(q * t * 0.25) / p;
    objects.push_back(axis);

    default_random_engine generator;
    normal_distribution<real> distribution(0.0, 0.2 / (real) width);

    cout << "P3" << endl;
    cout << width << " " << height << endl;
    cout << 255 << endl;
    for (int j = 0; j < height; ++j) {
        cerr << (j * 100) / height << "%" << endl;
        real y_ = 1 - 2 * (j / (real) height);
        for (int i = 0; i < width; ++i) {
            real x = 1 - 2 * (i / (real) width);
            x *= aspect_ratio;
            x += distribution(generator);
            real y = y_ + distribution(generator);
            quaternion target = {0, x*view_width, y*view_width, view_depth};
            target = 1 + target_transform * target / target_transform;
            target = target / norm(target);
            color pixel = raytrace_S3(camera_transform, target * camera_transform, 6, objects);
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
