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
#include "raytrace/power3d.h"
#include "raytrace/ray_traceable.h"
#include "raytrace/gradient_traceable.h"
#include "raytrace/sphere.h"
#include "raytrace/mandelbulb.h"
#include "raytrace/clifford_torus.h"
#include "raytrace/hyper_torus.h"
#include "raytrace/trace.h"
#include "raytrace/input_parser.h"

using namespace std;

const int N = 1000;

int main(int argc, char *argv[])
{
    int width = 1000;
    int height = 1000;

    real view_width = 1.5;
    real view_depth = 0.1;

    real aspect_ratio = (real)width / (real)height;

    Mandelbulb *obj = new Mandelbulb;
    obj->threshold = 0.5;

    default_random_engine generator;
    normal_distribution<real> distribution(0.0, 0.2 / (real) width);

    quaternion loop[N];
    quaternion directions[N];
    for (int i = 0; i < N; ++i) {
        real t = 2*M_PI * i / (real) N;
        loop[i] = {0, view_width * cos(t), view_width * sin(t), view_depth};
        directions[i] = {0, -cos(t), -sin(t), 0};
    }

    for (int j = 0; j < 10000; ++j) {
        for (int i = 0; i < N; ++i) {
            real dist = obj->s_distance(loop[i]);
            loop[i] = loop[i] + 0.0001 * dist * directions[i];
        }
    }

    real min_dist = numeric_limits<real>::infinity();
    real max_dist = -min_dist;
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
            quaternion location = {0, x*view_width, y*view_width, view_depth};
            //color pixel = pow3d(location, 6);
            real dist = obj->s_distance(location);
            if (dist > max_dist) max_dist = dist;
            if (dist < min_dist) min_dist = dist;
            color pixel = {0, dist, sqrt(dist + 1), log(dist + 2)};
            pixel = pixel * 0.2;
            for (int k = 0; k < N; k++) {
                if (norm(location - loop[k]) < 0.01) {
                    pixel = {0, 0.6, 0.8, 0.9};
                }
            }
            pixel = clip_color(pixel);
            cout << setw(4) << left << (int) (pixel.x * 255);
            cout << setw(4) << left << (int) (pixel.y * 255);
            cout << setw(4) << left << (int) (pixel.z * 255);
            cout << " ";
        }
        cout << endl;
    }
    cerr << "Min:" << min_dist << endl;
    cerr << "Max:" << max_dist << endl;
    return EXIT_SUCCESS;
}
