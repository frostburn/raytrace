#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <vector>
#include <memory>
#include <random>
#include <cassert>

#include "raytrace/quaternion.h"

using namespace std;

const real ACCURACY = 0.01;
const int RAY_MARCH_ITERATIONS = 5000;
const int BISECT_ITERATIONS = 100;
const int NEWTON_ITERATIONS = 10;
const real NEWTON_EPSILON = 1e-5;

const real EPSILON = 1e-12;

class Raytraceable {
public:
    quaternion scale {1, 1, 1, 1};
    quaternion left_transform {1, 0, 0, 0};
    quaternion right_transform {1, 0, 0, 0};
    quaternion location {0, 0, 0, 0};
    bool reflective {false};
    virtual real s_distance(quaternion) = 0;
    virtual quaternion normal(quaternion) = 0;
    virtual color get_color(quaternion, quaternion) = 0;
    virtual tuple<quaternion, quaternion> trace(quaternion, quaternion);
    virtual tuple<quaternion, quaternion, real> trace_S3(quaternion, quaternion);
    quaternion project(quaternion);
    quaternion inverse_project(quaternion);
};

quaternion Raytraceable::project(quaternion location) {
    location = (quaternion){
        location.r * this->scale.r,
        location.x * this->scale.x,
        location.y * this->scale.y,
        location.z * this->scale.z,
    };
    return this->left_transform * location * this->right_transform;
}

quaternion Raytraceable::inverse_project(quaternion location) {
    location = (1.0 / this->left_transform) * location / this->right_transform;
    return (quaternion){
        location.r / this->scale.r,
        location.x / this->scale.x,
        location.y / this->scale.y,
        location.z / this->scale.z,
    };
}

tuple<quaternion, quaternion> Raytraceable::trace(quaternion source, quaternion target) {
    quaternion direction = target - source;
    direction = direction / norm(direction);
    target = source + direction * ACCURACY;

    source = this->inverse_project(source - this->location);
    target = this->inverse_project(target - this->location);

    real a = this->s_distance(source);
    real b = this->s_distance(target);
    // Ray march
    for (int i = 0; i < RAY_MARCH_ITERATIONS; ++i) {
        if (a*b <= 0) {
            break;
        }
        quaternion temp = target;
        target = 2.0 * target - source;
        source = temp;
        a = b;
        b = this->s_distance(target);
    }
    if (a*b > 0) {
        return {infinity(), direction};
    }
    // Ray bisect
    for (int i = 0; i < BISECT_ITERATIONS; ++i) {
        quaternion half_way = 0.5 * (source + target);
        real c = this->s_distance(half_way);
        if (a*c >= 0) {
            a = c;
            source = half_way;
        } else if (b*c >= 0) {
            b = c;
            target = half_way;
        } else {
            cerr << "Bijection collapse" << endl;
            break;
        }
    }

    return {this->project(0.5 * (source + target)) + this->location, direction};
}

tuple<quaternion, quaternion, real> Raytraceable::trace_S3(quaternion source, quaternion target) {
    target = cross_align(source, target);
    // Ray geodesic: source * cos(t) + target * sin(t)

    quaternion sp = this->inverse_project(source);
    quaternion tp = this->inverse_project(target);
    quaternion c = this->inverse_project(-this->location);

    // Ray march
    quaternion v = sp + c;
    real a = this->s_distance(v);
    real b;
    real t = 0;
    real dt = 2 * M_PI / (real) RAY_MARCH_ITERATIONS;
    for (int i = 0; i < RAY_MARCH_ITERATIONS; ++i) {
        // cout << a << "@" << t << endl;
        t += dt;
        v = sp*cos(t) + tp*sin(t) +  c;
        b = this->s_distance(v);
        if (a*b <= 0) {
            break;
        }
        a = b;
    }
    if (a*b > 0) {
        return {infinity(), infinity(), numeric_limits<real>::infinity()};
    }
    // Ray bisect
    for (int i = 0; i < BISECT_ITERATIONS; ++i) {
        dt *= 0.5;
        real half_way = t + dt;
        v = sp*cos(t) + tp*sin(t) +  c;
        real c = this->s_distance(v);
        if (a*c >= 0) {
            a = c;
            t = half_way;
        } else if (b*c >= 0) {
            b = c;
            t = half_way;
        } else {
            cerr << "Bijection collapse" << endl;
            break;
        }
    }

    return {source*cos(t) + target*sin(t), target*cos(t) - source*sin(t), t};
}

class GradientTraceable: public Raytraceable {
public:
    virtual quaternion gradient(quaternion) = 0;
    virtual tuple<quaternion, quaternion, real> trace_S3(quaternion, quaternion);
};

class Sphere: public GradientTraceable {
public:
    real s_distance(quaternion);
    quaternion normal(quaternion);
    quaternion gradient(quaternion);
    color get_color(quaternion, quaternion);
    tuple<quaternion, quaternion> trace(quaternion, quaternion);
};

real Sphere::s_distance(quaternion location) {
    return norm(location) - 1;
}

quaternion Sphere::normal(quaternion location) {
    return location / norm(location);
}

quaternion Sphere::gradient(quaternion location) {
    return this->normal(location);
}

color Sphere::get_color(quaternion location, quaternion ray) {
    quaternion n = this->normal(location);
    real angle = -dot(n, ray);
    if (angle < 0.1) {
        angle = 0.1 - 0.1 * fabs(angle);
    }
    return (color) {0, fabs(angle), 0, 0};
}

tuple<quaternion, quaternion> Sphere::trace(quaternion source, quaternion target) {
    // Solve for t
    // |v|^2 = 1
    // v = source + t * (target - source)
    // ~ |u + t * w|^2 = 1
    // (u + t * w)*(u + t * w)' = 1
    // |u|^2 + t * (u*w' + w*u') + t^2*|w|^2 = 1
    // |w|^2*t^2 + 2*dot(u, w) * t + |u|^2 - 1 = 0
    // t = (-dot(u, w) +- sqrt(dot(u, w)^2+|w|^2*(1-|u|^2))) / |w|^2

    quaternion direction = target - source;
    direction = direction / norm(direction);

    source = this->inverse_project(source - this->location);
    target = this->inverse_project(target - this->location);

    quaternion u = source;
    quaternion w = target - source;
    w = w / norm(w);
    real b = dot(u, w);
    real d = b*b + 1.0 - dot(u, u);
    if (d < 0) {
        return {infinity(), direction};
    }
    d = sqrt(d);
    real t = -b + d;
    real t1 = -b - d;
    if (t1 < t && t1 > 0) {
        t = t1;
    }
    if (t < 0) {
        return {infinity(), direction};
    }
    return {this->project(u + t * w) + this->location, direction};
}

class CliffordTorus: public GradientTraceable {
public:
    real s_distance(quaternion);
    quaternion normal(quaternion);
    quaternion gradient(quaternion);
    color get_color(quaternion, quaternion);
};

real CliffordTorus::s_distance(quaternion location) {
    real rx = sqrt(location.r * location.r + location.x * location.x);
    real yz = sqrt(location.y * location.y + location.z * location.z);
    if (rx > yz) {
        return rx - 1;
    }
    return yz - 1;
}

quaternion CliffordTorus::normal(quaternion location) {
    real rx = sqrt(location.r * location.r + location.x * location.x);
    real yz = sqrt(location.y * location.y + location.z * location.z);
    if (rx > yz) {
        return (quaternion){location.r / rx, location.x / rx, 0, 0};
    }
    return (quaternion){0, 0, location.y / yz, location.z / yz};
}

quaternion CliffordTorus::gradient(quaternion location) {
    return this->normal(location);
}

color CliffordTorus::get_color(quaternion location, quaternion ray) {
    quaternion n = this->normal(location);
    real angle = -dot(n, ray);
    if (angle < 0.1) {
        angle = 0.1 - 0.1 * fabs(angle);
    }
    real pattern = cos(10*location.r)*cos(10*location.x)*cos(10*location.y)*cos(10*location.z);
    return (color) {0, (0.5 + 0.4 * pattern) * angle, 0, 0};
}

tuple<quaternion, quaternion, real> GradientTraceable::trace_S3(quaternion source, quaternion target) {
    target = cross_align(source, target);
    // Ray geodesic: source * cos(t) + target * sin(t)

    quaternion sp = this->inverse_project(source);
    quaternion tp = this->inverse_project(target);
    quaternion c = this->inverse_project(-this->location);

    // Signed distance estimate: s(source*cos(t) + target*sin(t))
    // ds / dt = dot(-source*sin(t) + target*cos(t), grad s)

    quaternion v;
    real t = 0;
    real cos_t = 1;
    real sin_t = 0;

    // Newton's method
    for (int i = 0; i < NEWTON_ITERATIONS; ++i) {
        cos_t = cos(t);
        sin_t = sin(t);
        // v = inverse_project(source*cost_t + target*sin_t - this->location)
        v = sp*cos_t + tp*sin_t + c;
        quaternion grad = this->gradient(v);
        real dt = dot(tp*cos_t - sp*sin_t, grad);
        if (fabs(dt) < EPSILON) {
            dt = ACCURACY;
        } else {
            real f = this->s_distance(v);
            dt = -f / dt;
        }
        if (dt > 0.5 * M_PI) {
            dt = 0.5 * M_PI;
        } else if (dt < -0.5 * M_PI) {
            dt = -0.5 * M_PI;
        }
        if (t + dt < 0) {
            dt = ACCURACY;
        } else if (t + dt > 2 * M_PI) {
            dt = -ACCURACY;
        }
        t += dt;
    }
    if (this->s_distance(v) > NEWTON_EPSILON) {
        return {infinity(), infinity(), numeric_limits<real>::infinity()};
    }

    return {source*cos_t + target*sin_t, target*cos_t - source*sin_t, t};
}

color raytrace(quaternion source, quaternion target, int depth, const vector<shared_ptr<Raytraceable>>& objects) {
    if (depth <= 0) {
        return (color){0, 0.05, 0, 0};
    }
    vector<shared_ptr<Raytraceable>>::const_iterator obj;
    color result = {0, 0.05, 0, 0};
    quaternion closest_surface = infinity();

    for (obj = objects.begin(); obj != objects.end(); ++obj) {
        auto [surface, ray] = (*obj)->trace(source, target);
        if (norm(surface - source) >= norm(closest_surface - source)) {
            continue;
        }
        closest_surface = surface;
        surface = (*obj)->inverse_project(surface - (*obj)->location);
        ray = (*obj)->inverse_project(ray);
        result = (*obj)->get_color(surface, ray / norm(ray));
        if ((*obj)->reflective) {
            quaternion normal = (*obj)->normal(surface);
            ray = ray - 2 * normal * dot(ray, normal);
            // surface = (*obj)->project(surface) + (*obj)->location;
            ray = (*obj)->project(ray);
            result = 0.1 * result + 0.9 * raytrace(
                closest_surface + ACCURACY * ray,
                surface + ray,
                depth - 1,
                objects
            );
        }
    }
    return result;
}

color raytrace_S3(quaternion source, quaternion target, int depth, const vector<shared_ptr<Raytraceable>>& objects) {
    if (depth <= 0) {
        return (color){0, 0.05, 0, 0};
    }
    vector<shared_ptr<Raytraceable>>::const_iterator obj;
    color result = {0, 0.05, 0, 0};
    quaternion closest_surface = infinity();
    real closest_distance = numeric_limits<real>::infinity();

    for (obj = objects.begin(); obj != objects.end(); ++obj) {
        auto [surface, ray, distance] = (*obj)->trace_S3(source, target);
        if (distance >= closest_distance) {
            continue;
        }
        closest_surface = surface;
        closest_distance = distance;
        surface = (*obj)->inverse_project(surface - (*obj)->location);
        ray = (*obj)->inverse_project(ray);
        result = (*obj)->get_color(surface, ray / norm(ray));
        if ((*obj)->reflective) {
            quaternion normal = (*obj)->normal(surface);
            ray = ray - 2 * normal * dot(ray, normal);
            // surface = (*obj)->project(surface) + (*obj)->location;
            ray = (*obj)->project(ray);
            source = closest_surface + ACCURACY * ray;
            target = source + ray;
            source = source / norm(source);
            target = target / norm(target);
            result = 0.1 * result + 0.9 * raytrace_S3(
                source,
                target,
                depth - 1,
                objects
            );
        }
    }
    return result;
}

int main_cartesian ()
{
    quaternion camera_pos = {0, 1, 0.5, -2.2};
    quaternion look_at = {0, 0, 0.3, 0};
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
    sphere->location = {0, 0, -0.2, 0};
    sphere->scale = sphere->scale * 1.4;
    sphere->scale.y = 0.3;
    sphere->reflective = true;

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

    vector<shared_ptr<Raytraceable>> objects;
    objects.push_back(sphere);
    objects.push_back(another_sphere);
    objects.push_back(cylinder);

    int width = 100;
    int height = 100;

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

int main_S3()
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

    vector<shared_ptr<Raytraceable>> objects;
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
                // quaternion q = {1, 0, 0, 0.2};
                // q = q / norm(q);
                // point->location = q * point->location;
            }
        }
    }
    int width = 2048;
    int height = 2048;

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

void test_trace() {
    quaternion source = {0, 0, 0, -2};
    quaternion target = {0, 0, 0, 0};
    Sphere *sphere = new Sphere;
    auto [surface, direction] = sphere->trace(source, target);
    assert(norm(surface - (quaternion){0, 0, 0, -1}) < ACCURACY);
}

void test_scale() {
    quaternion source = {0, 0, 3, 0};
    quaternion target = {0, 0, 0, 0};
    Sphere *sphere = new Sphere;
    sphere->scale = sphere->scale * 0.5;
    auto [surface, direction] = sphere->trace(source, target);
    assert(norm(surface - (quaternion){0, 0, 0.5, 0}) < ACCURACY);
}

void test_location() {
    quaternion source = {0, 1, 0.2, 0};
    quaternion target = {0, 0, 0, 0.1};
    Sphere *sphere = new Sphere;
    sphere->location = {0, 0, 0.2, 0.1};
    sphere->scale = sphere->scale * 0.5;
    auto [surface, direction] = sphere->trace(source, target);
    assert(norm(surface - (quaternion){0, 0.486928, 0.0973857, 0.0513072}) < ACCURACY);
    real radius = norm(surface - sphere->location);
    assert(radius > 0.5 - ACCURACY);
    assert(radius < 0.5 + ACCURACY);
}

void test_trace_S3() {
    default_random_engine generator;
    generator.seed(random_device()());
    normal_distribution<real> distribution(0.0, 1);

    quaternion source = {-1, 0, 0, 0};
    quaternion target = {0, 1, 0, 0};
    Sphere *sphere = new Sphere;
    sphere->location = (quaternion){1, 0, 0, 0};
    auto [surface, direction, distance] = sphere->trace_S3(source, target);
    real radius = norm(surface - sphere->location);
    assert(radius < 1 + ACCURACY);
    assert(radius > 1 - ACCURACY);
    assert(distance < M_PI);
    assert(direction.r > 0);
    assert(direction.x < 0);
}

void test_trace_S3_random() {
    default_random_engine generator;
    generator.seed(random_device()());
    normal_distribution<real> distribution(0.0, 1);

    Sphere *sphere = new Sphere;
    quaternion scale = sphere->scale;
    for (int i = 0; i < 100; ++i) {
        quaternion source = {
            distribution(generator),
            distribution(generator),
            distribution(generator),
            distribution(generator)
        };
        quaternion target = {
            distribution(generator),
            distribution(generator),
            distribution(generator),
            distribution(generator)
        };
        sphere->location = {
            distribution(generator),
            distribution(generator),
            distribution(generator),
            distribution(generator)
        };
        source = source / norm(source);
        target = target / norm(target);
        sphere->scale = scale * fabs(distribution(generator));

        auto [surface, direction, distance] = sphere->trace_S3(source, target);
        if (norm(target - sphere->location) < sphere->scale.r - ACCURACY) {
            quaternion v = sphere->inverse_project(target) + sphere->inverse_project(-sphere->location);
            assert(norm(v) < 1 + ACCURACY);

            assert(distance <= 2*M_PI);
            real radius = norm(surface - sphere->location);
            assert(radius < sphere->scale.r + ACCURACY);
            assert(radius > sphere->scale.r - ACCURACY);
        }
    }
};

int main() {
    test_trace();
    test_scale();
    test_location();
    test_trace_S3();
    // test_trace_S3_random();  // TODO: Fix
    // return main_cartesian();
    return main_S3();
}
