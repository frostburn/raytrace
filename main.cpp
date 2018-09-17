#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <vector>
#include <memory>
#include <random>
using namespace std;

typedef double real;

const real ACCURACY = 0.005;
const int PATIENCE = 2000;
const int BISECT_ITERATIONS = 100;

struct quaternion
{
    real r;
    real x;
    real y;
    real z;
};

struct color
{
    real r;
    real g;
    real b;
};

std::ostream& operator<<(std::ostream& os, const quaternion& v)
{
    return os << "(quaternion){" <<  v.r << ", " << v.x << ", " << v.y << ", " << v.z << "}"; 
}

real norm(quaternion v) {
    return sqrt(v.x*v.x + v.y*v.y + v.z*v.z + v.r*v.r);
}

quaternion operator +(const quaternion& a, const quaternion& b) {
    return (quaternion){a.r + b.r, a.x + b.x, a.y + b.y, a.z + b.z};
}

quaternion operator -(const quaternion& a, const quaternion& b) {
    return (quaternion){a.r - b.r, a.x - b.x, a.y - b.y, a.z - b.z};
}

quaternion operator *(const quaternion& a, const quaternion& b) {
    return (quaternion){
        a.r*b.r - a.x*b.x - a.y*b.y - a.z*b.z,
        a.r*b.x + a.x*b.r + a.y*b.z - a.z*b.y,
        a.r*b.y - a.x*b.z + a.y*b.r + a.z*b.x,
        a.r*b.z + a.x*b.y - a.y*b.x + a.z*b.r
    };
}

quaternion operator *(const quaternion& a, const real b) {
    return (quaternion){a.r * b, a.x * b, a.y * b, a.z * b};
}

quaternion operator *(const real b, const quaternion& a) {
    return (quaternion){a.r * b, a.x * b, a.y * b, a.z * b};
}

quaternion operator /(const quaternion& a, const real b) {
    return a * (1.0 / b);
}

quaternion infinity () {
    real inf = numeric_limits<real>::infinity();
    return (quaternion){inf, inf, inf, inf};
}

real dot(const quaternion& a, const quaternion& b) {
    return a.r * b.r + a.x * b.x + a.y * b.y + a.z * b.z;
}

quaternion project(const quaternion& a, const quaternion& b) {
    return dot(a, b) * a;
}

class Raytraceable {
public:
    quaternion location;
    bool reflective;
    virtual real s_distance(quaternion) = 0;
    virtual quaternion normal(quaternion) = 0;
    virtual color get_color(quaternion, quaternion) = 0;
    virtual tuple<quaternion, quaternion> trace(quaternion, quaternion);
    tuple<quaternion, quaternion, real> trace_S3(quaternion, quaternion);
};

class Sphere: public Raytraceable {
public:
    real radius;
    real s_distance(quaternion);
    quaternion normal(quaternion);
    color get_color(quaternion, quaternion);
    tuple<quaternion, quaternion> trace(quaternion, quaternion);
};

real Sphere::s_distance(quaternion location) {
    return norm(this->location - location) - this->radius;
}

quaternion Sphere::normal(quaternion location) {
    quaternion n = location - this->location;
    return n / norm(n);
}

color Sphere::get_color(quaternion location, quaternion ray) {
    quaternion n = this->normal(location);
    real angle = -dot(n, ray);
    if (angle < 0.1) {
        angle = 0.1;
    }
    return (color) {angle, 0, 0};
}

class Box: public Raytraceable {
public:
    real length;
    real width;
    real height;
    real depth;
    real s_distance(quaternion);
    quaternion normal(quaternion);
    color get_color(quaternion, quaternion);
};

real Box::s_distance(quaternion location) {
    real temp;
    location = location - this->location;
    real distance = fabs(location.r) - 0.5 * this->length;
    temp = fabs(location.x) - 0.5 * this->width;
    if (temp > distance) {
        distance = temp;
    }
    temp = fabs(location.y) - 0.5 * this->height;
    if (temp > distance) {
        distance = temp;
    }
    temp = fabs(location.z) - 0.5 * this->depth;
    if (temp > distance) {
        distance = temp;
    }
    return distance;
}

quaternion Box::normal(quaternion location) {
    real temp;
    location = location - this->location;
    quaternion n = {location.r, 0, 0, 0};
    real distance = fabs(location.r) - 0.5 * this->length;
    temp = fabs(location.x) - 0.5 * this->width;
    if (temp > distance) {
        distance = temp;
        n = {0, location.x, 0, 0};
    }
    temp = fabs(location.y) - 0.5 * this->height;
    if (temp > distance) {
        distance = temp;
        n = {0, 0, location.y, 0};
    }
    temp = fabs(location.z) - 0.5 * this->depth;
    if (temp > distance) {
        distance = temp;
        n = {0, 0, 0, location.z};
    }
    return n / norm(n);
}

color Box::get_color(quaternion location, quaternion ray) {
    location = location * 5;
    real result = 1;
    // real result *= location.r - floor(location.r + 0.5);
    result *= location.x - floor(location.x + 0.5);
    result *= location.y - floor(location.y + 0.5);
    result *= location.z - floor(location.z + 0.5);
    return (color) {0.2 + 0.2 * (result > 0), 0, 0};
}

tuple<quaternion, quaternion> Raytraceable::trace(quaternion source, quaternion target) {
    quaternion direction = target - source;
    direction = direction / norm(direction);
    target = source + direction * ACCURACY;
    real a = this->s_distance(source);
    real b = this->s_distance(target);
    // Ray march
    for (int i = 0; i < PATIENCE; ++i) {
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
    return {0.5 * (source + target), direction};
}

tuple<quaternion, quaternion> Sphere::trace(quaternion source, quaternion target) {
    // Solve for t
    // |v - this->location|^2 = this->radius^2
    // v = source + t * (target - source)
    // ~ |u + t * w|^2 = r^2
    // (u + t * w)*(u + t * w)' = r^2
    // |u|^2 + t * (u*w' + w*u') + t^2*|w|^2 = r^2
    // |w|^2*t^2 + 2*dot(u, w) * t + |u|^2 - r^2 = 0
    // t = (-dot(u, w) +- sqrt(dot(u, w)^2+|w|^2*(r^2-|u|^2))) / |w|^2
    quaternion u = source - this->location;
    quaternion w = target - source;
    w = w / norm(w);
    real b = dot(u, w);
    real d = b*b + this->radius * this->radius - dot(u, u);
    if (d < 0) {
        return {infinity(), w};
    }
    d = sqrt(d);
    real t = -b + d;
    real t1 = -b - d;
    if (t1 < t && t1 > 0) {
        t = t1;
    }
    if (t < 0) {
        return {infinity(), w};
    }
    return {source + t * w, w};
}

tuple<quaternion, quaternion, real> Raytraceable::trace_S3(quaternion source, quaternion target) {
    quaternion direction = target - source;
    direction = direction / norm(direction);
    target = source + direction * ACCURACY;
    target = target / norm(target);
    real a = this->s_distance(source);
    real b = this->s_distance(target);
    real distance = 0;
    // Ray march
    for (int i = 0; i < PATIENCE; ++i) {
        if (a*b <= 0) {
            break;
        }
        quaternion temp = target;
        target = 2.0 * target - source;
        target = target / norm(target);
        source = temp;
        distance += norm(source - target);
        a = b;
        b = this->s_distance(target);
    }
    direction = target - source;
    direction = direction / norm(direction);
    if (a*b > 0) {
        return {infinity(), direction, numeric_limits<real>::infinity()};
    }
    // Ray bisect
    quaternion bisect_start = source;
    for (int i = 0; i < BISECT_ITERATIONS; ++i) {
        quaternion half_way = 0.5 * (source + target);
        half_way = half_way / norm(half_way);
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
    target = 0.5 * (source + target);
    distance += norm(bisect_start - target);
    // Not recalculating direction here due to floating point instability.
    return {target, direction, distance};
}

color raytrace(quaternion source, quaternion target, int depth, const vector<shared_ptr<Raytraceable>>& objects) {
    if (depth <= 0) {
        return (color){0, 0, 0};
    }
    vector<shared_ptr<Raytraceable>>::const_iterator obj;
    color result = {0, 0, 0};
    quaternion closest_surface = infinity();

    for (obj = objects.begin(); obj != objects.end(); ++obj) {
        auto [surface, ray] = (*obj)->trace(source, target);
        if (norm(surface - source) >= norm(closest_surface - source)) {
            continue;
        }
        closest_surface = surface;
        if ((*obj)->reflective) {
            quaternion normal = (*obj)->normal(surface);
            ray = ray - 2 * normal * dot(ray, normal);
            result = raytrace(surface + ACCURACY * ray, surface + ray, depth - 1, objects);
        } else {
            result = (*obj)->get_color(surface, ray / norm(ray));
        }
    }
    return result;
}

color raytrace_S3(quaternion source, quaternion target, int depth, const vector<shared_ptr<Raytraceable>>& objects) {
    if (depth <= 0) {
        return (color){0, 0, 0};
    }
    vector<shared_ptr<Raytraceable>>::const_iterator obj;
    color result = {0, 0, 0};
    quaternion closest_surface = infinity();
    real closest_distance = numeric_limits<real>::infinity();

    for (obj = objects.begin(); obj != objects.end(); ++obj) {
        auto [surface, ray, distance] = (*obj)->trace_S3(source, target);
        if (distance >= closest_distance) {
            continue;
        }
        closest_surface = surface;
        closest_distance = distance;
        if ((*obj)->reflective) {
            quaternion normal = (*obj)->normal(surface);
            ray = ray - 2 * normal * dot(ray, normal);
            result = raytrace_S3(surface + ACCURACY * ray, surface + ray, depth - 1, objects);
        } else {
            result = (*obj)->get_color(surface, ray / norm(ray));
        }
    }
    return result;
}

int main_cartesian ()
{
    quaternion camera_pos = {0, 1.232, 1, -3};
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
    sphere->location = {0, 0, 0.4, 0};
    sphere->radius = 0.4;
    sphere->reflective = true;

    shared_ptr<Sphere> another_sphere = make_shared<Sphere>();
    another_sphere->location = {0, 0.6, 0.3, 0};
    another_sphere->radius = 0.2;
    another_sphere->reflective = false;

    shared_ptr<Box> box = make_shared<Box>();
    box->location = {0, 0, -1, 0};
    box->length = 0;
    box->width = 1;
    box->height = 2;
    box->depth = 1;
    box->reflective = false;

    vector<shared_ptr<Raytraceable>> objects;
    objects.push_back(sphere);
    objects.push_back(another_sphere);
    objects.push_back(box);

    int width = 200;
    int height = 200;

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
            quaternion target = look_at + x * vleft + y * up;
            color pixel = raytrace(camera_pos, target, 10, objects);
            cout << setw(4) << left << (int) (pixel.r * 255);
        }
        cout << endl;
    }
    return EXIT_SUCCESS;
}

int main_S3()
{
    quaternion camera_pos = {1, 0, 0, 0};

    shared_ptr<Sphere> sphere = make_shared<Sphere>();
    sphere->location = {1, 0.02, 0.1, 0.5};
    sphere->location = sphere->location / norm(sphere->location);
    sphere->radius = 0.2;
    sphere->reflective = true;

    shared_ptr<Sphere> sphere2 = make_shared<Sphere>();
    sphere2->location = {1, -0.03, -0.5, 0.7};
    sphere2->location = sphere2->location / norm(sphere2->location);
    sphere2->radius = 0.1;
    sphere2->reflective = false;

    shared_ptr<Box> box = make_shared<Box>();
    box->location = {1, 0.82, -0.2, 0.5};
    box->location = box->location / norm(box->location);
    box->length = 0.2;
    box->width = 0.31;
    box->height = 0.31;
    box->depth = 0.23;
    box->reflective = false;

    vector<shared_ptr<Raytraceable>> objects;
    objects.push_back(sphere);
    objects.push_back(sphere2);
    objects.push_back(box);

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
            quaternion target = {1, x, y, 1};
            target = target / norm(target);
            color pixel = raytrace_S3(camera_pos, target, 6, objects);
            cout << setw(4) << left << (int) (pixel.r * 255);
        }
        cout << endl;
    }
    return EXIT_SUCCESS;
}

int main() {
    return main_cartesian();
}
