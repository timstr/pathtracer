#include <Geometry.hpp>
#include <RandomNumberGenerator.hpp>

namespace {
    constexpr float pi = 3.141592653589793f;
    constexpr float epsilon = 1e-6f;
}

Triangle::Triangle(Pos _a, Pos _b, Pos _c) noexcept
    : a(_a)
    , b(_b)
    , c(_c) {

}

Vec Triangle::normal() const noexcept {
    return ((b - a) ^ (c - a)).unit();
}

float Triangle::area() const noexcept {
    return ((b - a) ^ (c - a)).norm() * 0.5f;
}

Sphere::Sphere(Pos _center, float _radius) noexcept
    : center(_center)
    , radius(_radius) {

}

float Sphere::surfaceArea() const noexcept {
    return 4.0f * pi * radius * radius;
}

float Sphere::volume() const noexcept {
    return (4.0f / 3.0f) * pi * radius * radius * radius;
}

Vec Sphere::normal(Pos p) const noexcept {
    return (p - center).unit();
}

std::optional<float> intersect(const Ray& ray, const Triangle& triangle) noexcept {
    // https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
    const auto edge1 = triangle.b - triangle.a;
    const auto edge2 = triangle.c - triangle.b;
    const auto h = ray.dir ^ edge2;
    const auto a = edge1 * h;
    if (std::abs(a) < epsilon) {
        return std::nullopt;
    }
    const auto f = 1.0f / a;
    const auto s = ray.pos - triangle.a;
    const auto u = f * (s * h);
    if (u < 0.0f || u > 1.0f) {
        return std::nullopt;
    }
    const auto q = s ^ edge1;
    const auto v = f * (ray.dir * q);
    if (v < 0.0f || u + v > 1.0f) {
        return std::nullopt;
    }
    const auto t = f * (edge2 * q);
    if (t > epsilon) {
        return t;
    }
    return std::nullopt;
}

std::optional<float> intersect(const Ray& ray, const Sphere& sphere) noexcept {
    // https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
    const auto a = ray.dir.normSquared();
    const auto b = 2.0f * (ray.dir * (ray.pos - sphere.center));
    const auto c = (ray.pos - sphere.center).normSquared() - (sphere.radius * sphere.radius);
    const auto disc = b * b - 4.0f * a * c;
    if (disc <= 0.0) {
        return std::nullopt;
    }
    const auto sqrtDisc = std::sqrt(disc);
    const auto t0 = (-b - sqrtDisc) / (2.0f * a);
    const auto t1 = (-b + sqrtDisc) / (2.0f * a);
    const auto epsilon = 1e-3f;
    if (t0 > epsilon) {
        if (t1 > epsilon) {
            return std::min(t0, t1);
        }
        return t0;
    } else {
        if (t1 > epsilon) {
            return t1;
        }
        return std::nullopt;
    }
}

Vec randomHemisphereVectorUniform(Vec normal) noexcept {
    const auto dist = std::uniform_real_distribution<float>{-1.0f, 1.0f};
    while (true) {
        auto x = dist(randomEngine());
        auto y = dist(randomEngine());
        auto z = dist(randomEngine()) * 0.5f + 0.5f;
        const auto l = std::sqrt(x * x + y * y + z * z);
        if (l > 1.0f || l < epsilon) {
            continue;
        }
        x /= l;
        y /= l;
        z /= l;

        const auto normalParallelToX = std::abs(1.0f - normal * Vec(1.0f, 0.0f, 0.0f)) < epsilon;
        Vec v1 = normal ^ (normalParallelToX ? Vec(0.0f, 1.0f, 0.0f) : Vec(1.0f, 0.0f, 0.0f));
        Vec v2 = v1 ^ normal;
        return x * v1 + y * v2 + z * normal;
    }
}

Ray::Ray(Vec _direction, Pos _position) noexcept
    : dir(_direction)
    , pos(_position) {
}
