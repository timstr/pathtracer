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

Box::Box(Pos corner, Pos oppositeCorner) noexcept
    : center(corner + 0.5f * (oppositeCorner - corner))
    , halfSize(Vec{
        0.5f * std::abs(corner.x - oppositeCorner.x),
        0.5f * std::abs(corner.y - oppositeCorner.y),
        0.5f * std::abs(corner.z - oppositeCorner.z)
    }){
}

Box::Box(Pos _center, Vec _halfSize) noexcept
    : center(_center)
    , halfSize(_halfSize) {

}

float Box::surfaceArea() const noexcept {
    return 8.0f * (
        halfSize.x * halfSize.y +
        halfSize.x * halfSize.z +
        halfSize.y * halfSize.z
    );
}

float Box::volume() const noexcept {
    return halfSize.x * halfSize.y * halfSize.z * 8.0f;
}

Box boxContaining(const Box& a, const Box& b) noexcept {
    const auto aBegin = a.center - a.halfSize;
    const auto aEnd = a.center + a.halfSize;
    const auto bBegin = b.center - b.halfSize;
    const auto bEnd = b.center + b.halfSize;
    const auto p0 = Pos{
        std::min(aBegin.x, bBegin.x),
        std::min(aBegin.y, bBegin.y),
        std::min(aBegin.z, bBegin.z),
    };
    const auto p1 = Pos{
        std::max(aEnd.x, bEnd.x),
        std::max(aEnd.y, bEnd.y),
        std::max(aEnd.z, bEnd.z),
    };
    return Box{p0, p1};
}

Vec bounce(const Vec& inbound, const Vec& normal) noexcept {
    return inbound - (2.0f * ((inbound * normal) * normal));
}

std::optional<float> intersect(const Ray& ray, const Triangle& triangle) noexcept {
    // https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
    const auto edge1 = triangle.b - triangle.a;
    const auto edge2 = triangle.c - triangle.a;
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
    const auto eps= 1e-3f;
    if (t0 > eps) {
        if (t1 > eps) {
            return std::min(t0, t1);
        }
        return t0;
    } else {
        if (t1 > eps) {
            return t1;
        }
        return std::nullopt;
    }
}

std::optional<float> intersect(const Ray& ray, const Box& box) noexcept {
    auto closestT = std::optional<float>{};
    
    const auto recordClosest = [&](float t) {
        closestT = closestT.has_value() ? std::min(*closestT, t) : t;
    };

    // negative-facing corner of box, relative to ray origin
    const auto bBegin = box.center - box.halfSize - ray.pos;

    // positive-facing corner of box, relative to ray origin
    const auto bEnd = box.center + box.halfSize - ray.pos;

    using Dimension = const float Vec::*;

    const auto otherDimWithinBounds = [&](Dimension dim, float t) {
        const auto v = ray.dir.*dim * t;
        return v >= bBegin.*dim && v <= bEnd.*dim;
    };

    const auto projectAlongDim = [&](Dimension dim, Vec corner) {
        if (std::abs(ray.dir.*dim) < epsilon) {
            return;
        }
        const auto t = corner.*dim / ray.dir.*dim;

        if (
            t > 0.0f &&
            (dim == &Vec::x || otherDimWithinBounds(&Vec::x, t)) &&
            (dim == &Vec::y || otherDimWithinBounds(&Vec::y, t)) &&
            (dim == &Vec::z || otherDimWithinBounds(&Vec::z, t))
        ) {
            recordClosest(t);
        }
    };

    projectAlongDim(&Vec::x, bBegin);
    projectAlongDim(&Vec::x, bEnd);
    projectAlongDim(&Vec::y, bBegin);
    projectAlongDim(&Vec::y, bEnd);
    projectAlongDim(&Vec::z, bBegin);
    projectAlongDim(&Vec::z, bEnd);

    return closestT;
}

bool inside(const Pos& p, const Sphere& s) noexcept {
    return (s.center - p).normSquared() <= (s.radius * s.radius);
}

bool inside(const Pos& p, const Box& b) noexcept {
    const auto q = p - b.center;
    return
        (std::abs(q.x) < b.halfSize.x) &&
        (std::abs(q.y) < b.halfSize.y) &&
        (std::abs(q.z) < b.halfSize.z);
}

std::pair<float, float> randomPointInCircle() noexcept {
    const auto dist = std::uniform_real_distribution<float>{-1.0f, 1.0f};
    while (true) {
        auto x = dist(randomEngine());
        auto y = dist(randomEngine());
        const auto l = x * x + y * y;
        if (l <= 1.0f) {
            return {x, y};
        }
    }
}

Vec randomPointOnHemisphereUniform(Vec normal) noexcept {
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
        return (x * v1 + y * v2 + z * normal).unit();
    }
}

Ray::Ray(Vec _direction, Pos _position) noexcept
    : dir(_direction)
    , pos(_position) {
}
