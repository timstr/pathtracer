#include <Geometry.hpp>
#include <RandomNumberGenerator.hpp>

#include <cassert>

namespace {
    constexpr float pi = 3.141592653589793f;
    constexpr float epsilon = 1e-3f;
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

Sphere::Sphere(float _radius) noexcept
    : radius(_radius) {
    assert(radius > 0.0f);
}

float Sphere::surfaceArea() const noexcept {
    return 4.0f * pi * radius * radius;
}

float Sphere::volume() const noexcept {
    return (4.0f / 3.0f) * pi * radius * radius * radius;
}

float Sphere::signedDistance(Pos p) const noexcept {
    return p.toVec().norm() - radius;
}

Vec Sphere::normal(Pos p) const noexcept {
    return p.toVec().unit();
}

Rectangle::Rectangle(Vec _halfSize) noexcept
    : halfSize(_halfSize) {
    assert(halfSize.x >= epsilon);
    assert(halfSize.y >= epsilon);
    assert(halfSize.z >= epsilon);
}

float Rectangle::surfaceArea() const noexcept {
    return 8.0f * (
        halfSize.x * halfSize.y +
        halfSize.x * halfSize.z +
        halfSize.y * halfSize.z
    );
}

float Rectangle::volume() const noexcept {
    return halfSize.x * halfSize.y * halfSize.z * 8.0f;
}

float Rectangle::signedDistance(Pos p) const noexcept {
    const auto dx = std::abs(p.x) - this->halfSize.x;
    const auto dy = std::abs(p.y) - this->halfSize.y;
    const auto dz = std::abs(p.z) - this->halfSize.z;
    return Vec(
        std::max(dx, 0.0f),
        std::max(dy, 0.0f),
        std::max(dz, 0.0f)
    ).norm() + std::min(std::max(std::max(dx, dy), dz), 0.0f);
}

Vec Rectangle::normal(Pos p) const noexcept {
    const auto ax = std::abs(p.x / halfSize.x);
    const auto ay = std::abs(p.y / halfSize.y);
    const auto az = std::abs(p.z / halfSize.z);

    const auto sign = [](float v) {
        return v > 0.0f ? 1.0f : -1.0f;
    };

    if (ax > ay) {
        if (ax > az) {
            // ax > ay and ax > az
            return Vec{sign(p.x), 0.0f, 0.0f};
        } else {
            // az >= ax > ay
            return Vec{0.0f, 0.0f, sign(p.z)};
        }
    } else {
        if (ay > az) {
            // ay >= ax and ay > az
            return Vec{0.0f, sign(p.y), 0.0f};
        } else {
            // az >= ay >= ax
            return Vec{0.0f, 0.0f, sign(p.z)};
        }
    }
}

AxisAlignedBox boxContaining(const AxisAlignedBox& a, const AxisAlignedBox& b) noexcept {
    assert(a.halfSize.x > 0.0f);
    assert(a.halfSize.y > 0.0f);
    assert(a.halfSize.z > 0.0f);
    assert(b.halfSize.x > 0.0f);
    assert(b.halfSize.y > 0.0f);
    assert(b.halfSize.z > 0.0f);
    const auto aBegin = a.center - a.halfSize;
    const auto aEnd = a.center + a.halfSize;
    const auto bBegin = b.center - b.halfSize;
    const auto bEnd = b.center + b.halfSize;
    assert(aBegin.x < aEnd.x);
    assert(aBegin.y < aEnd.y);
    assert(aBegin.z < aEnd.z);
    assert(bBegin.x < bEnd.x);
    assert(bBegin.y < bEnd.y);
    assert(bBegin.z < bEnd.z);
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
    const auto bb = AxisAlignedBox{p0, p1};
    // assert(inside(aBegin, bb));
    // assert(inside(aEnd,   bb));
    // assert(inside(bBegin, bb));
    // assert(inside(bEnd,   bb));
    // assert(!inside(aBegin - 0.1f * a.halfSize, bb) || !inside(bBegin - 0.1f * b.halfSize, bb));
    // assert(!inside(aEnd   + 0.1f * a.halfSize, bb) || !inside(bEnd   + 0.1f * b.halfSize, bb));
    return bb;
}

Vec bounce(const Vec& inbound, const Vec& normal) noexcept {
    return inbound - (2.0f * ((inbound * normal) * normal));
}

std::optional<Pos> intersect(const Ray& ray, const Triangle& triangle) noexcept {
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
        return ray.pos + t * ray.dir;
    }
    return std::nullopt;
}

std::optional<Pos> intersect(const Ray& ray, const Sphere& sphere) noexcept {
    // https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
    const auto a = ray.dir.normSquared();
    const auto b = 2.0f * (ray.dir * ray.pos.toVec());
    const auto c = ray.pos.toVec().normSquared() - (sphere.radius * sphere.radius);
    const auto disc = b * b - 4.0f * a * c;
    if (disc <= 0.0) {
        return std::nullopt;
    }
    const auto sqrtDisc = std::sqrt(disc);
    const auto t0 = (-b - sqrtDisc) / (2.0f * a);
    const auto t1 = (-b + sqrtDisc) / (2.0f * a);
    const auto eps= epsilon;
    auto t = float{};
    if (t0 > eps) {
        if (t1 > eps) {
            t = std::min(t0, t1);
        } else {
            t = t0;
        }
    } else {
        if (t1 > eps) {
            t = t1;
        } else {
            return std::nullopt;
        }
    }
    return ray.pos + t * ray.dir;
}

std::optional<Pos> intersect(const Ray& ray, const Rectangle& box) noexcept {
    auto closestT = std::optional<float>{};

    const auto recordClosest = [&](float t) {
        closestT = closestT.has_value() ? std::min(*closestT, t) : t;
    };

    // negative-facing corner of box, relative to ray origin
    const auto bBegin = -ray.pos.toVec() - box.halfSize;

    // positive-facing corner of box, relative to ray origin
    const auto bEnd = -ray.pos.toVec() + box.halfSize;

    // ray is now at origin

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

    if (!closestT.has_value()) {
        return std::nullopt;
    }
    return ray.pos + ((*closestT) * ray.dir);
}

std::optional<Pos> intersect(const Ray& ray, const AxisAlignedBox& b) noexcept {
    auto r = Ray(ray.dir, ray.pos - b.center.toVec());
    auto p = intersect(r, Rectangle{b.halfSize});
    if (!p.has_value()) {
        return std::nullopt;
    }
    return (*p) + b.center.toVec();
}

bool inside(const Pos& p, const Sphere& s) noexcept {
    return p.toVec().normSquared() <= (s.radius * s.radius);
}

bool inside(const Pos& p, const Rectangle& b) noexcept {
    return
        (std::abs(p.x) <= b.halfSize.x) &&
        (std::abs(p.y) <= b.halfSize.y) &&
        (std::abs(p.z) <= b.halfSize.z);
}

bool inside(const Pos& p, const AxisAlignedBox& b) noexcept {
    return inside(p - b.center.toVec(), Rectangle{b.halfSize});
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

std::pair<Vec, Vec> orthogonalPair(Vec v) noexcept {
    v = v.unit();
    const auto xAxis = Vec(1.0f, 0.0f, 0.0f);
    const auto yAxis = Vec(0.0f, 1.0f, 0.0f);
    const auto alignedToX = ((v ^ xAxis).normSquared() < 0.1);
    const auto w = alignedToX ? yAxis : xAxis;
    const auto v1 = w ^ v;
    const auto v2 = v1 ^ v;
    return std::make_pair(v1.unit(), v2.unit()); // HACK: why are these vectors slightly not unit?
}

Vec randomPointOnHemisphereUniform(Vec normal) noexcept {
    normal = normal.unit();
    const auto dist = std::uniform_real_distribution<float>{-1.0f, 1.0f};
    while (true) {
        auto x = dist(randomEngine());
        auto y = dist(randomEngine());
        auto z = dist(randomEngine());
        const auto lSqr = (x * x) + (y * y) + (z * z);
        if (lSqr > 1.0f || lSqr < epsilon) {
            continue;
        }
        auto v = Vec(x, y, z);
        if ((v * normal) <= 0.0) {
            continue;
        }
        return v.unit();
    }
}

Vec randomPointOnHemisphereCosine(Vec normal) noexcept {
    normal = normal.unit();
    const auto dist = std::uniform_real_distribution<float>{-1.0f, 1.0f};
    while (true) {
        auto x = dist(randomEngine());
        auto y = dist(randomEngine());
        const auto lSqr = (x * x) + (y * y);
        if (lSqr > 1.0f) {
            continue;
        }
        const auto z = std::sqrt(1.0f - lSqr);
        const auto [v1, v2] = orthogonalPair(normal);
        assert(std::abs(normal * v1) < epsilon);
        assert(std::abs(normal * v2) < epsilon);
        assert(std::abs(v1 * v2) < epsilon);
        assert(std::abs(normal.norm() - 1.0f) < epsilon);
        assert(std::abs(v1.norm() - 1.0f) < epsilon);
        assert(std::abs(v2.norm() - 1.0f) < epsilon);
        assert(std::abs(Vec(x, y, z).norm() - 1.0f) < epsilon);
        return (z * normal) + (x * v1) + (y * v2);
    }
}

Ray::Ray(Vec _direction, Pos _position) noexcept
    : dir(_direction)
    , pos(_position) {
}

AxisAlignedBox::AxisAlignedBox(Pos corner, Pos oppositeCorner) noexcept
    : center(corner + 0.5f * (oppositeCorner - corner))
    , halfSize(Vec{
        0.5f * std::abs(corner.x - oppositeCorner.x),
        0.5f * std::abs(corner.y - oppositeCorner.y),
        0.5f * std::abs(corner.z - oppositeCorner.z)
    }) {
    assert(halfSize.x >= epsilon);
    assert(halfSize.y >= epsilon);
    assert(halfSize.z >= epsilon);
}

AxisAlignedBox::AxisAlignedBox(Pos _center, Vec _halfSize) noexcept
    : center(_center)
    , halfSize(_halfSize) {
    assert(halfSize.x >= epsilon);
    assert(halfSize.y >= epsilon);
    assert(halfSize.z >= epsilon);
}
