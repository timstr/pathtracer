#pragma once

#include <LinearAlgebra.hpp>


class Triangle {
public:
    Pos a;
    Pos b;
    Pos c;

    Triangle(Pos _a, Pos _b, Pos _c) noexcept;

    Vec normal() const noexcept;

    float area() const noexcept;
};

class Sphere {
public:
    Pos center;
    float radius;

    Sphere(Pos _center, float _radius) noexcept;

    float surfaceArea() const noexcept;

    float volume() const noexcept;

    Vec normal(Pos p) const noexcept;
};

class AxisAlignedBox {
public:
    Pos center;
    Vec halfSize;
};


class Ray {
public:
    Ray() noexcept = default;
    Ray(Vec _direction, Pos _position) noexcept;

    Vec dir;
    Pos pos;
};

std::optional<float> intersect(const Ray&, const Triangle&) noexcept;

std::optional<float> intersect(const Ray&, const Sphere&) noexcept;

Vec randomHemisphereVectorUniform(Vec normal) noexcept;