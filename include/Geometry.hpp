#pragma once

#include <LinearAlgebra.hpp>


class Triangle {
public:
    Pos a;
    Pos b;
    Pos c;

    constexpr Triangle(Pos _a, Pos _b, Pos _c) noexcept;

    constexpr Vec normal() const noexcept;

    constexpr float area() const noexcept;
};

class Sphere {
public:
    Pos center;
    float radius;

    constexpr Sphere(Pos _center, float _radius) noexcept;

    constexpr float surfaceArea() const noexcept;

    constexpr float volume() const noexcept;
};

class AxisAlignedBox {
public:
    Pos center;
    Vec halfSize;
};


class Ray {
public:
    Vec dir;
    Pos pos;
};

std::optional<float> intersect(const Ray&, const Triangle&) noexcept;

std::optional<float> intersect(const Ray&, const Sphere&) noexcept;

Vec randomHemisphereVectorUniform(Vec normal) noexcept;