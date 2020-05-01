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

class Box {
public:
    Pos center;
    Vec halfSize;

    Box(Pos corner, Pos oppositeCorner) noexcept;

    Box(Pos _center, Vec _halfSize) noexcept;

    float surfaceArea() const noexcept;

    float volume() const noexcept;

    Vec normal(Pos p) const noexcept;
};

class Ray {
public:
    Ray() noexcept = default;
    Ray(Vec _direction, Pos _position) noexcept;

    Vec dir;
    Pos pos;
};

Box boxContaining(const Box& a, const Box& b) noexcept;

Vec bounce(const Vec& inbound, const Vec& normal) noexcept;

std::optional<float> intersect(const Ray&, const Triangle&) noexcept;

std::optional<float> intersect(const Ray&, const Sphere&) noexcept;

std::optional<float> intersect(const Ray&, const Box&) noexcept;

bool inside(const Pos&, const Sphere&) noexcept;

bool inside(const Pos&, const Box&) noexcept;

std::pair<float, float> randomPointInCircle() noexcept;

Vec randomPointOnHemisphereUniform(Vec normal) noexcept;