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
    float radius;

    Sphere(float _radius) noexcept;

    float surfaceArea() const noexcept;

    float volume() const noexcept;

    Vec normal(Pos p) const noexcept;
};

class Rectangle {
public:
    Vec halfSize;

    Rectangle(Vec halfSize) noexcept;

    float surfaceArea() const noexcept;

    float volume() const noexcept;

    Vec normal(Pos p) const noexcept;
};

class AxisAlignedBox {
public:
    Pos center;
    Vec halfSize;

    AxisAlignedBox(Pos corner, Pos oppositeCorner) noexcept;

    AxisAlignedBox(Pos center, Vec halfSize) noexcept;
};

class Ray {
public:
    Ray() noexcept = default;
    Ray(Vec _direction, Pos _position) noexcept;

    Vec dir;
    Pos pos;
};

AxisAlignedBox boxContaining(const AxisAlignedBox& a, const AxisAlignedBox& b) noexcept;

Vec bounce(const Vec& inbound, const Vec& normal) noexcept;

std::optional<Pos> intersect(const Ray&, const Triangle&) noexcept;

std::optional<Pos> intersect(const Ray&, const Sphere&) noexcept;

std::optional<Pos> intersect(const Ray&, const Rectangle&) noexcept;

std::optional<Pos> intersect(const Ray&, const AxisAlignedBox&) noexcept;

bool inside(const Pos&, const Sphere&) noexcept;

bool inside(const Pos&, const Rectangle&) noexcept;

bool inside(const Pos&, const AxisAlignedBox&) noexcept;

std::pair<float, float> randomPointInCircle() noexcept;

Vec randomPointOnHemisphereUniform(Vec normal) noexcept;