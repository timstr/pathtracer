#pragma once

#include <optional>

class Vec {
public:
    float x;
    float y;
    float z;

    constexpr Vec(float _x = 0.0f, float _y = 0.0f, float _z = 0.0f) noexcept;
};

class Pos {
public:
    float x;
    float y;
    float z;

    constexpr Pos(float _x = 0.0f, float _y = 0.0f, float _z = 0.0f) noexcept;
};

class Ray {
public:
    Vec dir;
    Pos pos;
};

constexpr Vec operator*(const Vec& v, float a) noexcept;
constexpr Vec operator*(float a, Vec& v) noexcept;

constexpr Vec operator+(const Vec& u, const Vec& v) noexcept;
constexpr Vec operator-(const Vec& u, const Vec& v) noexcept;

constexpr Pos operator+(const Pos& p, const Vec& v) noexcept;
constexpr Pos operator+(const Vec& v, const Pos& p) noexcept;

constexpr Pos operator-(const Pos& p, const Vec& v) noexcept;
constexpr Pos operator-(const Vec& v, const Pos& p) noexcept;

constexpr Vec operator-(const Pos& p, const Pos& q) noexcept;

constexpr float operator*(const Vec& u, const Vec& v) noexcept;

constexpr Vec operator^(const Vec& u, const Vec& v) noexcept;

class Triangle {
public:
    Pos a;
    Pos b;
    Pos c;
};

class Sphere {
public:
    Pos center;
    float radius;
};


class Color {
public:
    float r;
    float g;
    float b;
};

class Object {
public:
    virtual ~Object() noexcept = default;

    /**
     * If the given ray intersects the object, a new ray is returned whose
     * position is the nearest positive point of intersection, and whose
     * direction has been bounced according to the object's surface and material
     */
    virtual std::optional<Ray> bounceRay(const Ray& ray, Color& luminance) const noexcept;
};

class BasicGlossyMaterial {
public:
    float diffuseReflection;
    float specularReflection;
    Color color;
};

class TriangleObject : public Object {
    // TODO
};

class SphereObject : public Object {
    // TODO
};