#pragma once

#include <array>
#include <cstdint>
#include <optional>

class Vec {
public:
    float x;
    float y;
    float z;

    constexpr Vec(float _x = 0.0f, float _y = 0.0f, float _z = 0.0f) noexcept;

    constexpr float norm() const noexcept;

    constexpr float normSquared() const noexcept;

    constexpr Vec unit() const noexcept;
};

class Pos {
public:
    float x;
    float y;
    float z;

    constexpr Pos(float _x = 0.0f, float _y = 0.0f, float _z = 0.0f) noexcept;
};

class Linear {
public:
    static Linear Identity() noexcept;

    static Linear RotationX(float angle) noexcept;
    static Linear RotationY(float angle) noexcept;
    static Linear RotationZ(float angle) noexcept;

    static Linear Rotation(Vec axis, float angle) noexcept;

    static Linear ScaleX(float k);
    static Linear ScaleY(float k);
    static Linear ScaleZ(float k);

    static Linear Scale(float k);

    constexpr Linear() noexcept;

    constexpr Linear(Vec column1, Vec column2, Vec column3) noexcept;

    constexpr Linear(float a, float b, float c, float d, float e, float f, float g, float h, float i) noexcept;

    constexpr Vec row(int i) const noexcept;
    constexpr Vec column(int i) const noexcept;

    constexpr float operator[](int idx) const noexcept;
    constexpr float& operator[](int idx) noexcept;

    constexpr float operator()(int i, int j) const noexcept;
    constexpr float& operator()(int i, int j) noexcept;

    constexpr Linear transpose() const noexcept;

    constexpr float determinant() const noexcept;

    constexpr std::optional<Linear> inverse() const noexcept;

private:
    std::array<float, 9> m_data;
};

class Affine {
public:
    constexpr Affine(Linear _linear = {}, Vec _translation = {}) noexcept;

    Linear linear;
    Vec translation;
};

// vector negation
constexpr Vec operator-(const Vec& v) noexcept;

// vector-scalar multiplication
constexpr Vec operator*(const Vec& v, float t) noexcept;
constexpr Vec operator*(float t, Vec& v) noexcept;

// vector-scalar division
constexpr Vec operator/(const Vec& v, float t) noexcept;

// vector-vector addition and subtraction
constexpr Vec operator+(const Vec& u, const Vec& v) noexcept;
constexpr Vec operator-(const Vec& u, const Vec& v) noexcept;

// vector-point addition
constexpr Pos operator+(const Pos& p, const Vec& v) noexcept;
constexpr Pos operator+(const Vec& v, const Pos& p) noexcept;

// vector-point subtraction
constexpr Pos operator-(const Pos& p, const Vec& v) noexcept;

// position-position subtraction
constexpr Vec operator-(const Pos& p, const Pos& q) noexcept;

// vector dot product
constexpr float operator*(const Vec& u, const Vec& v) noexcept;

// vector cross product
constexpr Vec operator^(const Vec& u, const Vec& v) noexcept;

// linear negation
constexpr Linear operator-(const Linear& A) noexcept;

// linear-linear addition and subtraction
constexpr Linear operator+(const Linear& A, const Linear& B) noexcept;
constexpr Linear operator-(const Linear& A, const Linear& B) noexcept;

// linear-scalar multiplication
constexpr Linear operator*(const Linear& A, float t) noexcept;
constexpr Linear operator*(float t, const Linear& A) noexcept;

// linear-scalar division
constexpr Linear operator/(const Linear& A, float t) noexcept;

// linear-vector multiplication
constexpr Vec operator*(const Linear& A, const Vec& v) noexcept;

// linear-linear multiplication
constexpr Linear operator*(const Linear& A, const Linear& B) noexcept;

// affine-vector multiplication
constexpr Vec operator*(const Affine& T, const Vec& v);

// affine-point multiplication
constexpr Pos operator*(const Affine& T, const Pos& p);