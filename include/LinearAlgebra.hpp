#pragma once

#include <array>
#include <cstdint>
#include <optional>

class Vec {
public:
    float x;
    float y;
    float z;

    Vec(float _x = 0.0f, float _y = 0.0f, float _z = 0.0f) noexcept;

    float norm() const noexcept;

    float normSquared() const noexcept;

    Vec unit() const noexcept;
};

class Pos {
public:
    float x;
    float y;
    float z;

    Pos(float _x = 0.0f, float _y = 0.0f, float _z = 0.0f) noexcept;
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

    Linear() noexcept;

    Linear(Vec column1, Vec column2, Vec column3) noexcept;

    Linear(float a, float b, float c, float d, float e, float f, float g, float h, float i) noexcept;

    Vec row(int i) const noexcept;
    Vec column(int i) const noexcept;

    float operator[](int idx) const noexcept;
    float& operator[](int idx) noexcept;

    float operator()(int i, int j) const noexcept;
    float& operator()(int i, int j) noexcept;

    Linear transpose() const noexcept;

    float determinant() const noexcept;

    std::optional<Linear> inverse() const noexcept;

private:
    std::array<float, 9> m_data;
};

class Affine {
public:
    Affine(Linear _linear = {}, Vec _translation = {}) noexcept;

    Linear linear;
    Vec translation;
};

// vector negation
Vec operator-(const Vec& v) noexcept;

// vector-scalar multiplication
Vec operator*(const Vec& v, float t) noexcept;
Vec operator*(float t, const Vec& v) noexcept;

// vector-scalar division
Vec operator/(const Vec& v, float t) noexcept;

// vector-vector addition and subtraction
Vec operator+(const Vec& u, const Vec& v) noexcept;
Vec operator-(const Vec& u, const Vec& v) noexcept;

// vector-point addition
Pos operator+(const Pos& p, const Vec& v) noexcept;
Pos operator+(const Vec& v, const Pos& p) noexcept;

// vector-point subtraction
Pos operator-(const Pos& p, const Vec& v) noexcept;

// position-position subtraction
Vec operator-(const Pos& p, const Pos& q) noexcept;

// vector dot product
float operator*(const Vec& u, const Vec& v) noexcept;

// vector cross product
Vec operator^(const Vec& u, const Vec& v) noexcept;

// linear negation
Linear operator-(const Linear& A) noexcept;

// linear-linear addition and subtraction
Linear operator+(const Linear& A, const Linear& B) noexcept;
Linear operator-(const Linear& A, const Linear& B) noexcept;

// linear-scalar multiplication
Linear operator*(const Linear& A, float t) noexcept;
Linear operator*(float t, const Linear& A) noexcept;

// linear-scalar division
Linear operator/(const Linear& A, float t) noexcept;

// linear-vector multiplication
Vec operator*(const Linear& A, const Vec& v) noexcept;

// linear-point multiplication
Pos operator*(const Linear& A, const Pos& v) noexcept;

// linear-linear multiplication
Linear operator*(const Linear& A, const Linear& B) noexcept;

// affine-vector multiplication
Vec operator*(const Affine& T, const Vec& v);

// affine-point multiplication
Pos operator*(const Affine& T, const Pos& p);