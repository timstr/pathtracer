#pragma once

#include <array>
#include <cstddef>
#include <optional>

class Pos;

class Vec {
public:
    float x;
    float y;
    float z;

    explicit Vec(float _x = 0.0f, float _y = 0.0f, float _z = 0.0f) noexcept;

    float norm() const noexcept;

    float normSquared() const noexcept;

    Vec unit() const noexcept;

    Vec abs() const noexcept;


    Pos toPos() const noexcept;

    Vec& operator-=(const Vec& other) noexcept;

    Vec& operator+=(const Vec& other) noexcept;
};

class Pos {
public:
    float x;
    float y;
    float z;

    explicit Pos(float _x = 0.0f, float _y = 0.0f, float _z = 0.0f) noexcept;

    Pos& operator+=(const Vec& v) noexcept;
    Pos& operator-=(const Vec& v) noexcept;

    Vec toVec() const noexcept;
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

    Vec row(std::size_t i) const noexcept;
    Vec column(std::size_t i) const noexcept;

    float operator[](std::size_t idx) const noexcept;
    float& operator[](std::size_t idx) noexcept;

    float operator()(std::size_t i, std::size_t j) const noexcept;
    float& operator()(std::size_t i, std::size_t j) noexcept;

    Linear transpose() const noexcept;

    float determinant() const noexcept;

    std::optional<Linear> inverse() const noexcept;

    Linear& operator*=(const Linear& other) noexcept;

private:
    std::array<float, 9> m_data;
};

class Affine {
public:
    Affine() noexcept = default;
    explicit Affine(const Linear& _linear) noexcept;
    explicit Affine(const Linear& _linear, const Vec& _translation) noexcept;

    static Affine Translation(const Vec&) noexcept;
    static Affine Translation(float x, float y, float z) noexcept;

    std::optional<Affine> inverse() const noexcept;

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

// vector-pos addition
Pos operator+(const Pos& p, const Vec& v) noexcept;
Pos operator+(const Vec& v, const Pos& p) noexcept;

// vector-pos subtraction
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

// linear-pos multiplication
Pos operator*(const Linear& A, const Pos& v) noexcept;

// linear-linear multiplication
Linear operator*(const Linear& A, const Linear& B) noexcept;

// affine-vector multiplication
Vec operator*(const Affine& T, const Vec& v);

// affine-pos multiplication
Pos operator*(const Affine& T, const Pos& p);

// affine-linear multiplication
Affine operator*(const Affine& A, const Linear& L) noexcept;
Affine operator*(const Linear& L, const Affine& A) noexcept;

// affine-affine multiplication
Affine operator*(const Affine& A, const Affine& B) noexcept;
