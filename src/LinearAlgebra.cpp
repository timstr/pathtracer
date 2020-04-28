#include <LinearAlgebra.hpp>

#include <cmath>
#include <cassert>

namespace {
    constexpr float epsilon = 1e-6f;
}

Vec::Vec(float _x, float _y, float _z) noexcept
    : x(_x)
    , y(_y)
    , z(_z) {
    
}

float Vec::norm() const noexcept {
    return std::sqrt(normSquared());
}

float Vec::normSquared() const noexcept {
    return x * x + y * y + z * z;
}

Vec Vec::unit() const noexcept {
    const auto l = norm();
    assert(l >= epsilon);
    return Vec(x / l, y / l, z / l);
}

Pos::Pos(float _x, float _y, float _z) noexcept
    : x(_x)
    , y(_y)
    , z(_z) {

}

Linear Linear::Identity() noexcept {
    return Linear(
        1.0f, 0.0f, 0.0f,
        0.0f, 1.0f, 0.0f,
        0.0f, 0.0f, 1.0f
    );
}

Linear Linear::RotationX(float angle) noexcept {
    const auto c = std::cos(angle);
    const auto s = std::sin(angle);
    return Linear(
        1.0f, 0.0f, 0.0f,
        0.0f, c,    -s,
        0.0f, s,    c
    );
}

Linear Linear::RotationY(float angle) noexcept {
    const auto c = std::cos(angle);
    const auto s = std::sin(angle);
    return Linear(
        c,    0.0f, s,
        0.0f, 1.0f, 0.0f,
        -s,   0.0f, c
    );
}

Linear Linear::RotationZ(float angle) noexcept {
    const auto c = std::cos(angle);
    const auto s = std::sin(angle);
    return Linear(
        c,    -s,   0.0f,
        s,    c,    0.0f,
        0.0f, 0.0f, 1.0f
    );
}

Linear Linear::Rotation(Vec u, float a) noexcept {
    // https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
    const auto c = std::cos(a);
    const auto s = std::sin(a);
    const auto omc = 1.0f - c;
    return Linear(
        c + u.x * u.x * omc,        u.x * u.y * omc - u.z * s,  u.x * u.z * omc + u.y * s,
        u.y * u.x * omc + u.z * s,  c + u.y * u.y * omc,        u.y * u.z * omc - u.x * s,
        u.z * u.x * omc - u.y * s,  u.z * u.y * omc + u.x * s,  c + u.z * u.z * omc
    );
}

Linear::Linear() noexcept
    : m_data({
        1.0f, 0.0f, 0.0f,
        0.0f, 1.0f, 0.0f,
        0.0f, 0.0f, 1.0f
    }) {

}

Linear::Linear(Vec v1, Vec v2, Vec v3) noexcept
    : m_data({
        v1.x, v2.x, v3.x,
        v1.y, v2.y, v3.y,
        v1.z, v2.z, v3.z
    }) {

}

Linear::Linear(float a, float b, float c, float d, float e, float f, float g, float h, float i) noexcept 
    : m_data({
        a, b, c,
        d, e, f,
        g, h, i
    }) {
}

Vec Linear::row(int i) const noexcept {
    assert(i <= 3);
    const auto& self = *this;
    return Vec(self(i, 0), self(i, 1), self(i, 2));
}

Vec Linear::column(int i) const noexcept {
    assert(i <= 3);
    const auto& self = *this;
    return Vec(self(0, i), self(1, i), self(2, i));
}

float Linear::operator[](int idx) const noexcept {
    assert(idx <= 9);
    return m_data[idx];
}

float& Linear::operator[](int idx) noexcept {
    assert(idx <= 9);
    return m_data[idx];
}

float Linear::operator()(int i, int j) const noexcept {
    assert(i <= 3);
    assert(j <= 3);
    return m_data[3 * j + i];
}

float& Linear::operator()(int i, int j) noexcept {
    assert(i <= 3);
    assert(j <= 3);
    return m_data[3 * j + i];
}

Linear Linear::transpose() const noexcept {
    const auto& [a, b, c, d, e, f, g, h, i] = m_data;
    return Linear(a, d, g, b, e, h, c, f, i);
}

float Linear::determinant() const noexcept {
    const auto& [a, b, c, d, e, f, g, h, i] = m_data;
    return (a * e * i) + (b * f * g) + (c * d * h)
        - (c * e * g) - (b * d * i) - (a * f * h);
}

std::optional<Linear> Linear::inverse() const noexcept {
    // https://en.wikipedia.org/wiki/Invertible_matrix#Inversion_of_3_%C3%97_3_matrices

    const auto det = determinant();

    if (det < epsilon) {
        return std::nullopt;
    }

    const auto& [a, b, c, d, e, f, g, h, i] = m_data;

    const auto A = e * i - f * h;
    const auto B = f * g - d * i;
    const auto C = d * h - e * g;
    const auto D = c * h - b * i;
    const auto E = a * i - c * g;
    const auto F = b * g - a * h;
    const auto G = b * f - c * e;
    const auto H = c * d - a * f;
    const auto I = a * e - b * d;

    return Linear(A, D, G, B, E, H, C, F, I) / det;
}

Affine::Affine(Linear _linear, Vec _translation) noexcept 
    : linear(_linear)
    , translation(_translation) {
}

Vec operator-(const Vec& v) noexcept {
    return Vec(-v.x, -v.y, -v.z);
}

Vec operator*(const Vec& v, float t) noexcept {
    return Vec(v.x * t, v.y * t, v.z * t);
}

Vec operator*(float t, const Vec& v) noexcept {
    return v * t;
}

Vec operator/(const Vec& v, float t) noexcept {
    assert(std::abs(t) > epsilon);
    return v * (1.0f / t);
}

Vec operator+(const Vec& u, const Vec& v) noexcept {
    return Vec(u.x + v.x, u.y + v.y, u.z + v.z);
}

Vec operator-(const Vec& u, const Vec& v) noexcept {
    return Vec(u.x - v.x, u.y - v.y, u.z - v.z);
}

Pos operator+(const Pos& p, const Vec& v) noexcept {
    return Pos(p.x + v.x, p.y + v.y, p.z + v.z);
}

Pos operator+(const Vec& v, const Pos& p) noexcept {
    return p + v;
}

Pos operator-(const Pos& p, const Vec& v) noexcept {
    return Pos(p.x - v.x, p.y - v.y, p.z - v.z);
}

Vec operator-(const Pos& p, const Pos& q) noexcept {
    return Vec(p.x - q.x, p.y - q.y, p.z - q.z);
}

float operator*(const Vec& u, const Vec& v) noexcept {
    return u.x * v.x + u.y * v.y + u.z * v.z;
}

Vec operator^(const Vec& u, const Vec& v) noexcept {
    return Vec(u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x);
}

Linear operator-(const Linear& A) noexcept {
    return Linear(
        -A[0],
        -A[1],
        -A[2],
        -A[3],
        -A[4],
        -A[5],
        -A[6],
        -A[7],
        -A[8]
    );
}

Linear operator+(const Linear& A, const Linear& B) noexcept {
    return Linear(
        A[0] + B[0],
        A[1] + B[1],
        A[2] + B[2],
        A[3] + B[3],
        A[4] + B[4],
        A[5] + B[5],
        A[6] + B[6],
        A[7] + B[7],
        A[8] + B[8]
    );
}

Linear operator-(const Linear& A, const Linear& B) noexcept {
    return Linear(
        A[0] - B[0],
        A[1] - B[1],
        A[2] - B[2],
        A[3] - B[3],
        A[4] - B[4],
        A[5] - B[5],
        A[6] - B[6],
        A[7] - B[7],
        A[8] - B[8]
    );
}

Linear operator*(const Linear& A, float t) noexcept {
    return Linear(
        t * A[0],
        t * A[1],
        t * A[2],
        t * A[3],
        t * A[4],
        t * A[5],
        t * A[6],
        t * A[7],
        t * A[8]
    );
}

Linear operator*(float t, const Linear& A) noexcept {
    return A * t;
}

Linear operator/(const Linear& A, float t) noexcept {
    assert(std::abs(t) > epsilon);
    return A * (1.0f / t);
}

Vec operator*(const Linear& A, const Vec& v) noexcept {
    return A.column(0) * v.x + A.column(1) * v.y + A.column(2) * v.z;
}

Pos operator*(const Linear& A, const Pos& p) noexcept {
    const auto v = A * Vec(p.x, p.y, p.z);
    return Pos(v.x, v.y, v.z);
}

Linear operator*(const Linear& A, const Linear& B) noexcept {
    return Linear(
        A.row(0) * B.column(0), A.row(0) * B.column(1), A.row(0) * B.column(2),
        A.row(1) * B.column(0), A.row(1) * B.column(1), A.row(1) * B.column(2),
        A.row(2) * B.column(0), A.row(2) * B.column(1), A.row(2) * B.column(2)
    );
}

Vec operator*(const Affine& T, const Vec& v) {
    return T.linear * v;
}

Pos operator*(const Affine& T, const Pos& p) {
    return T.linear * p + T.translation;
}

Linear Linear::ScaleX(float k) {
    return Linear(k, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f);
}

Linear Linear::ScaleY(float k) {
    return Linear(1.0f, 0.0f, 0.0f, 0.0f, k, 0.0f, 0.0f, 0.0f, 1.0f);
}

Linear Linear::ScaleZ(float k) {
    return Linear(1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, k);
}

Linear Linear::Scale(float k) {
    return Linear(k, 0.0f, 0.0f, 0.0f, k, 0.0f, 0.0f, 0.0f, k);
}
