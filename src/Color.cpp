#include <Color.hpp>

Color::Color(float _r, float _g, float _b) noexcept
    : r(_r)
    , g(_g)
    , b(_b) {

}

Color& Color::operator+=(const Color& other) noexcept {
    this->r += other.r;
    this->g += other.g;
    this->b += other.b;
    return *this;
}

Color operator+(const Color& l, const Color& r) noexcept {
    return Color(
        l.r + r.r,
        l.g + r.g,
        l.b + r.b
    );
}

Color operator*(float k, const Color& c) noexcept {
    return Color(
        k * c.r,
        k * c.g,
        k * c.b
    );
}
