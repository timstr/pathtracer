#pragma once

// TODO: rename this to Colour because we're in Canada
class Color {
public:
    explicit Color(float _r = 0.0f, float _g = 0.0f, float _b = 0.0f) noexcept;

    Color& operator+=(const Color& other) noexcept;

    float r;
    float g;
    float b;
};

Color operator+(const Color& l, const Color& r) noexcept;

Color operator*(float k, const Color& c) noexcept;
