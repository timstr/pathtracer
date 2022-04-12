#pragma once

#include <LinearAlgebra.hpp>
#include <Color.hpp>

// TODO: rename to ColourBounce
class ColorBounce {
public:
    ColorBounce(
        Color emitted,
        Color attenuation,
        Vec rayDirection,
        Vec normal
    ) noexcept;

    Color emitted;
    Color attenuation;
    Vec rayDirection;
    Vec normal;
};
