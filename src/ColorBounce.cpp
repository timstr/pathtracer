#include <ColorBounce.hpp>

#include <cassert>
#include <cmath>

ColorBounce::ColorBounce(
    Color emitted,
    Color attenuation,
    Vec rayDirection,
    Vec normal
) noexcept
    : emitted(emitted)
    , attenuation(attenuation)
    , rayDirection(rayDirection)
    , normal(normal) {
    assert(std::abs(1.0f - rayDirection.norm()) < 1e-5);
    assert(std::abs(1.0f - normal.norm()) < 1e-5);
    assert(attenuation.r >= 0.0f && attenuation.r <= 1.0f);
    assert(attenuation.g >= 0.0f && attenuation.g <= 1.0f);
    assert(attenuation.b >= 0.0f && attenuation.b <= 1.0f);
}
