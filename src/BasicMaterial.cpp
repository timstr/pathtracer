#include <BasicMaterial.hpp>

#include <Geometry.hpp>
#include <RandomNumberGenerator.hpp>

#include <cassert>
#include <random>

BasicMaterial::BasicMaterial() noexcept
    : m_diffuseReflection(0.3f)
    , m_specularReflection(0.3f)
    , m_specularSharpness(0.9f)
    , m_reflectedAbsorption(Color{1.0f, 1.0f, 1.0f})
    , m_emittedLuminance(Color{0.0f, 0.0f, 0.0f})
    , m_transmittance(0.0f)
    , m_indexOfRefraction(1.5f)
    , m_internalAbsorption(Color{0.9f, 0.9f, 0.9f})
    {

}

float BasicMaterial::diffuseReflection() const noexcept {
    return m_diffuseReflection;
}

float BasicMaterial::specularReflection() const noexcept {
    return m_specularReflection;
}

float BasicMaterial::specularSharpness() const noexcept {
    return m_specularSharpness;
}

Color BasicMaterial::reflectedAbsorption() const noexcept {
    return m_reflectedAbsorption;
}

Color BasicMaterial::emittedLuminance() const noexcept {
    return m_emittedLuminance;
}

float BasicMaterial::transmittance() const noexcept {
    return m_transmittance;
}

float BasicMaterial::indexOfRefraction() const noexcept {
    return m_indexOfRefraction;
}

Color BasicMaterial::internalAbsorption() const noexcept {
    return m_internalAbsorption;
}

void BasicMaterial::setDiffuseReflection(float r) noexcept {
    assert(r >= 0.0f && r <= 1.0f);
    m_diffuseReflection = r;
}

void BasicMaterial::setSpecularReflection(float r) noexcept {
    assert(r >= 0.0f && r <= 1.0f);
    m_specularReflection = r;
}

void BasicMaterial::setSpecularSharpness(float r) noexcept {
    assert(r >= 0.0f && r <= 1.0f);
    m_specularSharpness = r;
}

void BasicMaterial::setReflectedAbsorption(Color c) noexcept {
    assert(c.r >= 0.0f && c.r <= 1.0f);
    assert(c.g >= 0.0f && c.g <= 1.0f);
    assert(c.b >= 0.0f && c.b <= 1.0f);
    m_reflectedAbsorption = c;
}

void BasicMaterial::setEmittedLuminance(Color c) noexcept {
    m_emittedLuminance = c;
}

void BasicMaterial::setTransmittance(float r) noexcept {
    assert(r >= 0.0f && r <= 1.0f);
    m_transmittance = r;
}

void BasicMaterial::setIndexOfRefraction(float i) noexcept {
    assert(i >= 1.0f);
    m_indexOfRefraction = i;
}

void BasicMaterial::setInternalAbsorption(Color c) noexcept {
    assert(c.r >= 0.0f && c.r <= 1.0f);
    assert(c.g >= 0.0f && c.g <= 1.0f);
    assert(c.b >= 0.0f && c.b <= 1.0f);
    m_internalAbsorption = c;
}

ColorBounce BasicMaterial::deflect(const Vec& inbound, const Vec& normal) const noexcept {
    const auto sinTheta = inbound * normal;

    if (sinTheta >= 0.0f) {
        // Collision from inside
        const auto v = (inbound + ((inbound * normal) * (1.0f - m_indexOfRefraction) * normal)).unit();
        if (v * normal >= 0.0f) {
            // Refraction to outside
            // TODO: diffuse/specular scatter
            return ColorBounce {
                Color{0.0f, 0.0f, 0.0f}, // TODO: emitted color?
                Color{1.0f, 1.0f, 1.0f}, // TODO: internal absorption (requires knowing distance traveled)
                v,
                normal
            };
        } else {
            // Total internal reflection
            return ColorBounce{
                Color{0.0f, 1.0f, 1.0f}, // HACK: yellow if total internal reflection
                Color{0.0f, 0.0f, 0.0f}, // TODO: internal absorption (requires knowing distance traveled)
                bounce(inbound, -normal),
                normal
            };
        }
    }

    const auto reflection = m_diffuseReflection + m_specularReflection;
    const auto options = reflection + m_transmittance;
    const auto dist = std::uniform_real_distribution<float>{0.0f, options};
    const auto which = dist(randomEngine());
    if (which < reflection){
        // Reflection
        if (which < m_diffuseReflection){
            // Diffuse reflection
            const auto v = randomPointOnHemisphereCosine(normal);
            return ColorBounce{
                m_emittedLuminance,
                m_reflectedAbsorption,
                v,
                normal
            };
        } else {
            // Specular reflection
            const auto s = (1.0f - m_specularSharpness) * randomPointOnHemisphereCosine(normal);
            const auto v = (bounce(inbound, normal) + s).unit();
            return ColorBounce{
                m_emittedLuminance,
                m_reflectedAbsorption,
                v,
                normal
            };
        }
    } else {
        // Transmittance
        // TODO: specular/diffuse scattering
        const auto v = (inbound + ((inbound * normal) * (1.0f - 1.0f / m_indexOfRefraction) * normal)).unit();
        return ColorBounce{
            m_emittedLuminance,
            m_reflectedAbsorption,
            v,
            normal
        };
    }


    //const auto sinTheta = inbound * normal;
    //const auto cosTheta = std::sqrt(1.0f - sinTheta * sinTheta);
    //const auto b = bounce(inbound, normal);



    /*
    assert(diffuseness >= 0.0f && diffuseness <= 1.0f);
    const auto n = (inbound * normal) >= 0.0f ? -normal : normal;
    const auto diffuseBounce = randomPointOnHemisphereUniform(n);
    const auto specularBounce = bounce(inbound, n);
    const auto ray = (specularBounce + diffuseness * (diffuseBounce - specularBounce)).unit();
    return {
        emittedRadiance,
        reflectedColor,
        ray
    };*/
}
