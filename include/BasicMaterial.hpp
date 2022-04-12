#pragma once

#include <Color.hpp>
#include <ColorBounce.hpp>
#include <LinearAlgebra.hpp>

class BasicMaterial {
public:
    BasicMaterial() noexcept;

    float diffuseReflection() const noexcept;
    float specularReflection() const noexcept;
    float specularSharpness() const noexcept;
    Color reflectedAbsorption() const noexcept;
    Color emittedLuminance() const noexcept;
    float transmittance() const noexcept;
    float indexOfRefraction() const noexcept;
    Color internalAbsorption() const noexcept;

    void setDiffuseReflection(float) noexcept;
    void setSpecularReflection(float) noexcept;
    void setSpecularSharpness(float) noexcept;
    void setReflectedAbsorption(Color) noexcept;
    void setEmittedLuminance(Color) noexcept;
    void setTransmittance(float) noexcept;
    void setIndexOfRefraction(float) noexcept;
    void setInternalAbsorption(Color) noexcept;

    ColorBounce deflect(const Vec& inbound, const Vec& normal) const noexcept;

private:
    float m_diffuseReflection;
    float m_specularReflection;
    float m_specularSharpness;
    Color m_reflectedAbsorption;
    Color m_emittedLuminance;
    float m_transmittance;
    float m_indexOfRefraction;
    Color m_internalAbsorption;
    // TODO:
    // float m_internalDensity;
    // float m_internalScatterSharpness;
};
