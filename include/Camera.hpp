#pragma once

#include <Geometry.hpp>
#include <LinearAlgebra.hpp>

class Camera {
public:
    Camera(Affine transform, float aspectRatio = 1.0f, float fieldOfView = 30.0f) noexcept;

    const Affine& transform() const noexcept;
    Affine& transform() noexcept;
    float aspectRatio() const noexcept;
    float fieldOfView() const noexcept;
    float focalDistance() const noexcept;
    float focalBlurRadius() const noexcept;

    void setTransform(Affine) noexcept;
    void setAspectRatio(float) noexcept;
    void setFieldOfView(float) noexcept;
    void setFocalDistance(float) noexcept;
    void setFocalBlurRadius(float) noexcept;

    Ray getViewRay(float screenX, float screenY) const noexcept;

private:
    Affine m_transform;
    float m_aspectRatio;
    float m_fieldOfView;
    float m_focalDistance;
    float m_focalBlurRadius;
};
