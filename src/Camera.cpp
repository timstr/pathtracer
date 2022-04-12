#include <Camera.hpp>

#include <RandomNumberGenerator.hpp>

#include <cassert>
#include <random>

Camera::Camera(Affine transform, float aspectRatio, float fieldOfView) noexcept
    : m_transform(transform)
    , m_aspectRatio(aspectRatio)
    , m_fieldOfView(fieldOfView)
    , m_focalDistance(10.0f)
    , m_focalBlurRadius(0.0f) {

}

const Affine& Camera::transform() const noexcept {
    return m_transform;
}

Affine& Camera::transform() noexcept {
    return m_transform;
}

float Camera::aspectRatio() const noexcept {
    return m_aspectRatio;
}

float Camera::fieldOfView() const noexcept {
    return m_fieldOfView;
}

float Camera::focalDistance() const noexcept {
    return m_focalDistance;
}

float Camera::focalBlurRadius() const noexcept {
    return m_focalBlurRadius;
}

void Camera::setTransform(Affine t) noexcept {
    m_transform = t;
}

void Camera::setAspectRatio(float r) noexcept {
    assert(r > 0.0f);
    m_aspectRatio = r;
}

void Camera::setFieldOfView(float fov) noexcept {
    assert(fov > -180.0f && fov < 180.0f);
    m_fieldOfView = fov;
}

void Camera::setFocalDistance(float d) noexcept {
    assert(d > 0.0f);
    m_focalDistance = d;
}

void Camera::setFocalBlurRadius(float r) noexcept {
    assert(r >= 0.0f);
    m_focalBlurRadius = r;
}

Ray Camera::getViewRay(float screenX, float screenY) const noexcept {
    const auto x = screenX * 2.0f - 1.0f;
    const auto y = screenY * 2.0f - 1.0f;
    const auto sp = m_aspectRatio > 1.0f ? Pos(x, y / m_aspectRatio) : Pos(x * m_aspectRatio, y);
    const auto dist = std::uniform_real_distribution<float>(0.0f, 1.0f);
    const auto blurAngle = 2.0f * 3.141592654f * dist(randomEngine());

    const auto [randX, randY] = randomPointInCircle();
    const auto blurRad = m_focalBlurRadius * std::max(m_aspectRatio, 1.0f / m_aspectRatio);
    const auto blurX = randX * blurRad;
    const auto blurY = randY * blurRad;
    const auto blurVec = Vec(blurX, blurY, 0.0f);
    const auto fovScale = std::tan(m_fieldOfView * 3.141592654f / 180.0f);
    const auto viewVec = Vec(fovScale * sp.x, fovScale * sp.y, 1.0f) + (blurVec / m_focalDistance);
    const auto dir = (transform() * viewVec).unit();
    const auto pos = transform() * (sp - blurVec);
    return { dir, pos };
}
