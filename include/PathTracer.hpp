#pragma once

#include <Geometry.hpp>
#include <ObjectTree.hpp>

#include <memory>
#include <optional>
#include <vector>


class Color {
public:
    explicit Color(float _r = 0.0f, float _g = 0.0f, float _b = 0.0f) noexcept;

    float r;
    float g;
    float b;
};

// TODO: rename this
// TODO: add:
// - current index of refraction
// - distance from last bounce
// TODO: replace color with scalar, add wavelength
class ColorBounce {
public:
    ColorBounce(Color _emitted, Color _attenuation, Vec _direction) noexcept;

    Color emitted;
    Color attenuation;
    Vec direction;
};

// TODO: add an affine transformation to each Object,
// cache its inverse and normal transformation, and
// use these to transform rays during collision detection

class Object {
public:
    virtual ~Object() noexcept = default;

    // If the object hits the ray, returns the t-value along the ray
    // at which the object hits. The point of collision is the ray's
    // position plus t times the ray's direction. This function need
    // not be deterministic, and may return randomly distributed values
    // to allow effects like subsurface scattering, volumetric rendering,
    // and other creative noisy effects
    virtual std::optional<float> hit(const Ray& ray) const noexcept = 0;

    // Given a ray whose position is given by a previous call to hit (see above),
    // returns:
    // - Ray   : The newly deflected light ray, which may be reflected or refracted,
    //           located at the point of collision and pointing in the new ray direction
    // - Color : The emitted radiance at the point of collision in the direction of the
    //           incoming ray
    // - float : The attenuation factor of the radiance from the next ray, computed according
    //           to the bidirectional reflection distribution function
    // This should be general enough to account for glossy surfaces, geometric primitives,
    // non-trivial geometric objects, refractive objects, partially transparent objects with
    // sub-surface scattering, volumetric objects like smoke, and weird light-deflecting media
    virtual ColorBounce deflect(const Ray& ray) const noexcept = 0;

    // Returns the (ideally smallest) axis-aligned rectangle that fully contains the object
    // This will be recomputed for every render
    virtual Box getBoundingBox() const noexcept = 0;
};

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

class TriangleObject : public Object {
public:
    Triangle geometry;
    BasicMaterial material;

    TriangleObject(Triangle _geometry, BasicMaterial _material = {});

    std::optional<float> hit(const Ray& ray) const noexcept override;

    ColorBounce deflect(const Ray& ray) const noexcept override;

    Box getBoundingBox() const noexcept override;
};

class SphereObject : public Object {
public:
    Sphere geometry;
    BasicMaterial material;
    
    SphereObject(Sphere _geometry, BasicMaterial _material = {}) noexcept;

    std::optional<float> hit(const Ray& ray) const noexcept override;

    ColorBounce deflect(const Ray& ray) const noexcept override;

    Box getBoundingBox() const noexcept override;
};

class BoxObject : public Object {
public:
    Box geometry;
    BasicMaterial material;

    BoxObject(Box _geometry, BasicMaterial _material = {}) noexcept;

    std::optional<float> hit(const Ray& ray) const noexcept override;

    ColorBounce deflect(const Ray& ray) const noexcept override;

    Box getBoundingBox() const noexcept override;
};

class FractalObject : public Object {
public:
    BasicMaterial material;
    
    std::optional<float> hit(const Ray& ray) const noexcept override;

    ColorBounce deflect(const Ray& ray) const noexcept override;

    Box getBoundingBox() const noexcept override;

private:
    float signedDistance(const Pos&) const noexcept;
};

class Scene {
public:
    Scene();

    void addObject(std::unique_ptr<Object> object) noexcept;

    Color trace(Ray ray, std::size_t depth) const noexcept;

    void updateGeometry();

private:
    std::vector<std::unique_ptr<Object>> m_objects;

    ObjectTree::Tree m_objectTree;
};

class Camera {
public:
    virtual ~Camera() noexcept = default;

    virtual Ray getViewRay(float screenX, float screenY) const noexcept = 0;
};

class PerspectiveCamera : public Camera {
public:
    PerspectiveCamera(Affine transform, float aspectRatio = 1.0f, float fieldOfView = 30.0f) noexcept;

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

    Ray getViewRay(float screenX, float screenY) const noexcept override final;

private:
    Affine m_transform;
    float m_aspectRatio;
    float m_fieldOfView;
    float m_focalDistance;
    float m_focalBlurRadius;
};

class Image {
public:
    Image(std::size_t width, std::size_t height);

    std::size_t width() const noexcept;
    std::size_t height() const noexcept;

    const Color& operator()(std::size_t x, std::size_t y) const noexcept;
    Color& operator()(std::size_t x, std::size_t y) noexcept;

private:
    const std::size_t m_width;
    const std::size_t m_height;
    std::vector<Color> m_data;
};

class Renderer {
public:
    Renderer(std::size_t width = 1024, std::size_t height = 1024) noexcept;

    std::size_t width() const noexcept;
    std::size_t height() const noexcept;
    std::size_t numBounces() const noexcept;
    std::size_t samplesPerPixel() const noexcept;
    std::size_t numThreads() const noexcept;

    void setWidth(std::size_t) noexcept;
    void setHeight(std::size_t) noexcept;
    void setNumBounces(std::size_t) noexcept;
    void setSamplesPerPixel(std::size_t) noexcept;
    void setNumThreads(std::size_t) noexcept;

    Image render(const Scene&, const Camera&) const;

private:
    std::size_t m_width;
    std::size_t m_height;
    std::size_t m_numBounces;
    std::size_t m_samplesPerPixel;
    std::size_t m_numThreads;
};

class ToneMapper {
public:
    virtual ~ToneMapper() noexcept = default;

    virtual Image operator()(const Image&) const noexcept = 0;
};

class ReinhardToneMapper : public ToneMapper {
public:
    virtual Image operator()(const Image&) const noexcept override final;
};