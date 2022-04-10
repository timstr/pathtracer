#pragma once

#include <Geometry.hpp>
#include <ObjectTree.hpp>

#include <barrier>
#include <memory>
#include <mutex>
#include <optional>
#include <vector>
#include <thread>

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

// TODO: rename this
// TODO: add:
// - current index of refraction
// - distance from last bounce
// TODO: replace color with scalar, add wavelength
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

// TODO: add an affine transformation to each Object,
// cache its inverse and normal transformation, and
// use these to transform rays during collision detection

class Object {
public:
    Object() noexcept = default;
    virtual ~Object() noexcept = default;

    const Affine& transformation() const noexcept;
    const Affine& inverseTransformation() const noexcept;
    void setTransformation(const Affine&) noexcept;

    // TODO: think about the normal matrix

    // If the object hits the ray, returns the t-value along the ray
    // at which the object hits. The point of collision is the ray's
    // position plus t times the ray's direction. This function need
    // not be deterministic, and may return randomly distributed values
    // to allow effects like subsurface scattering, volumetric rendering,
    // and other creative noisy effects
    std::optional<Pos> hitRay(const Ray& ray) const noexcept;

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
    ColorBounce deflectRay(const Ray& ray) const noexcept;

    AxisAlignedBox getBoundingBox() const noexcept;

protected:
    virtual std::optional<Pos> hitLocalRay(const Ray& ray) const noexcept = 0;

    virtual ColorBounce deflectLocalRay(const Ray& ray) const noexcept = 0;

    // Returns the (ideally smallest) axis-aligned rectangle that fully contains the object
    // This will be recomputed for every render
    virtual AxisAlignedBox getLocalBoundingBox() const noexcept = 0;

private:
    Affine m_transformation;
    Affine m_inverseTransformation;
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

    std::optional<Pos> hitLocalRay(const Ray& ray) const noexcept override;

    ColorBounce deflectLocalRay(const Ray& ray) const noexcept override;

    AxisAlignedBox getLocalBoundingBox() const noexcept override;
};

class SphereObject : public Object {
public:
    Sphere geometry;
    BasicMaterial material;

    SphereObject(Sphere _geometry, BasicMaterial _material = {}) noexcept;

private:
    std::optional<Pos> hitLocalRay(const Ray& ray) const noexcept override;

    ColorBounce deflectLocalRay(const Ray& ray) const noexcept override;

    AxisAlignedBox getLocalBoundingBox() const noexcept override;
};

class BoxObject : public Object {
public:
    Rectangle geometry;
    BasicMaterial material;

    BoxObject(Rectangle _geometry, BasicMaterial _material = {}) noexcept;

private:
    std::optional<Pos> hitLocalRay(const Ray& ray) const noexcept override;

    ColorBounce deflectLocalRay(const Ray& ray) const noexcept override;

    AxisAlignedBox getLocalBoundingBox() const noexcept override;
};

// Required member functions of Derived:
//  - float Derived::signedDistance(const Pos&) const noexcept;
//  - AxisAlignedBox Derived::localBoundingBox() const noexcept;
template<typename Derived>
class SDFObjectCRTP : public Object {
public:
    SDFObjectCRTP(BasicMaterial mat)
        : material(mat) {

    }

    BasicMaterial material;

private:
    inline std::optional<Pos> hitLocalRay(const Ray& ray) const noexcept override {
        const auto self = static_cast<const Derived*>(this);
        const auto bb = self->localBoundingBox();
        auto p0 = std::optional<Pos>{};
        if (inside(ray.pos, bb)) {
            p0 = ray.pos;
        } else {
            p0 = intersect(ray, bb);
        }
        if (!p0.has_value()) {
            return std::nullopt;
        }
        auto p = *p0;
        auto d = self->signedDistance(p);
        auto sign = d > 0.0f;
        for (std::size_t i = 0; i < 256; ++i) {
            auto d2 = self->signedDistance(p);
            if ((d2 > 0.0f) != sign) {
                d *= 0.5;
                p -= d * ray.dir;
                continue;
            }
            d = d2;
            if (std::abs(d) < 1e-4f) {
                return p;
            }
            p = p + d * ray.dir;
            if (!inside(p, bb)) {
                return std::nullopt;
            }
        }
        return std::nullopt;
    }

    inline ColorBounce deflectLocalRay(const Ray& ray) const noexcept override {
        const auto delta = 1e-3f;
        // TODO: could probably use 4 points in something like a tetrahedral
        // arrangement and some clever midpoint calculations, no?
        const auto dx = Vec{delta, 0.0f, 0.0f};
        const auto dy = Vec{0.0f, delta, 0.0f};
        const auto dz = Vec{0.0f, 0.0f, delta};
        const auto self = static_cast<const Derived*>(this);
        const auto n = (Vec{
            self->signedDistance(ray.pos + dx) - self->signedDistance(ray.pos - dx),
            self->signedDistance(ray.pos + dy) - self->signedDistance(ray.pos - dy),
            self->signedDistance(ray.pos + dz) - self->signedDistance(ray.pos - dz)
        } / delta).unit();

        return material.deflect(ray.dir,n);
    }

    inline AxisAlignedBox getLocalBoundingBox() const noexcept override {
        return static_cast<const Derived*>(this)->localBoundingBox();
    }
};

class FractalObject : public Object {
public:
    BasicMaterial material;

private:
    std::optional<Pos> hitLocalRay(const Ray& ray) const noexcept override;

    ColorBounce deflectLocalRay(const Ray& ray) const noexcept override;

    AxisAlignedBox getLocalBoundingBox() const noexcept override;

private:
    float signedDistance(const Pos&) const noexcept;
};

class Scene {
public:
    Scene();

    void addObject(std::unique_ptr<Object> object) noexcept;

    template<typename T, typename... Args>
    T& addObject(Args&&... args) {
        auto up = std::make_unique<T>(std::forward<Args>(args)...);
        auto& r = *up;
        addObject(std::move(up));
        return r;
    }

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

    void save(const std::string& path) const;
    static Image load(const std::string& path);

    const Color& operator()(std::size_t x, std::size_t y) const noexcept;
    Color& operator()(std::size_t x, std::size_t y) noexcept;

    void fill(const Color&) noexcept;

    Image& operator+=(const Image& other) noexcept;

private:
    size_t m_width;
    size_t m_height;
    std::vector<Color> m_data;
};

class Renderer {
public:
    Renderer(std::size_t width = 1024, std::size_t height = 1024) noexcept;
    ~Renderer();

    Renderer(const Renderer&) = delete;
    Renderer(Renderer&&) = delete;
    Renderer& operator=(const Renderer&) = delete;
    Renderer& operator=(Renderer&&) = delete;

    std::size_t width() const noexcept;
    std::size_t height() const noexcept;
    std::size_t numBounces() const noexcept;
    std::size_t samplesPerPixel() const noexcept;

    void setSize(size_t width, size_t height) noexcept;
    void setNumBounces(std::size_t) noexcept;
    void setSamplesPerPixel(std::size_t) noexcept;

    void startThreadPool(size_t numThreads = 0);

    void stopThreadPool();

    Image render(const Scene&, const Camera&);

private:
    std::size_t m_width;
    std::size_t m_height;
    std::size_t m_numBounces;
    std::size_t m_samplesPerPixel;
    mutable std::optional<std::barrier<>> m_renderBarrierMaybe;
    mutable std::atomic<size_t> m_nextTaskIndex;
    std::atomic<bool> m_timeToExit;
    mutable std::mutex m_sizeMutex;

    static const size_t s_pixelsPerTask = 8;

    std::vector<std::thread> m_threadPool;

    struct render_task {
        size_t xStart;
        size_t xEnd;
        size_t y;
    };

    static std::optional<render_task> makeTask(
        size_t taskIndex,
        size_t width,
        size_t height
    ) noexcept;

    struct RenderData {
        const Scene* scene;
        const Camera* camera;
        Image* image;
    };

    std::optional<RenderData> m_renderData;

    void doWork() const noexcept;
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