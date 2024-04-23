#pragma once

#include <BasicMaterial.hpp>
#include <ColorBounce.hpp>
#include <Geometry.hpp>
#include <LinearAlgebra.hpp>

#include <cmath>
#include <optional>

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

protected:
    Vec signedDistanceNormal(const Pos& pos) const noexcept {
        const auto delta = 1e-3f;
        // TODO: could probably use 4 points in something like a tetrahedral
        // arrangement and some clever midpoint calculations, no?
        const auto dx = Vec{delta, 0.0f, 0.0f};
        const auto dy = Vec{0.0f, delta, 0.0f};
        const auto dz = Vec{0.0f, 0.0f, delta};
        const auto self = static_cast<const Derived*>(this);
        return (Vec{
            self->signedDistance(pos + dx) - self->signedDistance(pos - dx),
            self->signedDistance(pos + dy) - self->signedDistance(pos - dy),
            self->signedDistance(pos + dz) - self->signedDistance(pos - dz)
        } / delta).unit();
    }

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
            p = p + std::abs(d) * ray.dir;
            if (!inside(p, bb)) {
                return std::nullopt;
            }
        }
        return std::nullopt;
    }

    inline ColorBounce deflectLocalRay(const Ray& ray) const noexcept override {
        const auto n = this->signedDistanceNormal(ray.pos);
        return material.deflect(ray.dir, n);
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
