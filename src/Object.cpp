#include <Object.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>

const Affine& Object::transformation() const noexcept {
    return m_transformation;
}

const Affine& Object::inverseTransformation() const noexcept {
    return m_inverseTransformation;
}

void Object::setTransformation(const Affine& t) noexcept {
    m_transformation = t;
    auto invMaybe = m_transformation.inverse();
    assert(invMaybe.has_value());
    m_inverseTransformation = *invMaybe;
}

std::optional<Pos> Object::hitRay(const Ray& ray) const noexcept {
    const auto& inv = this->inverseTransformation();
    auto localRay = Ray(
        inv * ray.dir,
        inv * ray.pos
    );
    auto localHit = this->hitLocalRay(localRay);
    if (!localHit.has_value()) {
        return std::nullopt;
    }
    return this->transformation() * (*localHit);
}

ColorBounce Object::deflectRay(const Ray& ray) const noexcept {
    const auto& inv = this->inverseTransformation();
    auto localRay = Ray(
        inv * ray.dir,
        inv * ray.pos
    );
    auto deflected = this->deflectLocalRay(localRay);
    deflected.rayDirection = this->transformation() * deflected.rayDirection;
    return deflected;
}

AxisAlignedBox Object::getBoundingBox() const noexcept {
    auto bb = getLocalBoundingBox();
    const auto& c = bb.center;
    const auto& hs = bb.halfSize;
    const auto& T = transformation();
    auto pmin = Pos{
        std::numeric_limits<float>::max(),
        std::numeric_limits<float>::max(),
        std::numeric_limits<float>::max()
    };
    auto pmax = Pos{
        std::numeric_limits<float>::lowest(),
        std::numeric_limits<float>::lowest(),
        std::numeric_limits<float>::lowest()
    };
    const auto consider = [&](Vec corner){
        auto u = T * (c + corner);
        pmin.x = std::min(pmin.x, u.x);
        pmin.y = std::min(pmin.y, u.y);
        pmin.z = std::min(pmin.z, u.z);
        pmax.x = std::max(pmax.x, u.x);
        pmax.y = std::max(pmax.y, u.y);
        pmax.z = std::max(pmax.z, u.z);
    };
    consider(Vec{ hs.x,  hs.y,  hs.z});
    consider(Vec{ hs.x,  hs.y, -hs.z});
    consider(Vec{ hs.x, -hs.y,  hs.z});
    consider(Vec{ hs.x, -hs.y, -hs.z});
    consider(Vec{-hs.x,  hs.y,  hs.z});
    consider(Vec{-hs.x,  hs.y, -hs.z});
    consider(Vec{-hs.x, -hs.y,  hs.z});
    consider(Vec{-hs.x, -hs.y, -hs.z});

    return AxisAlignedBox(pmin, pmax);
}

TriangleObject::TriangleObject(Triangle _geometry, BasicMaterial _material)
    : geometry(_geometry)
    , material(_material) {

}

std::optional<Pos> TriangleObject::hitLocalRay(const Ray& ray) const noexcept {
    return intersect(ray, geometry);
}

ColorBounce TriangleObject::deflectLocalRay(const Ray& ray) const noexcept {
    return material.deflect(ray.dir, geometry.normal());
}

AxisAlignedBox TriangleObject::getLocalBoundingBox() const noexcept {
    const auto p0 = Pos{
        std::min({geometry.a.x, geometry.b.x, geometry.c.x}) - 1e-3f,
        std::min({geometry.a.y, geometry.b.y, geometry.c.y}) - 1e-3f,
        std::min({geometry.a.z, geometry.b.z, geometry.c.z}) - 1e-3f,
    };
    const auto p1 = Pos{
        std::max({geometry.a.x, geometry.b.x, geometry.c.x}) + 1e-3f,
        std::max({geometry.a.y, geometry.b.y, geometry.c.y}) + 1e-3f,
        std::max({geometry.a.z, geometry.b.z, geometry.c.z}) + 1e-3f,
    };
    return AxisAlignedBox{p0, p1};
}

SphereObject::SphereObject(Sphere _geometry, BasicMaterial _material) noexcept
    : geometry(_geometry)
    , material(_material) {

}

std::optional<Pos> SphereObject::hitLocalRay(const Ray& ray) const noexcept {
    return intersect(ray, geometry);
}

ColorBounce SphereObject::deflectLocalRay(const Ray& ray) const noexcept {
    return material.deflect(ray.dir, geometry.normal(ray.pos));
}

AxisAlignedBox SphereObject::getLocalBoundingBox() const noexcept {
    return AxisAlignedBox(
        Pos{0.0f, 0.0f, 0.0f},
        Vec{
            geometry.radius + 1e-3f,
            geometry.radius + 1e-3f,
            geometry.radius + 1e-3f
        }
    );
}

BoxObject::BoxObject(Rectangle _geometry, BasicMaterial _material) noexcept
    : geometry(_geometry)
    , material(_material) {

}

std::optional<Pos> BoxObject::hitLocalRay(const Ray& ray) const noexcept {
    return intersect(ray, geometry);
}

ColorBounce BoxObject::deflectLocalRay(const Ray& ray) const noexcept {
    return material.deflect(ray.dir, geometry.normal(ray.pos));
}

AxisAlignedBox BoxObject::getLocalBoundingBox() const noexcept {
    return AxisAlignedBox{
        Pos{0.0f, 0.0f, 0.0f},
        Vec{
            geometry.halfSize.x + 1e-3f,
            geometry.halfSize.y + 1e-3f,
            geometry.halfSize.z + 1e-3f
        }
    };
}

std::optional<Pos> FractalObject::hitLocalRay(const Ray& ray) const noexcept {
    const auto p0 = intersect(ray, getBoundingBox());

    if (!p0.has_value()) {
        return std::nullopt;
    }
    auto p = *p0;
    for (std::size_t i = 0; i < 256; ++i) {
        auto d = signedDistance(p);
        if (d < 1e-3f) {
            return p;
        }
        p = p + d * ray.dir;
    }
    return std::nullopt;
}

ColorBounce FractalObject::deflectLocalRay(const Ray& ray) const noexcept {
    const auto delta = 1e-3f;
    const auto dx = Vec{delta, 0.0f, 0.0f};
    const auto dy = Vec{0.0f, delta, 0.0f};
    const auto dz = Vec{0.0f, 0.0f, delta};
    const auto n = (Vec{
        signedDistance(ray.pos + dx) - signedDistance(ray.pos - dx),
        signedDistance(ray.pos + dy) - signedDistance(ray.pos - dy),
        signedDistance(ray.pos + dz) - signedDistance(ray.pos - dz)
    } / delta).unit();

    return material.deflect(ray.dir,n);
}

AxisAlignedBox FractalObject::getLocalBoundingBox() const noexcept {
    return AxisAlignedBox{
        Pos{0.0f, 0.0f, 0.0f},
        Vec{2.0f, 2.0f, 2.0f}
    };
}

float FractalObject::signedDistance(const Pos& p) const noexcept {
    /*
    // Mandelbulb
    const auto power = 8.0f;

    auto z = p.toVec();
    auto dr = 1.0f;
    auto r = 0.0f;

    for (int i = 0; i < 64; ++i) {
        r = z.norm();
        if (r > 2.0f) {
            break;
        }

        const auto theta = std::acos(z.z / r) * power;
        const auto phi = std::atan2(z.y, z.x) * power;
        const auto zr = std::pow(r, power);
        dr = std::pow(r, power - 1.0f) * power * dr + 1.0f;

        z = zr * Vec{
            std::sin(theta) * std::cos(phi),
            std::sin(theta) * std::sin(phi),
            std::cos(theta)
        };
        z = z + p.toVec();
    }
    const auto d = 0.5f * std::log(r) * r / dr;

    return d;
    */
    // 5 x 5 x 5 spheres

    const auto f = [](float v) {
        const auto l = 0.5f;
        const auto r = 2.0f;
        if (v < -l) {
            return v + l;
        } else if (v > l) {
            return v - l;
        } else {
            v *= r;
            return (v - std::round(v)) / r;
        }
    };

    const auto v = Vec{
        f(p.x),
        f(p.y),
        f(p.z)
    };

    return v.norm() - 0.2f;

}

