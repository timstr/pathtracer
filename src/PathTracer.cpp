#include <PathTracer.hpp>
#include <RandomNumberGenerator.hpp>

#include <cassert>

namespace {
    constexpr float epsilon = 1e-6f;
}


ColorBounce BasicGlossyMaterial::deflect(const Vec& inbound, const Vec& normal) const noexcept {
    assert(diffuseness >= 0.0f && diffuseness <= 1.0f);
    const auto diffuseBounce = randomHemisphereVectorUniform(normal);
    const auto specularBounce = inbound - (2.0f * ((inbound * normal) * normal));
    const auto ray = (specularBounce + diffuseness * (diffuseBounce - specularBounce)).unit();
    return {
        emittedRadiance,
        reflectedColor,
        ray
    };
}


TriangleObject::TriangleObject(Triangle _geometry, BasicGlossyMaterial _material)
    : geometry(_geometry)
    , material(_material) {
}

std::optional<float> TriangleObject::hit(const Ray& ray) const noexcept {
    return intersect(ray, geometry);
}

ColorBounce TriangleObject::deflect(const Ray& ray) const noexcept {
    return material.deflect(ray.dir, geometry.normal());
}

SphereObject::SphereObject(Sphere _geometry, BasicGlossyMaterial _material) noexcept
    : geometry(_geometry)
    , material(_material) {

}

std::optional<float> SphereObject::hit(const Ray& ray) const noexcept {
    return intersect(ray, geometry);
}

ColorBounce SphereObject::deflect(const Ray& ray) const noexcept {
    return material.deflect(ray.dir, geometry.normal(ray.pos));
}

void Scene::addObject(std::unique_ptr<Object> o) noexcept {
    m_objects.push_back(std::move(o));
}

Color Scene::trace(Ray ray, std::size_t depth) const noexcept {
    auto color = Color{};
    auto attenuation = Color{1.0f, 1.0f, 1.0f};
    for (std::size_t i = 0; i < depth; ++i) {
        auto closestT = std::numeric_limits<float>::max();
        const Object* hitObj = nullptr;
        for (const auto& obj : m_objects) {
            if (auto t = obj->hit(ray)) {
                if (*t < closestT) {
                    hitObj = obj.get();
                    closestT = *t;
                }
            }
        }
        if (!hitObj) {
            return color;
        }
        // TODO: overload +=
        ray.pos = ray.pos + closestT * ray.dir;
        const auto bounce = hitObj->deflect(ray);
        assert(bounce.emitted.r >= 0.0f);
        assert(bounce.emitted.g >= 0.0f);
        assert(bounce.emitted.b >= 0.0f);
        assert(bounce.attenuation.r >= 0.0f && bounce.attenuation.r <= 1.0f);
        assert(bounce.attenuation.g >= 0.0f && bounce.attenuation.g <= 1.0f);
        assert(bounce.attenuation.b >= 0.0f && bounce.attenuation.b <= 1.0f);
        assert(std::abs(1.0f - bounce.direction.norm()) < epsilon);
        // TODO: overload operators for Color
        color.r += bounce.emitted.r * attenuation.r;
        color.g += bounce.emitted.g * attenuation.g;
        color.b += bounce.emitted.b * attenuation.b;
        attenuation.r *= bounce.attenuation.r;
        attenuation.g *= bounce.attenuation.g;
        attenuation.b *= bounce.attenuation.b;
        ray.dir = bounce.direction;
    }
    return color;
}

OrthographicCamera::OrthographicCamera(Affine _transform, float _aspectRatio) noexcept
    : transform(_transform)
    , aspectRatio(_aspectRatio) {

}

Ray OrthographicCamera::getViewRay(float screenX, float screenY) const noexcept {
    const auto x = screenX * 2.0f - 1.0f;
    const auto y = screenY * 2.0f - 1.0f;
    const auto sp = aspectRatio > 1.0f ? Pos(x, y / aspectRatio) : Pos(x / aspectRatio, y);
    const auto pos = transform * sp;
    const auto dir = transform * Vec(0.02f * x, 0.02f * y, 1.0f);
    // const auto dir = transform * Vec(0.0f, 0.0f, 1.0f);
    return { dir, pos };
}

Image::Image(std::size_t width, std::size_t height)
    : m_width(width)
    , m_height(height)
    , m_data(width * height, Color{}) {

}

std::size_t Image::width() const noexcept {
    return m_width;
}

std::size_t Image::height() const noexcept {
    return m_height;
}

const Color& Image::operator()(std::size_t x, std::size_t y) const noexcept {
    assert(x < m_width);
    assert(y < m_height);
    return m_data[y * m_width + x];
}

Color& Image::operator()(std::size_t x, std::size_t y) noexcept {
    return const_cast<Color&>(const_cast<const Image*>(this)->operator()(x, y));
}

Image Renderer::render(const Scene& scene, const Camera& camera, std::size_t width, std::size_t height, std::size_t numBounces, std::size_t samplesPerPixel) const {
    const auto maxWidth = static_cast<float>(width - 1);
    const auto maxHeight = static_cast<float>(height - 1);
    const auto sppInv = 1.0f / static_cast<float>(samplesPerPixel);

    auto img = Image(width, height);
    for (std::size_t x = 0; x < width; ++x) {
        const auto sx = static_cast<float>(x) / maxWidth;
        for (std::size_t y = 0; y < height; ++y) {
            const auto sy = static_cast<float>(y) / maxHeight;
            auto acc = Color{};
            for (std::size_t i = 0; i < samplesPerPixel; ++i){
                 const auto c = scene.trace(camera.getViewRay(sx, sy), numBounces);
                 // TODO: operators
                 acc.r += c.r;
                 acc.g += c.g;
                 acc.b += c.b;
            }
            // TODO: operators
            img(x, y).r = acc.r * sppInv;
            img(x, y).g = acc.g * sppInv;
            img(x, y).b = acc.b * sppInv;
        }
    }
    return img;
}

Image ReinhardToneMapper::operator()(const Image& img) const noexcept {
    // TODO: this doesn't work
    auto acc = 0.0f;
    const auto w = img.width();
    const auto h = img.height();
    for (std::size_t i = 0; i < w; ++i) {
        for (std::size_t j = 0; j < h; ++j) {
            const auto& c = img(i, j);
            const auto lum = 0.27f * c.r + 0.67f * c.b + 0.06f * c.b;
            acc += std::log(lum);
        }
    }
    const auto k = 0.5f * std::exp(acc / static_cast<float>(w * h));

    auto ret = Image(w, h);

    for (std::size_t i = 0; i < w; ++i) {
        for (std::size_t j = 0; j < h; ++j) {
            const auto c = img(i, j);
            // TODO: operator overloading
            const auto scaled = Color{
                c.r * k,
                c.g * k,
                c.b * k
            };
            ret(i, j) = Color{
                scaled.r / (1.0f + scaled.r),
                scaled.g / (1.0f + scaled.g),
                scaled.b / (1.0f + scaled.b),
            };
        }
    }
    return ret;
}
