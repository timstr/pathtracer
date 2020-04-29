#include <PathTracer.hpp>
#include <RandomNumberGenerator.hpp>

#include <cassert>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <string>
#include <thread>

namespace {
    constexpr float epsilon = 1e-6f;
}

ColorBounce::ColorBounce(Color _emitted, Color _attenuation, Vec _direction) noexcept
    : emitted(_emitted)
    , attenuation(_attenuation)
    , direction(_direction) {

}

BasicGlossyMaterial::BasicGlossyMaterial(float _diffuseness, Color _reflectedColor, Color _emittedRadiance) noexcept
    : diffuseness(_diffuseness)
    , reflectedColor(_reflectedColor)
    , emittedRadiance(_emittedRadiance) {

}

ColorBounce BasicGlossyMaterial::deflect(const Vec& inbound, const Vec& normal) const noexcept {
    assert(diffuseness >= 0.0f && diffuseness <= 1.0f);
    const auto n = (inbound * normal) >= 0.0f ? -normal : normal;
    const auto diffuseBounce = randomPointOnHemisphereUniform(n);
    const auto specularBounce = bounce(inbound, n);
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


PerspectiveCamera::PerspectiveCamera(Affine transform, float aspectRatio, float fieldOfView) noexcept
    : m_transform(transform)
    , m_aspectRatio(aspectRatio)
    , m_fieldOfView(fieldOfView) 
    , m_focalDistance(10.0f)
    , m_focalBlurRadius(0.0f) {

}

const Affine& PerspectiveCamera::transform() const noexcept {
    return m_transform;
}

Affine& PerspectiveCamera::transform() noexcept {
    return m_transform;
}

float PerspectiveCamera::aspectRatio() const noexcept {
    return m_aspectRatio;
}

float PerspectiveCamera::fieldOfView() const noexcept {
    return m_fieldOfView;
}

float PerspectiveCamera::focalDistance() const noexcept {
    return m_focalDistance;
}

float PerspectiveCamera::focalBlurRadius() const noexcept {
    return m_focalBlurRadius;
}

void PerspectiveCamera::setTransform(Affine t) noexcept {
    m_transform = t;
}

void PerspectiveCamera::setAspectRatio(float r) noexcept {
    assert(r > 0.0f);
    m_aspectRatio = r;
}

void PerspectiveCamera::setFieldOfView(float fov) noexcept {
    assert(fov > -180.0f && fov < 180.0f);
    m_fieldOfView = fov;
}

void PerspectiveCamera::setFocalDistance(float d) noexcept {
    assert(d > 0.0f);
    m_focalDistance = d;
}

void PerspectiveCamera::setFocalBlurRadius(float r) noexcept {
    assert(r >= 0.0f);
    m_focalBlurRadius = r;
}

Ray PerspectiveCamera::getViewRay(float screenX, float screenY) const noexcept {
    const auto x = screenX * 2.0f - 1.0f;
    const auto y = screenY * 2.0f - 1.0f;
    const auto sp = m_aspectRatio > 1.0f ? Pos(x, y / m_aspectRatio) : Pos(x / m_aspectRatio, y);
    const auto dist = std::uniform_real_distribution<float>(0.0f, 1.0f);
    const auto blurAngle = 2.0f * 3.141592654f * dist(randomEngine());

    const auto [randX, randY] = randomPointInCircle();
    const auto blurRad = m_focalBlurRadius * std::max(m_aspectRatio, 1.0f / m_aspectRatio);
    const auto blurX = randX * blurRad;
    const auto blurY = randY * blurRad;
    const auto blurVec = Vec(blurX, blurY, 0.0f);
    const auto fovScale = std::tan(m_fieldOfView * 3.141592654f / 180.0f);
    const auto viewVecOrigin = (Vec(fovScale * sp.x, fovScale * sp.y, 1.0f) + (blurVec / m_focalDistance)).unit();
    const auto dir = m_transform * viewVecOrigin;
    const auto pos = m_transform * (sp - blurVec);
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

Renderer::Renderer(std::size_t width, std::size_t height) noexcept
    : m_width(width)
    , m_height(height)
    , m_numBounces(16)
    , m_samplesPerPixel(32)
    , m_numThreads(std::thread::hardware_concurrency()) {
    assert(m_width > 0);
    assert(m_height > 0);
    assert(m_numThreads > 0);
}

std::size_t Renderer::width() const noexcept {
    return m_width;
}

std::size_t Renderer::height() const noexcept {
    return m_height;
}

std::size_t Renderer::numBounces() const noexcept {
    return m_numBounces;
}

std::size_t Renderer::samplesPerPixel() const noexcept {
    return m_samplesPerPixel;
}

std::size_t Renderer::numThreads() const noexcept {
    return m_numThreads;
}

void Renderer::setWidth(std::size_t w) noexcept {
    assert(w > 0);
    m_width = w;
}

void Renderer::setHeight(std::size_t h) noexcept {
    assert(h > 0);
    m_height = h;
}

void Renderer::setNumBounces(std::size_t n) noexcept {
    assert(n > 0);
    m_numBounces = n;
}

void Renderer::setSamplesPerPixel(std::size_t s) noexcept {
    assert(s > 0);
    m_samplesPerPixel = s;
}

void Renderer::setNumThreads(std::size_t n) noexcept {
    assert(n > 0);
    m_numThreads = n;
}

Image Renderer::render(const Scene& scene, const Camera& camera) const {
    const auto maxWidth = static_cast<float>(m_width - 1);
    const auto maxHeight = static_cast<float>(m_height - 1);
    const auto sppInv = 1.0f / static_cast<float>(m_samplesPerPixel);

    std::mutex printMutex;
    std::size_t numRowsCompleted = 0;
    const auto reportRowCompleted = [&]() {
        auto lock = std::lock_guard{printMutex};
        ++numRowsCompleted;
        const auto pct = static_cast<float>(numRowsCompleted) / static_cast<float>(m_height);
        const auto totalWidth = 80;
        const auto progressBarWidth = static_cast<int>(totalWidth * pct);
        const auto remainingWidth = totalWidth - progressBarWidth;
        std::cout
            << "\r["
            << std::string(progressBarWidth, '=')
            << std::string(remainingWidth, ' ')
            << "] "
            << std::fixed
            << std::setprecision(2)
            << (pct * 100.0f)
            << '%';
        std::cout.flush();
    };

    auto img = Image(m_width, m_height);

    const auto doWork = [&](std::size_t yBegin, std::size_t yEnd, std::size_t yStep) {
        const auto deltaX = 1.0f / maxWidth;
        const auto deltaY = 1.0f / maxHeight;
        const auto jitterX = std::uniform_real_distribution<float>{ -0.5f * deltaX, 0.5f * deltaX };
        const auto jitterY = std::uniform_real_distribution<float>{ -0.5f * deltaY, 0.5f * deltaY };
        for (std::size_t y = yBegin; y < yEnd; y += yStep) {
            const auto py = static_cast<float>(y) / maxHeight;
            for (std::size_t x = 0; x < m_width; ++x) {
                const auto px = static_cast<float>(x) / maxWidth;
                auto acc = Color{};
                for (std::size_t i = 0; i < m_samplesPerPixel; ++i){
                    const auto sx = px + jitterX(randomEngine());
                    const auto sy = py + jitterY(randomEngine());
                    const auto c = scene.trace(camera.getViewRay(sx, sy), m_numBounces);
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
            reportRowCompleted();
        }
    };

    std::vector<std::thread> workers;
    for (std::size_t i = 0; i < m_numThreads; ++i) {
        workers.emplace_back(doWork, i, m_height, m_numThreads);
    }
    for (std::size_t i = 0; i < m_numThreads; ++i) {
        workers[i].join();
    }
    
    std::cout << '\n';

    return img;
}

Image ReinhardToneMapper::operator()(const Image& img) const noexcept {
    // https://www.researchgate.net/publication/255682296_Parameter_Estimation_for_Photographic_Tone_Reproduction
    auto acc = 0.0f;
    const auto w = img.width();
    const auto h = img.height();
    auto minLum = std::numeric_limits<float>::max();
    auto maxLum = std::numeric_limits<float>::min();
    for (std::size_t i = 0; i < w; ++i) {
        for (std::size_t j = 0; j < h; ++j) {
            const auto& c = img(i, j);
            const auto lum = 0.27f * c.r + 0.67f * c.b + 0.06f * c.b;
            minLum = std::min(minLum, lum);
            maxLum = std::max(minLum, lum);
            acc += std::log(lum + epsilon);
        }
    }

    const auto avgLum = std::exp(acc / static_cast<float>(w * h));
    const auto logAvgLum = std::log2(avgLum);
    const auto logMinLum = std::log2(minLum + epsilon);
    const auto logMaxLum = std::log2(maxLum + epsilon);

    const auto alpha = 0.18f * std::pow(4.0f, (2.0f * logAvgLum - logMinLum - logMaxLum) / (logMaxLum - logMinLum));

    const auto k = alpha / avgLum;

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