#include <PathTracer.hpp>
#include <RandomNumberGenerator.hpp>

#include <cassert>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <numeric>
#include <string>
#include <thread>

namespace {
    constexpr float epsilon = 1e-6f;
}

Color::Color(float _r, float _g, float _b) noexcept
    : r(_r)
    , g(_g)
    , b(_b) {

}

ColorBounce::ColorBounce(Color _emitted, Color _attenuation, Vec _direction) noexcept
    : emitted(_emitted)
    , attenuation(_attenuation)
    , direction(_direction) {
    assert(std::abs(1.0f - direction.norm()) < epsilon);
    assert(attenuation.r >= 0.0f && attenuation.r <= 1.0f);
    assert(attenuation.g >= 0.0f && attenuation.g <= 1.0f);
    assert(attenuation.b >= 0.0f && attenuation.b <= 1.0f);
}

BasicMaterial::BasicMaterial() noexcept
    : m_diffuseReflection(0.3f)
    , m_specularReflection(0.3f)
    , m_specularSharpness(0.9f)
    , m_reflectedAbsorption(Color{1.0f, 1.0f, 1.0f})
    , m_emittedLuminance(Color{0.0f, 0.0f, 0.0f})
    , m_transmittance(0.2f)
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
            return ColorBounce{
                Color{0.0f, 0.0f, 0.0f}, // TODO: emitted color?
                Color{1.0f, 1.0f, 1.0f}, // TODO: internal absorption (requires knowing distance traveled)
                inbound
            };
        } else {
            // Total internal reflection
            return ColorBounce{
                Color{10.0f, 10.0f, 0.0f}, // HACK: yellow if total internal reflection
                Color{0.0f, 0.0f, 0.0f}, // TODO: internal absorption (requires knowing distance traveled)
                bounce(inbound, -normal)
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
            const auto v = randomPointOnHemisphereUniform(normal);
            return ColorBounce{
                m_emittedLuminance,
                m_reflectedAbsorption,
                v
            };
        } else {
            // Specular reflection
            const auto s = sinTheta * (1.0f - m_specularSharpness) * randomPointOnHemisphereUniform(normal);
            const auto v = (bounce(inbound, normal) + s).unit();
            return ColorBounce{
                m_emittedLuminance,
                m_reflectedAbsorption,
                v
            };
        }
    } else {
        // Transmittance
        const auto v = (inbound + ((inbound * normal) * (1.0f - 1.0f / m_indexOfRefraction) * normal)).unit();
        return ColorBounce{
            m_emittedLuminance,
            m_reflectedAbsorption,
            v
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


TriangleObject::TriangleObject(Triangle _geometry, BasicMaterial _material)
    : geometry(_geometry)
    , material(_material) {

}

std::optional<float> TriangleObject::hit(const Ray& ray) const noexcept {
    return intersect(ray, geometry);
}

ColorBounce TriangleObject::deflect(const Ray& ray) const noexcept {
    return material.deflect(ray.dir, geometry.normal());
}

Box TriangleObject::getBoundingBox() const noexcept {
    const auto p0 = Pos{
        std::min({geometry.a.x, geometry.b.x, geometry.c.x}),
        std::min({geometry.a.y, geometry.b.y, geometry.c.y}),
        std::min({geometry.a.z, geometry.b.z, geometry.c.z}),
    };
    const auto p1 = Pos{
        std::max({geometry.a.x, geometry.b.x, geometry.c.x}),
        std::max({geometry.a.y, geometry.b.y, geometry.c.y}),
        std::max({geometry.a.z, geometry.b.z, geometry.c.z}),
    };
    return Box{p0, p1};
}

SphereObject::SphereObject(Sphere _geometry, BasicMaterial _material) noexcept
    : geometry(_geometry)
    , material(_material) {

}

std::optional<float> SphereObject::hit(const Ray& ray) const noexcept {
    return intersect(ray, geometry);
}

ColorBounce SphereObject::deflect(const Ray& ray) const noexcept {
    return material.deflect(ray.dir, geometry.normal(ray.pos));
}

Box SphereObject::getBoundingBox() const noexcept {
    return Box(geometry.center, geometry.radius * Vec{1.0f, 1.0f, 1.0f});
}

void Scene::addObject(std::unique_ptr<Object> object) noexcept {
    m_objects.push_back(std::move(object));
}

Scene::Scene()
    : m_objectTree(std::vector<const Object*>{}) {

}

Color Scene::trace(Ray ray, std::size_t depth) const noexcept {
    auto color = Color{};
    auto attenuation = Color{1.0f, 1.0f, 1.0f};
    for (std::size_t i = 0; i < depth; ++i) {
        // Brute force search
        /*
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
        */

        
        // Tree search
        const auto hit = m_objectTree.hit(ray);
        if (!hit) {
            return color;
        }
        const auto& [closestT, hitObj] = *hit;
        

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

        if (attenuation.r + attenuation.b + attenuation.g < epsilon) {
            return color;
        }
    }
    return color;
}

void Scene::updateGeometry() {
    std::vector<const Object*> v;
    v.reserve(m_objects.size());
    std::transform(
        m_objects.begin(),
        m_objects.end(),
        std::back_inserter(v),
        [](const std::unique_ptr<Object>& o) {
            return o.get();
        }
    );
    m_objectTree = ObjectTree::Tree{v};
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
    const auto dir = (m_transform * viewVec).unit();
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

    const auto startTime = std::chrono::steady_clock::now();

    const auto prettyPrintSeconds = [](std::size_t totalSeconds) {
        const auto seconds = totalSeconds % 60;
        const auto minutes = totalSeconds / 60 % 60;
        const auto hours = totalSeconds / 3600;
        std::cout
            << std::setfill('0') << std::setw(2)
            << hours << ':' << std::setw(2) << minutes << ':' << std::setw(2) << seconds;
    };

    std::mutex printMutex;
    std::size_t numRowsCompleted = 0;
    const auto reportRowCompleted = [&]() {
        auto lock = std::lock_guard{printMutex};
        ++numRowsCompleted;
        const auto pct = static_cast<float>(numRowsCompleted) / static_cast<float>(m_height);
        const auto totalWidth = 80;
        const auto progressBarWidth = static_cast<int>(totalWidth * pct);
        const auto remainingWidth = totalWidth - progressBarWidth;
        const auto currTime = std::chrono::steady_clock::now();
        const auto totalSeconds = std::chrono::duration_cast<std::chrono::seconds>(currTime - startTime).count();
        
        std::cout
            << "\r["
            << std::string(progressBarWidth, '=')
            << std::string(remainingWidth, ' ')
            << "] "
            << std::fixed << std::setprecision(2)
            << std::setfill('0') << std::setw(2)
            << (pct * 100.0f) << "% - ";
        prettyPrintSeconds(totalSeconds);

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
    
    std::cout << "\nRender completed after ";
    const auto currTime = std::chrono::steady_clock::now();
    const auto totalSeconds = std::chrono::duration_cast<std::chrono::seconds>(currTime - startTime).count();
    prettyPrintSeconds(totalSeconds);
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
