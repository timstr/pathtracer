#include "..\include\PathTracer.hpp"
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

Color& Color::operator+=(const Color& other) noexcept {
    this->r += other.r;
    this->g += other.g;
    this->b += other.b;
    return *this;
}

Color operator+(const Color& l, const Color& r) noexcept {
    return Color(
        l.r + r.r,
        l.g + r.g,
        l.b + r.b
    );
}

Color operator*(float k, const Color& c) noexcept {
    return Color(
        k * c.r,
        k * c.g,
        k * c.b
    );
}

ColorBounce::ColorBounce(
    Color emitted,
    Color attenuation,
    Vec rayDirection,
    Vec normal
) noexcept
    : emitted(emitted)
    , attenuation(attenuation)
    , rayDirection(rayDirection)
    , normal(normal) {
    assert(std::abs(1.0f - rayDirection.norm()) < epsilon);
    assert(std::abs(1.0f - normal.norm()) < epsilon);
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
            return ColorBounce {
                Color{0.0f, 0.0f, 0.0f}, // TODO: emitted color?
                Color{1.0f, 1.0f, 1.0f}, // TODO: internal absorption (requires knowing distance traveled)
                inbound,
                normal
            };
        } else {
            // Total internal reflection
            return ColorBounce{
                Color{0.0f, 0.0f, 0.0f}, // HACK: yellow if total internal reflection
                Color{0.0f, 0.0f, 0.0f}, // TODO: internal absorption (requires knowing distance traveled)
                bounce(inbound, -normal),
                normal
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
                v,
                normal
            };
        } else {
            // Specular reflection
            const auto s = sinTheta * (1.0f - m_specularSharpness) * randomPointOnHemisphereUniform(normal);
            const auto v = (bounce(inbound, normal) + s).unit();
            return ColorBounce{
                m_emittedLuminance,
                m_reflectedAbsorption,
                v,
                normal
            };
        }
    } else {
        // Transmittance
        // TODO: specular/diffuse scattering
        const auto v = (inbound + ((inbound * normal) * (1.0f - 1.0f / m_indexOfRefraction) * normal)).unit();
        return ColorBounce{
            m_emittedLuminance,
            m_reflectedAbsorption,
            v,
            normal
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
        // auto closestT = std::numeric_limits<float>::max();
        // const Object* hitObj = nullptr;
        // for (const auto& obj : m_objects) {
        //     if (auto t = obj->hit(ray)) {
        //         if (*t < closestT) {
        //             hitObj = obj.get();
        //             closestT = *t;
        //         }
        //     }
        // }
        // if (!hitObj) {
        //     return color;
        // }

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
        assert(std::abs(1.0f - bounce.rayDirection.norm()) < epsilon);

        // TODO: overload operators for Color
        color.r += bounce.emitted.r * attenuation.r;
        color.g += bounce.emitted.g * attenuation.g;
        color.b += bounce.emitted.b * attenuation.b;
        attenuation.r *= bounce.attenuation.r;
        attenuation.g *= bounce.attenuation.g;
        attenuation.b *= bounce.attenuation.b;
        ray.dir = bounce.rayDirection;

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

void Image::fill(const Color& c) noexcept {
    for (auto& d : m_data) {
        d = c;
    }
}

Image& Image::operator+=(const Image& other) noexcept {
    assert(m_width == other.m_width);
    assert(m_height == other.m_height);
    for (size_t i = 0; i < m_data.size(); ++i) {
        m_data[i] += other.m_data[i];
    }
    return *this;
}

Renderer::Renderer(std::size_t width, std::size_t height) noexcept
    : m_width(width)
    , m_height(height)
    , m_numBounces(16)
    , m_samplesPerPixel(32)
    , m_renderBarrierMaybe(std::nullopt)
    , m_renderData(std::nullopt) {

    assert(m_width > 0);
    assert(m_height > 0);
}

Renderer::~Renderer() {
    this->stopThreadPool();
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

void Renderer::setNumBounces(std::size_t n) noexcept {
    assert(n > 0);
    m_numBounces = n;
}

void Renderer::setSamplesPerPixel(std::size_t s) noexcept {
    assert(s > 0);
    m_samplesPerPixel = s;
}

void Renderer::startThreadPool(size_t numThreads) {
    if (numThreads == 0) {
        numThreads = std::thread::hardware_concurrency();
    }
    std::cout << "Starting thread pool using " << numThreads << " threads"
        << std::endl;
    this->stopThreadPool();
    assert(!m_renderBarrierMaybe.has_value());
    assert(!m_renderData.has_value());
    assert(m_threadPool.size() == 0);
    // NOTE: the current thread must also arrive at the barrier to allow
    // phase transitions
    m_renderBarrierMaybe.emplace(numThreads + 1);
    m_timeToExit.store(false);
    for (std::size_t i = 0; i < numThreads; ++i) {
        m_threadPool.emplace_back([this]{ this->doWork(); });
    }
}

void Renderer::stopThreadPool() {
    if (m_threadPool.size() == 0) {
        assert(!m_renderBarrierMaybe.has_value());
        assert(!m_renderData.has_value());
        return;
    }
    assert(m_renderBarrierMaybe.has_value());
    assert(m_timeToExit.load() == false);
    m_timeToExit.store(true);
    m_renderBarrierMaybe->arrive_and_wait();
    for (auto& t : m_threadPool) {
        t.join();
    }
    m_threadPool.clear();
    m_renderBarrierMaybe.reset();
}

void Renderer::doWork() const noexcept {
    while (true) {
        // Phase 1: waiting for work
        assert(m_renderBarrierMaybe.has_value());
        m_renderBarrierMaybe->arrive_and_wait();
        if (m_timeToExit.load()) {
            return;
        }

        // Phase 2: doing the work
        assert(m_renderData.has_value());
        const auto& camera = *m_renderData->camera;
        const auto& scene = *m_renderData->scene;
        auto& img = *m_renderData->image;

        const auto maxWidth = static_cast<float>(m_width - 1);
        const auto maxHeight = static_cast<float>(m_height - 1);
        const auto sppInv = 1.0f / static_cast<float>(m_samplesPerPixel);
        const auto deltaX = 1.0f / maxWidth;
        const auto deltaY = 1.0f / maxHeight;
        const auto jitterX = std::uniform_real_distribution<float>{ -0.5f * deltaX, 0.5f * deltaX };
        const auto jitterY = std::uniform_real_distribution<float>{ -0.5f * deltaY, 0.5f * deltaY };
        while (true) {
            const auto taskIndex = m_nextTaskIndex.fetch_add(1);
            auto taskMaybe = this->make_task(taskIndex);
            if (!taskMaybe.has_value()) {
                break;
            }
            const auto& task = *taskMaybe;
            const auto y = task.y;
            const auto py = static_cast<float>(y) / maxHeight;
            for (std::size_t x = task.xStart; x < task.xEnd; ++x) {
                const auto px = static_cast<float>(x) / maxWidth;
                auto acc = Color{};
                for (std::size_t i = 0; i < m_samplesPerPixel; ++i){
                    const auto sx = px + jitterX(randomEngine());
                    const auto sy = py + jitterY(randomEngine());
                    const auto viewRay = camera.getViewRay(sx, sy);
                    const auto c = scene.trace(viewRay, m_numBounces);
                    // TODO: operators
                    acc.r += c.r;
                    acc.g += c.g;
                    acc.b += c.b;
                }
                // TODO: operators
                auto& pixel = img(x, y);
                pixel.r = acc.r * sppInv;
                pixel.g = acc.g * sppInv;
                pixel.b = acc.b * sppInv;
            }
        }

        assert(m_renderBarrierMaybe.has_value());
        m_renderBarrierMaybe->arrive_and_wait();
        if (m_timeToExit.load()) {
            return;
        }
    }
}

std::optional<Renderer::render_task> Renderer::make_task(
    size_t taskIndex
) const noexcept {
    assert(m_renderData.has_value());
    const auto& imgPtr = m_renderData->image;
    assert(imgPtr != nullptr);
    const auto tasksPerRow = imgPtr->width() / s_pixelsPerTask;
    const auto xStart = s_pixelsPerTask * (taskIndex % tasksPerRow);
    const auto y = taskIndex / tasksPerRow;
    if (y >= imgPtr->height()) {
        return std::nullopt;
    }
    const auto xEnd = std::min(
        xStart + s_pixelsPerTask,
        imgPtr->width()
    );
    return render_task {
        xStart,
        xEnd,
        y
    };
}

Image Renderer::render(
    const Scene& scene,
    const Camera& camera
) {
    if (m_threadPool.size() == 0) {
        this->startThreadPool();
    }

    auto img = Image(m_width, m_height);

    assert(!m_renderData.has_value());
    m_renderData.emplace(RenderData{
        &scene,
        &camera,
        &img
    });
    m_nextTaskIndex.store(0);
    assert(m_renderBarrierMaybe.has_value());
    assert(!m_timeToExit.load());
    // Phase 1: threads are all waiting at the barrier to start working.
    // This last thread arriving starts them
    m_renderBarrierMaybe->arrive_and_wait();
    // Phase 2: threads are all busy doing work.
    // When the last thread arrives, all work is done.
    m_renderBarrierMaybe->arrive_and_wait();

    m_renderData.reset();

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

BoxObject::BoxObject(Box _geometry, BasicMaterial _material) noexcept
    : geometry(_geometry)
    , material(_material) {

}

std::optional<float> BoxObject::hit(const Ray& ray) const noexcept {
    return intersect(ray, geometry);
}

ColorBounce BoxObject::deflect(const Ray& ray) const noexcept {
    return material.deflect(ray.dir, geometry.normal(ray.pos));
}

Box BoxObject::getBoundingBox() const noexcept {
    return geometry;
}

std::optional<float> FractalObject::hit(const Ray& ray) const noexcept {
    auto r = ray;
    auto t = 0.0f;
    auto pd = std::numeric_limits<float>::min();
    for (std::size_t i = 0; i < 128; ++i) {
        auto d = signedDistance(r.pos);
        if (d < pd && d < 1e-6f) {
            return t;
        }
        t += d;
        r.pos = r.pos + d * r.dir;
        pd = d;
    }
    return std::nullopt;
}

ColorBounce FractalObject::deflect(const Ray& ray) const noexcept {
    const auto delta = 1e-6f;
    const auto dx = Vec{delta, 0.0f, 0.0f};
    const auto dy = Vec{0.0f, delta, 0.0f};
    const auto dz = Vec{0.0f, 0.0f, delta};
    const auto n = Vec{
        signedDistance(ray.pos + dx) - signedDistance(ray.pos - dx),
        signedDistance(ray.pos + dy) - signedDistance(ray.pos - dy),
        signedDistance(ray.pos + dz) - signedDistance(ray.pos - dz)
    }.unit();

    return material.deflect(ray.dir, n);
}

Box FractalObject::getBoundingBox() const noexcept {
    return Box{
        Pos{0.0f, 0.0f, 0.0f},
        Vec{3.0f, 3.0f, 3.0f}
    };
}

float FractalObject::signedDistance(const Pos& p) const noexcept {
    // Mandelbulb
    const auto power = 8.0f;

    auto z = p.toVec();
    auto dr = 1.0f;
    auto r = 0.0f;

    for (int i = 0; i < 15; ++i) {
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
    return 0.5f * std::log(r) * r / dr;

    // 5 x 5 x 5 spheres
    /*
    const auto f = [](float v) {
        const auto l = 1.0f;
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
    */
}
