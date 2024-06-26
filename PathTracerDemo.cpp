// To silence warning C4996: 'localtime': This function or variable may be
// unsafe. Consider using localtime_s instead. To disable deprecation, use
// _CRT_SECURE_NO_WARNINGS.
#define _CRT_SECURE_NO_WARNINGS

#include <Camera.hpp>
#include <Geometry.hpp>
#include <Image.hpp>
#include <LinearAlgebra.hpp>
#include <Object.hpp>
#include <RandomNumberGenerator.hpp>
#include <Renderer.hpp>
#include <RenderSettings.hpp>
#include <Scene.hpp>
#include <ToneMapper.hpp>

#include <SFML/Graphics.hpp>

#include <cassert>
#include <ctime>
#include <functional>
#include <iostream>
#include <iomanip>
#include <sstream>

float triangle(float x) noexcept {
    return 2.0f * std::abs(x - std::floor(x + 0.5f));
}

template<typename T>
float hashf(T t) {
    const auto h = std::hash<T>{}(t) % std::size_t{0xff};
    const auto v = 3.141592654f / static_cast<float>(h);
    return v - std::floor(v);
}

const Affine& noiseTransform() noexcept {
    static auto t = []{
        auto dist = std::uniform_real_distribution<float>(0.0f, 1.0f);
        const auto t1 = dist(randomEngine());
        const auto t2 = dist(randomEngine());
        const auto t3 = dist(randomEngine());
        const auto r1 = dist(randomEngine());
        const auto r2 = dist(randomEngine());
        const auto r3 = dist(randomEngine());
        std::cout << t1 << ' ' << t2 << ' ' << t3 << std::endl;
        std::cout << r1 << ' ' << r2 << ' ' << r3 << std::endl;
        return Affine::Translation(Vec(t1, t2, t3))
            * Linear::RotationX(r1)
            * Linear::RotationY(r2)
            * Linear::RotationZ(r3);
    }();
    return t;
}

float noise(Pos location) noexcept {
    auto v = 0.0f;
    const auto& transform = noiseTransform();
    auto k = 0.1f;
    auto kk = 0.85f;
    for (int i = 0; i < 32; ++i) {
        location = transform * location;
        v += k * triangle(location.x);
        location.x *= 1.1f;
        location.y *= 1.1f;
        location.z *= 1.1f;
        k *= kk;
    }
    return v;
}

float smin(float a, float b, float k) noexcept {
    return -std::log2(
        std::exp2(-k * a) + std::exp2(-k * b)
    ) / k;
}

class RoughSphereObject : public SDFObjectCRTP<RoughSphereObject> {
public:
    RoughSphereObject(Sphere geometry, BasicMaterial mat1, BasicMaterial mat2)
        : SDFObjectCRTP(mat1)
        , m_geometry(geometry)
        , m_material2(mat2) {

    }

    float signedDistanceCraters(const Pos& pos) const noexcept {
        return std::max(0.0f, -1.0f + 3.0f * noise((1.0f * pos.toVec()).toPos()));
    }

    float signedDistanceBumps(const Pos& pos) const noexcept {
        return -0.15f + 0.5f * noise((4.0f * pos.toVec()).toPos());
    }

    float signedDistance(const Pos& pos) const noexcept {
        const auto r = 0.8f * m_geometry.radius;
        const auto sdBox = Rectangle(Vec(r, r, r)).signedDistance(pos);
        const auto sdSphere = m_geometry.signedDistance(pos);
        const auto sd = smin(sdBox, sdSphere, 8.0f);
        const auto craters = this->signedDistanceCraters(pos);
        const auto bumps = this->signedDistanceBumps(pos);
        // return smin(sd + craters, sd + craters + bumps, 32.0f);
        return std::min(sd + craters, sd + craters + bumps);
    }

    ColorBounce deflectLocalRay(const Ray& ray) const noexcept override {
        const auto n = this->signedDistanceNormal(ray.pos);
        auto m = this->material;
        auto noise1 = this->signedDistanceNormal((-ray.pos.toVec() * 1.724674f + Vec(3.474577f, 1.734674f, 2.713652f)).toPos());
        auto noise2 = this->signedDistanceNormal(( ray.pos.toVec() * 2.446744f + Vec(1.734674f, 2.713652f, 4.137371f)).toPos());
        // const auto noise1 = this->signedDistanceNormal(ray.pos);

        noise1 = Vec(0.5f, 0.5f, 0.5f) + 0.5f * noise1;
        noise2 = Vec(0.5f, 0.5f, 0.5f) + 0.5f * noise2;

        m.setReflectedAbsorption(Color(noise1.x, noise1.y, noise1.z));
        m.setDiffuseReflection(noise2.x);
        m.setSpecularReflection(1.0f - noise2.x);
        m.setSpecularSharpness(noise2.y);
        m.setEmittedLuminance(0.125f * noise2.z * m.reflectedAbsorption());
        return m.deflect(ray.dir, n);

        // const auto craters = this->signedDistanceCraters(ray.pos);
        // const auto bumps = this->signedDistanceBumps(ray.pos);
        // if (craters < bumps) {
        //     return this->material.deflect(ray.dir, n);
        // }
        // return this->m_material2.deflect(ray.dir, n);

        // return ColorBounce(
        //     Color(1.0f, 0.0f, 1.0f),
        //     Color(0.0f, 0.0f, 0.0f),
        //     bounce(ray.dir, n),
        //     n
        // );
    }

    AxisAlignedBox localBoundingBox() const noexcept {
        return AxisAlignedBox(
            Pos(0.0f, 0.0f, 0.0f),
            1.2f * m_geometry.radius * Vec(1.0f, 1.0f, 1.0f)
        );
    }

private:
    Sphere m_geometry;
    BasicMaterial m_material2;
};

class RoughBoxObject : public SDFObjectCRTP<RoughBoxObject> {
public:
    RoughBoxObject(Rectangle geometry, BasicMaterial mat)
        : SDFObjectCRTP(mat)
        , m_geometry(geometry) {

    }


    float signedDistance(const Pos& pos) const noexcept {
        // return m_geometry.signedDistance(pos) + 0.1f * noise((50.0f * pos.toVec()).toPos());
        return m_geometry.signedDistance(pos) - 0.01f;
    }

    AxisAlignedBox localBoundingBox() const noexcept {
        return AxisAlignedBox(
            Pos(0.0f, 0.0f, 0.0f),
            1.1f * m_geometry.halfSize
        );
    }

private:
    Rectangle m_geometry;
};

class InfiniteLightSource : public Object {
public:
    Color color;
    Vec direction;
    float focus;

    InfiniteLightSource(Color _color, Vec _direction, float _focus) : color(_color), direction(_direction), focus(_focus) {}

    std::optional<Pos> hitLocalRay(const Ray& ray) const noexcept override {
        return ray.pos + ray.dir * 1e10f;
    }

    ColorBounce deflectLocalRay(const Ray& ray) const noexcept override {
        const auto k = ray.dir * direction >= focus ? 1.0f : 0.0f;
        //const auto k = std::pow(std::max(0.0f, ray.dir * direction), focus);
        return ColorBounce {
            Color{
                color.r * k,
                color.g * k,
                color.b * k
            },
            Color{0.0f, 0.0f, 0.0f},
            ray.dir,
            ray.dir
        };
    }
};

Image scaleImage(const Image& img, float k) noexcept {
    auto ret = Image(img.width(), img.height());
    for (std::size_t j = 0; j < img.height(); ++j) {
        for (std::size_t i = 0; i < img.width(); ++i) {
            ret(i, j) = k * img(i, j);
        }
    }
    return ret;
}

void copyToSFImage(const Image& src, sf::Image& dst) {
    auto s = dst.getSize();
    if (s.x != src.width() || s.y != src.height()) {
        dst.create(
            static_cast<unsigned int>(src.width()),
            static_cast<unsigned int>(src.height())
        );
    }
    assert(dst.getSize().x == src.width());
    assert(dst.getSize().y == src.height());
    for (unsigned x = 0; x < src.width(); ++x) {
        for (unsigned y = 0; y < src.height(); ++y) {
            auto color = src(x, y);
            auto px = sf::Color(
                static_cast<std::uint8_t>(std::clamp(color.r, 0.0f, 1.0f) * 255.0f),
                static_cast<std::uint8_t>(std::clamp(color.g, 0.0f, 1.0f) * 255.0f),
                static_cast<std::uint8_t>(std::clamp(color.b, 0.0f, 1.0f) * 255.0f)
            );
            dst.setPixel(x, y, px);
        }
    }
}


int main() {

    auto s = Scene{};

    // Random spheres
    // {
    //     auto pdist = std::uniform_real_distribution<float>{0.0f, 1.0f};
    //     auto sdist = std::uniform_real_distribution<float>{-1.0f, 1.0f};
    //     for (std::size_t i = 0; i < 200; ++i) {
    //         const auto geom = Sphere{
    //             0.1f + 0.7f * pdist(randomEngine())
    //         };
    //         // const auto geom = 0.1f + 0.4f * pdist(randomEngine());
    //         auto m = BasicMaterial{};
    //         m.setDiffuseReflection(1.0f);
    //         m.setSpecularReflection(0.0f);
    //         m.setTransmittance(0.0f);
    //         const auto k = 30.0f * std::pow(pdist(randomEngine()), 3.0f);
    //         if (pdist(randomEngine()) > 0.9) {
    //             m.setEmittedLuminance(Color{
    //                 k * pdist(randomEngine()),
    //                 k * pdist(randomEngine()),
    //                 k * pdist(randomEngine())
    //             });
    //         }
    //         auto& sphere = s.addObject<SphereObject>(geom, m);
    //         sphere.setTransformation(Affine::Translation(Vec{
    //             10.0f * sdist(randomEngine()),
    //             10.0f * sdist(randomEngine()),
    //             10.0f * sdist(randomEngine())
    //         }));
    //     }
    // }

    // Cornell box
    {
        auto matte = BasicMaterial();
        matte.setDiffuseReflection(1.0f);
        // matte.setEmittedLuminance(Color(0.1f, 0.1f, 0.1f));
        matte.setEmittedLuminance(Color(0.0f, 0.0f, 0.0f));
        matte.setSpecularReflection(0.0f);
        matte.setTransmittance(0.0f);

        auto matteRed = matte;
        auto matteGreen = matte;
        auto matteGray = matte;
        auto matteWhite = matte;
        matteRed.setReflectedAbsorption(Color(0.8f, 0.1f, 0.1f));
        matteGreen.setReflectedAbsorption(Color(0.1f, 0.8f, 0.1f));
        matteGray.setReflectedAbsorption(Color(0.2f, 0.2f, 0.3f));
        matteWhite.setReflectedAbsorption(Color(1.0f, 1.0f, 1.0f));

        auto glowing = matte;
        auto antiGlowing = matte;
        glowing.setEmittedLuminance(Color(50.0f, 40.0f, 30.0f));
        antiGlowing.setEmittedLuminance(Color(-25.0f, -20.0f, -15.0f));

        auto mirror = matte;
        mirror.setDiffuseReflection(0.05f);
        mirror.setSpecularReflection(0.95f);
        mirror.setSpecularSharpness(0.95f);
        mirror.setReflectedAbsorption(Color(0.9f, 0.9f, 0.9f));

        const auto d = 0.01f;

        auto& lightSource = s.addObject<BoxObject>(Rectangle(Vec{1.0f, d, 1.0f}), glowing);
        lightSource.setTransformation(Affine::Translation(Vec{0.0f, -5.0f + 2.0f * d, 0.0f}));

        // auto& lightSink = s.addObject<BoxObject>(Rectangle(Vec{1.0f, d, 1.0f}), antiGlowing);
        // lightSink.setTransformation(Affine::Translation(Vec{0.0f, 5.0f - 2.0f * d, 0.0f}));

        auto& rearWall = s.addObject<BoxObject>(Rectangle(Vec{5.0f, 5.0f, d}), matteWhite);
        rearWall.setTransformation(Affine::Translation(Vec{0.0f, 0.0f, 5.0f}));

        auto& floor = s.addObject<BoxObject>(Rectangle(Vec{5.0f, d, 5.0f}), matteWhite);
        floor.setTransformation(Affine::Translation(Vec{0.0f, 5.0f, 0.0f}));

        // ceiling
        auto& ceiling = s.addObject<BoxObject>(Rectangle(Vec{5.0f, d, 5.0f}), matteWhite);
        ceiling.setTransformation(Affine::Translation(Vec{0.0f, -5.0f, 0.0f}));

        // left wall
        auto& leftWall = s.addObject<BoxObject>(Rectangle(Vec{d, 5.0f, 5.0f}), matteRed);
        leftWall.setTransformation(Affine::Translation(Vec{-5.0f, 0.0f, 0.0f}));

        // right wall
        auto& rightWall = s.addObject<BoxObject>(Rectangle(Vec{d, 5.0f, 5.0f}), matteGreen);
        rightWall.setTransformation(Affine::Translation(Vec{5.0f, 0.0f, 0.0f}));

        // front box
        // auto& box = s.addObject<BoxObject>(Rectangle(Vec{1.5f, 2.0, 1.5f}), matteWhite);
        // box.setTransformation(Affine::Translation(Vec{-1.5f, 3.0f, 1.0f}) * Linear::RotationY(0.4136f));
        // box.material.setReflectedAbsorption(Color(0.5f, 0.5f, 1.0f));

        auto& whiteSphere1 = s.addObject<SphereObject>(Sphere(1.5f), matteWhite);
        whiteSphere1.setTransformation(Affine::Translation(Vec{-3.0f, 3.5f, -2.5f}));

        auto& whiteSphere2 = s.addObject<SphereObject>(Sphere(1.5f), matteWhite);
        whiteSphere2.setTransformation(Affine::Translation(Vec{3.0f, -3.5f, -2.5f}));

        auto& mirrorSphere1 = s.addObject<SphereObject>(Sphere(1.5f), mirror);
        mirrorSphere1.setTransformation(Affine::Translation(Vec{3.0f, 3.5f, 2.5f}));

        auto& mirrorSphere2 = s.addObject<SphereObject>(Sphere(1.5f), mirror);
        mirrorSphere2.setTransformation(Affine::Translation(Vec{-3.0f, -3.5f, 2.5f}));

        // auto& fractal = s.addObject<FractalObject>();
        // fractal.material = matteGray;
        // fractal.material.setEmittedLuminance(Color(0.05f, 0.1f, 0.2f));
        // fractal.material.setReflectedAbsorption(Color(0.2f, 0.4f, 0.8f));
        // fractal.material.setDiffuseReflection(1.0f);
        // fractal.material.setSpecularReflection(1.0f);
        // fractal.setTransformation(
        //     Affine::Translation(Vec{-1.5f, -1.0f, 1.0f})
        //     // * Linear::RotationY(0.6f)
        //     // * Linear::RotationX(0.3f)
        //     // * Linear::Scale(1.5f)
        // );
    }

    // Globe of cubes (or glob of spheres)
    // {
    //     const auto res = std::size_t{16};
    //     const auto size = 3.0f;
    //     const auto dsize = size / static_cast<float>(res);
    //     const auto mapToSpace = [&](std::size_t idx) {
    //         return (static_cast<float>(idx) / static_cast<float>(res - 1) * 2.0f - 1.0f) * size;
    //     };


    //     auto dist = std::uniform_real_distribution<float>{0.0f, 1.0f};
    //     auto sdist = std::uniform_real_distribution<float>{-1.0f, 1.0f};

    //     for (std::size_t i = 0; i < res; ++i){
    //         const auto x = mapToSpace(i);
    //         for (std::size_t j = 0; j < res; ++j){
    //             const auto y = mapToSpace(j);
    //             for (std::size_t k = 0; k < res; ++k){
    //                 const auto z = mapToSpace(k);
    //                 if (x * x + y * y + z * z > size) {
    //                     continue;
    //                 }

    //                 // if (dist(randomEngine()) < 0.6f) {
    //                 //     continue;
    //                 // }

    //                 const auto g = Rectangle(dsize * Vec(
    //                     (0.8f + 0.2f * dist(randomEngine())),
    //                     (0.8f + 0.2f * dist(randomEngine())),
    //                     (0.8f + 0.2f * dist(randomEngine()))
    //                 ));
    //                 // const auto g = Sphere(0.0f + 0.3f * dist(randomEngine()));

    //                 auto m = BasicMaterial{};
    //                 if (dist(randomEngine()) < 0.9f) {
    //                     m.setDiffuseReflection(0.05f);// * dist(randomEngine()));
    //                     m.setSpecularReflection(0.95f);//1.0f * dist(randomEngine()));
    //                     m.setSpecularSharpness(0.9f);//1.0f * dist(randomEngine()));
    //                     m.setTransmittance(0.0f);
    //                     m.setReflectedAbsorption(Color{1.0f, 0.9f, 0.2f});
    //                     m.setEmittedLuminance(Color{0.0f, 0.0f, 0.0f});
    //                 } else {
    //                     m.setDiffuseReflection(0.9f);
    //                     m.setSpecularReflection(0.1f);
    //                     m.setSpecularSharpness(0.95f);
    //                     m.setTransmittance(0.0f);
    //                     m.setIndexOfRefraction(1.55f);
    //                     m.setReflectedAbsorption(Color{1.0f, 1.0f, 1.0f});
    //                     // m.setEmittedLuminance(Color{0.10f, 0.25f, 0.0f});
    //                     m.setEmittedLuminance(Color{
    //                         dist(randomEngine()),
    //                         dist(randomEngine()),
    //                         dist(randomEngine()),
    //                     });
    //                 }
    //                 // auto& b = s.addObject<BoxObject>(g, m);
    //                 auto& b = s.addObject<RoughBoxObject>(g, m);
    //                 // auto& b = s.addObject<RoughSphereObject>(g, m);
    //                 b.setTransformation(
    //                     Affine::Translation(Vec(
    //                         x + 0.2f * sdist(randomEngine()),
    //                         y + 0.2f * sdist(randomEngine()),
    //                         z + 0.2f * sdist(randomEngine())
    //                     ))
    //                     * Linear::RotationX(0.05f * sdist(randomEngine()))
    //                     * Linear::RotationY(0.05f * sdist(randomEngine()))
    //                     * Linear::RotationZ(0.05f * sdist(randomEngine()))
    //                 );
    //             }
    //         }
    //     }
    // }

    // Just one sphere (or box)
    // {
    //     const auto g = Sphere(4.0f);
    //     // const auto g = Rectangle(Vec(1.0f, 1.0f, 1.0f));
    //     // const auto g = 1.0f;
    //     auto m1 = BasicMaterial{};
    //     m1.setDiffuseReflection(1.0f);
    //     m1.setSpecularReflection(0.0f);
    //     m1.setSpecularSharpness(1.0f);
    //     m1.setTransmittance(0.0f);
    //     m1.setReflectedAbsorption(Color{1.0f, 1.0f, 1.0f});
    //     m1.setEmittedLuminance(Color{0.0f, 0.0f, 0.0f});

    //     auto m2 = BasicMaterial{};
    //     m2.setDiffuseReflection(0.05f);// * dist(randomEngine()));
    //     m2.setSpecularReflection(0.95f);//1.0f * dist(randomEngine()));
    //     m2.setSpecularSharpness(0.9f);//1.0f * dist(randomEngine()));
    //     m2.setTransmittance(0.0f);
    //     m2.setReflectedAbsorption(Color{1.0f, 0.9f, 0.2f});
    //     m2.setEmittedLuminance(Color{0.0f, 0.0f, 0.0f});

    //     auto& b = s.addObject<RoughSphereObject>(g, m1, m2);
    //     //auto& b = s.addObject<SphereObject>(g, m);
    //     // auto& b = s.addObject<BoxObject>(g, m);
    // }

    // // light source
    // {
    //     auto m = BasicMaterial{};
    //     m.setDiffuseReflection(1.0f);
    //     m.setSpecularReflection(0.0f);
    //     m.setTransmittance(0.0f);
    //     m.setReflectedAbsorption(Color{1.0f, 1.0f, 1.0f});
    //     m.setEmittedLuminance(Color{1.0f, 1.0f, 1.0f});

    //     auto& b = s.addObject<BoxObject>(Rectangle(Vec(100.0f, 0.1f, 100.0f)), m);
    //     b.setTransformation(Affine::Translation(Vec{0.0f, -50.0f, 0.0f}));
    // }

    // Fractal
    /*{
        auto o = std::make_unique<FractalObject>();
        o->material.setDiffuseReflection(1.0f);
        o->material.setSpecularReflection(0.3f);
        o->material.setSpecularSharpness(0.97f);
        o->material.setReflectedAbsorption(Color{1.0f, 1.0f, 1.0f});
        o->material.setEmittedLuminance(Color{0.01f, 0.03f, 0.1f});

        s.addObject(std::move(o));
    }*/

    // Test sphere
    /*
    {
        auto o = std::make_unique<SphereObject>(Sphere{Pos{-1.6f, -0.8f, 0.0f}, 0.2f});
        o->material.setDiffuseReflection(1.0f);
        o->material.setSpecularReflection(0.3f);
        o->material.setSpecularSharpness(0.97f);
        o->material.setReflectedAbsorption(Color{1.0f, 1.0f, 1.0f});
        o->material.setEmittedLuminance(Color{0.01f, 0.03f, 0.1f});
        s.addObject(std::move(o));
    }*/

    // Glass sphere
    /*{
        const auto geom = Sphere(
            Pos(0.0f, -1.0f, 0.0f),
            1.5f
        );
        auto o = std::make_unique<SphereObject>(geom);
        o->material.setDiffuseReflection(0.05f);
        o->material.setSpecularReflection(0.1f);
        o->material.setSpecularSharpness(1.0f);
        o->material.setTransmittance(0.9f);
        o->material.setIndexOfRefraction(1.55f);
        o->material.setReflectedAbsorption(Color{0.9f, 0.9f, 0.9f});
        o->material.setEmittedLuminance(Color{0.1f, 0.1f, 0.1f});

        s.addObject(std::move(o));
    }*/


    // Ground plane
    /*{
        auto o = std::make_unique<BoxObject>(
            Box{Pos{0.0f, 2.0f, 0.0f}, Vec{5.0f, 0.1f, 5.0f}}
        );
        o->material.setDiffuseReflection(1.0f);
        o->material.setSpecularReflection(0.0f);
        o->material.setSpecularSharpness(0.4f);
        o->material.setReflectedAbsorption(Color{1.0f, 1.0f, 1.0f});
        o->material.setTransmittance(0.0f);
        o->material.setEmittedLuminance(Color{0.0f, 0.0f, 0.0f});
        s.addObject(std::move(o));
    }*/


    s.updateGeometry();

    auto camera = Camera(Affine{});
    camera.setAspectRatio(1.0f);
    camera.setFieldOfView(10.0f);
    camera.setFocalDistance(30.0f);
    camera.setFocalBlurRadius(0.0f);//0.05f);
    const auto T = Affine::Translation(0.0f, 0.0f, -30.0f) * Linear::Scale(0.01f);
    camera.setTransform(T);

    auto renderSettings = RenderSettings(256, 256);
    renderSettings.setNumBounces(8);
    renderSettings.setSamplesPerPixel(1);

    auto r = Renderer{};
    r.startThreadPool();

    auto window = sf::RenderWindow(
        sf::VideoMode(
            static_cast<unsigned int>(renderSettings.width()),
            static_cast<unsigned int>(renderSettings.height())
        ),
        "Path Tracer"
    );

    auto img = sf::Image();
    auto tex = sf::Texture();
    tex.create(
        static_cast<unsigned int>(renderSettings.width()),
        static_cast<unsigned int>(renderSettings.height())
    );


    auto done = std::atomic<bool>{false};
    auto renderReset = std::atomic<bool>{true};
    auto textureMutex = std::mutex{};
    auto renderMutex = std::mutex{};
    auto cameraSettingsMutex = std::mutex{};

    auto renderThread = std::thread([&]{
        auto count = size_t{0};
        auto localCamera = Camera({});
        auto localSettings = RenderSettings{};
        {
            auto lock = std::lock_guard{cameraSettingsMutex};
            localCamera = camera;
            localSettings = renderSettings;
        }
        auto acc = Image(localSettings.width(), localSettings.height());
        while (!done.load()) {
            auto lock = std::lock_guard{renderMutex};
            if (renderReset.load()) {
                {
                    auto settingsLock = std::lock_guard{cameraSettingsMutex};
                    localCamera = camera;
                    localSettings = renderSettings;
                }
                {
                    auto texLock = std::lock_guard{textureMutex};
                    tex.create(
                        static_cast<unsigned>(localSettings.width()),
                        static_cast<unsigned>(localSettings.height())
                    );
                    acc = Image(localSettings.width(), localSettings.height());
                }
                count = 0;
                renderReset.store(false);
            }
            auto rendered = r.render(s, localCamera, localSettings);
            acc += rendered;
            count += 1;
            {
                auto texLock = std::lock_guard{textureMutex};
                auto scaled = scaleImage(acc, 1.0f / static_cast<float>(count));
                auto mapped = FilmicToneMapper{}(scaled);
                copyToSFImage(scaled, img);
                tex.loadFromImage(img);
            }
        }
    });

    auto running = true;
    while (running) {
        auto event = sf::Event{};
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) {
                running = false;
                r.cancelRender();
            } else if (event.type == sf::Event::Resized) {
                const auto widthF = static_cast<float>(event.size.width);
                const auto heightF = static_cast<float>(event.size.height);
                window.setView(sf::View(sf::FloatRect(
                    0.0f, 0.0f, widthF, heightF
                )));
                {
                    auto settingsLock = std::lock_guard{cameraSettingsMutex};
                    camera.setAspectRatio(
                        widthF / heightF
                    );
                    renderSettings.setSize(
                        event.size.width,
                        event.size.height
                    );
                }
                renderReset.store(true);
                r.cancelRender();
            } else if (event.type == sf::Event::KeyPressed) {
                auto delta = Vec(0.0f, 0.0f, 0.0f);
                switch (event.key.code) {
                    case sf::Keyboard::Key::A: delta.x -= 1.0f; break;
                    case sf::Keyboard::Key::D: delta.x += 1.0f; break;
                    case sf::Keyboard::Key::S: delta.z -= 1.0f; break;
                    case sf::Keyboard::Key::W: delta.z += 1.0f; break;
                    case sf::Keyboard::Key::Q: delta.y -= 1.0f; break;
                    case sf::Keyboard::Key::E: delta.y += 1.0f; break;
                    case sf::Keyboard::Key::Dash: {
                        {
                            auto settingsLock = std::lock_guard{cameraSettingsMutex};
                            if (event.key.control) {
                                camera.setFocalBlurRadius(
                                    std::max(0.0f, camera.focalBlurRadius() - 0.1f)
                                );
                            } else {
                                camera.setFieldOfView(camera.fieldOfView() + 1.0f);
                            }
                        }
                        renderReset.store(true);
                        r.cancelRender();
                        break;
                    }
                    case sf::Keyboard::Key::Equal: {
                        {
                            auto settingsLock = std::lock_guard{cameraSettingsMutex};
                            if (event.key.control) {
                                camera.setFocalBlurRadius(camera.focalBlurRadius() + 0.1f);
                            } else {
                                camera.setFieldOfView(camera.fieldOfView() - 1.0f);
                            }
                        }
                        renderReset.store(true);
                        r.cancelRender();
                        break;
                    }
                    case sf::Keyboard::Key::Enter: {
                        auto now = std::time(nullptr);
                        auto tm = *std::localtime(&now);
                        auto ss = std::ostringstream{};
                        // TODO: use std::filesystem and create the output
                        // directory if it doesn't already exist
                        // or change working directory in CMake
                        ss << "../output/image "
                            << std::put_time(&tm, "%Y-%m-%d %H-%M-%S")
                            << ".png";
                        const auto path = ss.str();
                        {
                            // TODO: avoid mutex contention here
                            auto lock = std::lock_guard{textureMutex};
                            img.saveToFile(path);
                        }
                        std::cout << "Saved image to \"" << path
                            << '\"' << std::endl;
                    }
                }
                if (delta.norm() > 1e-3f) {
                    auto settingsLock = std::lock_guard{cameraSettingsMutex};
                    auto& t = camera.transform();
                    if (event.key.control) {
                        camera.setFocalDistance(
                            camera.focalDistance() + 0.25f * delta.z
                        );
                    } else if (event.key.shift) {
                        const auto k = 3.141592654f * 0.025f;
                        t.linear *= (
                            Linear::RotationX(k * delta.x)
                            * Linear::RotationY(k * delta.y)
                            * Linear::RotationZ(k * delta.z)
                        );
                    } else {
                        t.translation += 10.0f * t.linear * delta;
                    }
                    renderReset.store(true);
                    r.cancelRender();
                }
            }
        }
        if (!running) {
            break;
        }

        window.clear();
        {
            auto lock = std::lock_guard{textureMutex};
            const auto spr = sf::Sprite(tex);
            window.draw(spr);
        }
        window.display();
    }

    done.store(true);
    renderThread.join();

    return 0;
}
