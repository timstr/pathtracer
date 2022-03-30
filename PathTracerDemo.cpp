#include <PathTracer.hpp>
#include <RandomNumberGenerator.hpp>

#include <SFML/Graphics.hpp>

#include <iostream>

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
            ray.dir
        };
    }
};

void saveImage(const Image& rendered, const std::string& path) {
    auto img = sf::Image();
    img.create(
        static_cast<unsigned int>(rendered.width()),
        static_cast<unsigned int>(rendered.height())
    );
    for (unsigned x = 0; x < rendered.width(); ++x) {
        for (unsigned y = 0; y < rendered.height(); ++y) {
            auto color = rendered(x, y);
            auto px = sf::Color(
                static_cast<std::uint8_t>(std::clamp(color.r, 0.0f, 1.0f) * 255.0f),
                static_cast<std::uint8_t>(std::clamp(color.g, 0.0f, 1.0f) * 255.0f),
                static_cast<std::uint8_t>(std::clamp(color.b, 0.0f, 1.0f) * 255.0f)
            );
            img.setPixel(x, y, px);
        }
    }
    img.saveToFile(path);
}

int main() {

    auto s = Scene{};

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

        auto& lightSink = s.addObject<BoxObject>(Rectangle(Vec{1.0f, d, 1.0f}), antiGlowing);
        lightSink.setTransformation(Affine::Translation(Vec{0.0f, 5.0f - 2.0f * d, 0.0f}));

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

    // Globe of cubes
    /*
    {
        const auto res = std::size_t{32};
        const auto size = 3.0f;
        const auto dsize = size / static_cast<float>(res);
        const auto mapToSpace = [&](std::size_t idx) {
            return (static_cast<float>(idx) / static_cast<float>(res - 1) * 2.0f - 1.0f) * size;
        };


        const auto dist = std::uniform_real_distribution<float>{0.0f, 1.0f};

        for (std::size_t i = 0; i < res; ++i){
            const auto x = mapToSpace(i);
            for (std::size_t j = 0; j < res; ++j){
                const auto y = mapToSpace(j);
                for (std::size_t k = 0; k < res; ++k){
                    const auto z = mapToSpace(k);
                    if (x * x + y * y + z * z > size) {
                        continue;
                    }

                    if (dist(randomEngine()) < 0.1f) {
                        continue;
                    }

                    auto b = Box{Pos{}, Vec{dsize, dsize, dsize}};
                    b.halfSize.x *= (1.0f + 1.0f * dist(randomEngine()));
                    b.halfSize.y *= (1.0f + 1.0f * dist(randomEngine()));
                    b.halfSize.z *= (1.0f + 1.0f * dist(randomEngine()));
                    b.center = Pos{ x, y, z };

                    auto m = BasicMaterial{};
                    if (dist(randomEngine()) < 0.7f) {
                        m.setDiffuseReflection(1.0f);
                        m.setSpecularReflection(0.05f);
                        m.setSpecularSharpness(0.9f);
                        m.setTransmittance(0.0f);
                        m.setReflectedAbsorption(Color{1.0f, 1.0f, 1.0f});
                        m.setEmittedLuminance(Color{0.0f, 0.0f, 0.0f});
                    } else {
                        m.setDiffuseReflection(0.02f);
                        m.setSpecularReflection(0.2f);
                        m.setSpecularSharpness(0.99f);
                        m.setTransmittance(0.95f);
                        m.setIndexOfRefraction(1.55f);
                        m.setReflectedAbsorption(Color{0.8f, 0.9f, 1.0f});
                        m.setEmittedLuminance(Color{0.01f, 0.02f, 0.04f});
                    }
                    s.addObject(std::make_unique<BoxObject>(b, m));
                }
            }
        }
    }
    */

    // light source
    /*{
        auto m = BasicMaterial{};
        m.setDiffuseReflection(1.0f);
        m.setSpecularReflection(0.0f);
        m.setTransmittance(0.0f);
        m.setReflectedAbsorption(Color{1.0f, 1.0f, 1.0f});
        m.setEmittedLuminance(Color{1.0f, 1.0f, 1.0f});

        s.addObject(std::make_unique<BoxObject>(Box{Pos{0.0f, -5.0f, 0.0f}, Vec{10.0f, 0.1f, 10.0f}}, m));
    }*/

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

    auto c = PerspectiveCamera(Affine{});
    c.setAspectRatio(1.0f);
    c.setFieldOfView(10.0f);
    c.setFocalDistance(30.0f);
    c.setFocalBlurRadius(0.05f);
    const auto T = Affine::Translation(0.0f, 0.0f, -30.0f) * Linear::Scale(1.0);
    c.setTransform(T);

    auto r = Renderer{512, 512};
    // r.setNumThreads(1);
    r.setNumBounces(16);
    r.setSamplesPerPixel(1024);

    std::cout << "Using " << r.numThreads() << " threads\n";

    /*
    {
        auto lrr = r;
        lrr.setWidth(r.width() / 8);
        lrr.setHeight(r.height() / 8);
        lrr.setSamplesPerPixel(std::max(r.samplesPerPixel(), std::size_t{1}));
        std::cout << "Low resolution version\n";

        auto rendered = lrr.render(s, c);
        saveImage(rendered, "output low res.png");
    }
    std::cout << "High resolution version\n";
    */

    {
        auto rendered = r.render(s, c);

        //auto tm = ReinhardToneMapper();
        //auto toneMapped = tm(rendered);

        saveImage(rendered, "output.png");
    }


    return 0;
}