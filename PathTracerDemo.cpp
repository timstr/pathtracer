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

    std::optional<float> hit(const Ray&) const noexcept override {
        return 1e10f;
    }

    ColorBounce deflect(const Ray& ray) const noexcept override {
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

    // Random spheres
    /*
    {
        const auto pdist = std::uniform_real_distribution<float>{0.0f, 1.0f};
        const auto sdist = std::uniform_real_distribution<float>{-1.0f, 1.0f};
        for (std::size_t i = 0; i < 30; ++i) {
            const auto geom = Sphere{
                Pos{
                    5.0f * sdist(randomEngine()),
                    -4.0f - 1.0f * sdist(randomEngine()),
                    5.0f * sdist(randomEngine())
                },
                0.1f + 0.4f * sdist(randomEngine())
            };
            auto o = std::make_unique<SphereObject>(geom);
            o->material.setDiffuseReflection(1.0f);
            o->material.setSpecularReflection(0.0f);
            o->material.setTransmittance(0.0f);
            o->material.setEmittedLuminance(Color{
                5.0f * pdist(randomEngine()),
                5.0f * pdist(randomEngine()),
                5.0f * pdist(randomEngine())
            });
            s.addObject(std::move(o));
        }
    }
    */

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
    {
        auto m = BasicMaterial{};
        m.setDiffuseReflection(1.0f);
        m.setSpecularReflection(0.0f);
        m.setTransmittance(0.0f);
        m.setReflectedAbsorption(Color{1.0f, 1.0f, 1.0f});
        m.setEmittedLuminance(Color{1.0f, 1.0f, 1.0f});

        s.addObject(std::make_unique<BoxObject>(Box{Pos{0.0f, -2.0f, 0.0f}, Vec{1.0f, 0.1f, 1.0f}}, m));
    }

    // Fractal
    {
        auto o = std::make_unique<FractalObject>();
        o->material.setDiffuseReflection(1.0f);
        o->material.setSpecularReflection(0.3f);
        o->material.setSpecularSharpness(0.97f);
        o->material.setReflectedAbsorption(Color{1.0f, 1.0f, 1.0f});
        o->material.setEmittedLuminance(Color{0.01f, 0.03f, 0.1f});
        
        s.addObject(std::move(o));
    }

    // Test sphere
    {
        auto o = std::make_unique<SphereObject>(Sphere{Pos{-1.6f, -0.8f, 0.0f}, 0.2f});
        o->material.setDiffuseReflection(1.0f);
        o->material.setSpecularReflection(0.3f);
        o->material.setSpecularSharpness(0.97f);
        o->material.setReflectedAbsorption(Color{1.0f, 1.0f, 1.0f});
        o->material.setEmittedLuminance(Color{0.01f, 0.03f, 0.1f});
        s.addObject(std::move(o));
    }

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
    c.setFocalDistance(5.0f);
    c.setFocalBlurRadius(0.0f);
    const auto translation = Affine::Translation(0.0f, 0.0f, -5.0f);
    const auto scale = Linear::RotationX(3.141592654f * -0.1f) * Linear::RotationY(0.4f);
    const auto T = scale * translation;
    c.setTransform(T);

    auto r = Renderer{256, 256};
    // r.setNumThreads(1);
    r.setNumBounces(16);
    r.setSamplesPerPixel(8);

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