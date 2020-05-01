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
    auto pDist = std::uniform_real_distribution<float>{0.0f, 1.0f};
    auto pnDist = std::uniform_real_distribution<float>{-1.0f, 1.0f};
    for (std::size_t i = 0; i < 500; ++i) {
        const auto x = pDist(randomEngine());
        const auto y = pDist(randomEngine());
        const auto geom = Sphere(
            Pos(3.5f * (2.0f * x - 1.0f), 3.5f * (2.0f * y - 1.0f), 1.0f * pnDist(randomEngine())),
            0.5f * pDist(randomEngine())
        );
        auto o = std::make_unique<SphereObject>(geom);
        o->material.setDiffuseReflection(pDist(randomEngine()));
        o->material.setSpecularReflection(pDist(randomEngine()));
        o->material.setSpecularSharpness(pDist(randomEngine()));
        o->material.setTransmittance(pDist(randomEngine()));
        o->material.setIndexOfRefraction(1.0f + 0.8f * pDist(randomEngine()));
        o->material.setReflectedAbsorption(Color{
            0.3f + 0.7f * x,
            0.3f + 0.7f * y,
            0.3f + 0.7f * pDist(randomEngine())
        });
        o->material.setEmittedLuminance(Color{
            0.2f * x,
            0.2f * y,
            0.1f * pDist(randomEngine())
        });

        s.addObject(std::move(o));
    }

    // one sphere
    {
        const auto geom = Sphere(
            Pos(0.0f, 0.1f, 0.0f),
            0.1f
        );
        auto o = std::make_unique<SphereObject>(geom);
        o->material.setDiffuseReflection(1.0f);
        o->material.setSpecularReflection(1.0f);
        o->material.setSpecularSharpness(0.8f);
        o->material.setTransmittance(0.0f);
        o->material.setReflectedAbsorption(Color{0.2f, 0.4f, 0.9f});
        o->material.setEmittedLuminance(Color{0.0f, 0.7f, 0.2f});

        s.addObject(std::move(o));
    }

    // light source
    {
        const auto geom = Sphere(
            Pos(0.0f, -100.0f, 0.0f),
            40.0f
        );
        auto o = std::make_unique<SphereObject>(geom);
        o->material.setDiffuseReflection(1.0f);
        o->material.setSpecularReflection(0.0f);
        o->material.setTransmittance(0.0f);
        o->material.setReflectedAbsorption(Color{1.0f, 1.0f, 1.0f});
        o->material.setEmittedLuminance(Color{30.0f, 20.0f, 10.0f});

        s.addObject(std::move(o));
    }


    s.updateGeometry();

    auto c = PerspectiveCamera(Affine{});
    c.setAspectRatio(1.0f);
    c.setFieldOfView(5.0f);
    c.setFocalDistance(50.0f);
    c.setFocalBlurRadius(0.0f);// 0.1f);
    const auto translation = Affine::Translation(0.0f, 0.0f, -500.0f);
    const auto scale = Linear::RotationX(3.141592654f * 0.0f) * Linear::Scale(0.1f);
    const auto T = scale * translation;
    c.setTransform(T);

    auto r = Renderer{1024, 1024};
    r.setNumBounces(32);
    r.setSamplesPerPixel(1024);

    std::cout << "Using " << r.numThreads() << " threads\n";

    {
        auto lrr = r;
        lrr.setWidth(r.width() / 8);
        lrr.setHeight(r.height() / 8);
        lrr.setSamplesPerPixel(std::max(r.samplesPerPixel() / 2, std::size_t{1}));
        std::cout << "Low resolution version\n";

        auto rendered = lrr.render(s, c);
        saveImage(rendered, "output low res.png");
    }
    std::cout << "High resolution version\n";

    {
        auto rendered = r.render(s, c);

        // auto tm = ReinhardToneMapper();
        // auto toneMapped = tm(rendered);

        saveImage(rendered, "output before.png");
    }


    // Glass sphere
    {
        const auto geom = Sphere(
            Pos(0.0f, 0.0f, -49.5f),
            0.1f
        );
        auto o = std::make_unique<SphereObject>(geom);
        o->material.setDiffuseReflection(0.05f);
        o->material.setSpecularReflection(0.1f);
        o->material.setSpecularSharpness(1.0f);
        o->material.setTransmittance(0.9f);
        o->material.setIndexOfRefraction(1.2f);
        o->material.setReflectedAbsorption(Color{0.9f, 0.9f, 0.9f});
        o->material.setEmittedLuminance(Color{0.1f, 0.1f, 0.1f});

        s.addObject(std::move(o));
        s.updateGeometry();
    }

    {
        auto rendered = r.render(s, c);

        // auto tm = ReinhardToneMapper();
        // auto toneMapped = tm(rendered);

        saveImage(rendered, "output after.png");
    }


    return 0;
}