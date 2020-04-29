#include <PathTracer.hpp>
#include <RandomNumberGenerator.hpp>

#include <SFML/Graphics.hpp>

#include <iostream>

class GlassSphere : public Object {
public:
    Sphere geometry;

    GlassSphere(Sphere _geometry) noexcept : geometry(_geometry) {}

    std::optional<float> hit(const Ray& ray) const noexcept override {
        return intersect(ray, geometry);
    }

    ColorBounce deflect(const Ray& ray) const noexcept override {
        // TODO aaaaaaaa
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

    //doTest();
    //return 0;

    auto s = Scene{};

    /*
    // just a triangle
    {
        const auto triangle = Triangle{
            Pos{ 0.0f, -0.45f,  5.0f},
            Pos{-0.4f,  0.23f, 5.0f},
            Pos{ 0.4f,  0.23f, 5.0f}
        };
        auto o = std::make_unique<TriangleObject>(triangle);
        o->material.reflectedColor = Color{0.2f, 0.5f, 0.2f};
        o->material.emittedRadiance = Color{0.0f, 0.1f, 0.0f};
        o->material.diffuseness = 0.8f;

        s.addObject(std::move(o));
    }

    // Big shiny spheres in the back
    {
        auto o1 = std::make_unique<SphereObject>(Sphere{Pos{-1.0f,  0.0f,  8.0f}, 0.9f});
        auto o2 = std::make_unique<SphereObject>(Sphere{Pos{ 1.0f,  0.0f,  8.0f}, 0.9f});
        auto o3 = std::make_unique<SphereObject>(Sphere{Pos{ 0.0f, -1.0f,  9.0f}, 0.9f});
        auto o4 = std::make_unique<SphereObject>(Sphere{Pos{ 0.0f,  1.0f,  9.0f}, 0.9f});
        o1->material.reflectedColor = Color{1.0f, 1.0f, 1.0f};
        o2->material.reflectedColor = Color{1.0f, 1.0f, 1.0f};
        o3->material.reflectedColor = Color{1.0f, 1.0f, 1.0f};
        o4->material.reflectedColor = Color{1.0f, 1.0f, 1.0f};
        o1->material.emittedRadiance = Color{0.1f, 0.1f, 0.0f };
        o2->material.emittedRadiance = Color{0.0f, 0.1f, 0.1f};
        o3->material.emittedRadiance = Color{0.1f, 0.1f, 0.0f };
        o4->material.emittedRadiance = Color{0.1f, 0.1f, 0.1f};
        o1->material.diffuseness = 0.0f;
        o2->material.diffuseness = 0.0f;
        o3->material.diffuseness = 0.0f;
        o4->material.diffuseness = 0.0f;
        s.addObject(std::move(o1));
        s.addObject(std::move(o2));
        s.addObject(std::move(o3));
        s.addObject(std::move(o4));
    }

    // Coloured spheres
    {
        auto o1 = std::make_unique<SphereObject>(Sphere{Pos{-0.6f, -0.6f, 5.0f}, 0.4f});
        auto o2 = std::make_unique<SphereObject>(Sphere{Pos{-0.6f,  0.6f, 5.0f}, 0.4f});
        auto o3 = std::make_unique<SphereObject>(Sphere{Pos{ 0.6f, -0.6f, 5.0f}, 0.4f});
        auto o4 = std::make_unique<SphereObject>(Sphere{Pos{ 0.6f,  0.6f, 5.0f}, 0.4f});
        o1->material.reflectedColor = Color{1.0f, 0.5f, 0.5f};
        o2->material.reflectedColor = Color{0.5f, 1.0f, 0.5f};
        o3->material.reflectedColor = Color{0.5f, 0.5f, 1.0f};
        o4->material.reflectedColor = Color{1.0f, 1.0f, 0.5f};
        o1->material.emittedRadiance = Color{0.1f, 0.0f, 0.0f};
        o2->material.emittedRadiance = Color{0.0f, 0.1f, 0.0f};
        o3->material.emittedRadiance = Color{0.0f, 0.0f, 0.1f};
        o4->material.emittedRadiance = Color{0.1f, 0.1f, 0.0f};
        o1->material.diffuseness = 0.05f;
        o2->material.diffuseness = 1.0f / 3.0f;
        o3->material.diffuseness = 2.0f / 3.0f;
        o4->material.diffuseness = 1.0f;
        s.addObject(std::move(o1));
        s.addObject(std::move(o2));
        s.addObject(std::move(o3));
        s.addObject(std::move(o4));
    }

    // Glowing sphere in center
    {
        auto o = std::make_unique<SphereObject>(Sphere{Pos{ 0.0f,  0.0f, 5.0f}, 0.2f});
        o->material.reflectedColor = Color{1.0f, 1.0f, 1.0f};
        o->material.emittedRadiance = Color{10.0f, 10.0f, 10.0f};
        o->material.diffuseness = 1.0f;
        s.addObject(std::move(o));
    }

    // smaller staggered glowing spheres between coloured spheres
    {
        auto o1 = std::make_unique<SphereObject>(Sphere{Pos{-0.8f,  0.0f, 2.0f}, 0.1f});
        auto o2 = std::make_unique<SphereObject>(Sphere{Pos{ 0.8f,  0.0f, 4.0f}, 0.1f});
        auto o3 = std::make_unique<SphereObject>(Sphere{Pos{ 0.0f, -0.8f, 6.0f}, 0.1f});
        auto o4 = std::make_unique<SphereObject>(Sphere{Pos{ 0.0f,  0.8f, 8.0f}, 0.1f});
        o1->material.reflectedColor = Color{1.0f, 1.0f, 1.0f};
        o2->material.reflectedColor = Color{1.0f, 1.0f, 1.0f};
        o3->material.reflectedColor = Color{1.0f, 1.0f, 1.0f};
        o4->material.reflectedColor = Color{1.0f, 1.0f, 1.0f};
        o1->material.emittedRadiance = Color{3.0f, 1.0f, 1.0f};
        o2->material.emittedRadiance = Color{1.0f, 3.0f, 1.0f};
        o3->material.emittedRadiance = Color{1.0f, 1.0f, 3.0f};
        o4->material.emittedRadiance = Color{3.0f, 3.0f, 1.0f};
        o1->material.diffuseness = 0.5f;
        o2->material.diffuseness = 0.5f;
        o3->material.diffuseness = 0.5f;
        o4->material.diffuseness = 0.5f;
        s.addObject(std::move(o1));
        s.addObject(std::move(o2));
        s.addObject(std::move(o3));
        s.addObject(std::move(o4));
    }

    // tiny coloured glowing spheres in corners
    {
        auto o1 = std::make_unique<SphereObject>(Sphere{Pos{-0.95f, -0.95f, 4.5f}, 0.02f});
        auto o2 = std::make_unique<SphereObject>(Sphere{Pos{-0.95f,  0.95f, 4.5f}, 0.02f});
        auto o3 = std::make_unique<SphereObject>(Sphere{Pos{ 0.95f, -0.95f, 4.5f}, 0.02f});
        auto o4 = std::make_unique<SphereObject>(Sphere{Pos{ 0.95f,  0.95f, 4.5f}, 0.02f});
        o1->material.reflectedColor = Color{0.0f, 0.0f, 0.0f};
        o2->material.reflectedColor = Color{0.0f, 0.0f, 0.0f};
        o3->material.reflectedColor = Color{0.0f, 0.0f, 0.0f};
        o4->material.reflectedColor = Color{0.0f, 0.0f, 0.0f};
        o1->material.emittedRadiance = Color{  0.0f,   0.0f, 100.0f};
        o2->material.emittedRadiance = Color{100.0f, 100.0f,   0.0f};
        o3->material.emittedRadiance = Color{100.0f,   0.0f,   0.0f};
        o4->material.emittedRadiance = Color{  0.0f, 100.0f,   0.0f};
        o1->material.diffuseness = 1.0f;
        o2->material.diffuseness = 1.0f;
        o3->material.diffuseness = 1.0f;
        o4->material.diffuseness = 1.0f;
        s.addObject(std::move(o1));
        s.addObject(std::move(o2));
        s.addObject(std::move(o3));
        s.addObject(std::move(o4));
    }
    */

    /*auto pDist = std::uniform_real_distribution<float>{0.0f, 1.0f};
    auto pnDist = std::uniform_real_distribution<float>{-1.0f, 1.0f};
    for (std::size_t i = 0; i < 500; ++i) {
        const auto geom = Sphere(
            Pos(pnDist(randomEngine()), pnDist(randomEngine()), pnDist(randomEngine())),
            0.1f * pDist(randomEngine())
        );
        auto o = std::make_unique<SphereObject>(geom);
        o->material.diffuseness = pDist(randomEngine());
        o->material.emittedRadiance.r = pDist(randomEngine());
        o->material.emittedRadiance.g = pDist(randomEngine());
        o->material.emittedRadiance.b = pDist(randomEngine());
        o->material.reflectedColor.r = pDist(randomEngine());
        o->material.reflectedColor.g = pDist(randomEngine());
        o->material.reflectedColor.b = pDist(randomEngine());

        s.addObject(std::move(o));
    }*/

    // Ground plane triangles
    {
        auto y = 0.0f;
        auto bl = Pos{ -10.0f, y,  10.0f, };
        auto br = Pos{  10.0f, y,  10.0f, };
        auto fl = Pos{ -10.0f, y, -10.0f, };
        auto fr = Pos{  10.0f, y, -10.0f, };
        auto cc = Pos{  0.0f, y,    0.0f, };
        auto o1 = std::make_unique<TriangleObject>(Triangle{br, bl, cc});
        auto o2 = std::make_unique<TriangleObject>(Triangle{bl, fl, cc});
        auto o3 = std::make_unique<TriangleObject>(Triangle{fl, fr, cc});
        auto o4 = std::make_unique<TriangleObject>(Triangle{fr, br, cc});
        o1->material.diffuseness = 0.7f;
        o2->material.diffuseness = 0.7f;
        o3->material.diffuseness = 0.7f;
        o4->material.diffuseness = 0.7f;
        o1->material.emittedRadiance = Color{ 0.0f, 0.0f, 0.0f };
        o2->material.emittedRadiance = Color{ 0.0f, 0.0f, 0.0f };
        o3->material.emittedRadiance = Color{ 0.0f, 0.0f, 0.0f };
        o4->material.emittedRadiance = Color{ 0.0f, 0.0f, 0.0f };
        o1->material.reflectedColor = Color{ 1.0f, 0.3f, 0.4f };
        o2->material.reflectedColor = Color{ 1.0f, 1.0f, 0.4f };
        o3->material.reflectedColor = Color{ 1.0f, 0.3f, 0.4f };
        o4->material.reflectedColor = Color{ 1.0f, 1.0f, 0.4f };
        s.addObject(std::move(o1));
        s.addObject(std::move(o2));
        s.addObject(std::move(o3));
        s.addObject(std::move(o4));
    }
    // Light source sphere
    {
        auto o = std::make_unique<SphereObject>(Sphere{Pos{0.0f, -100.0f, 0.0f}, 50.0f});
        o->material.diffuseness = 0.01f;
        o->material.emittedRadiance = Color{3.0f, 3.0f, 3.0f};
        o->material.reflectedColor = Color{0.8f, 0.8f, 0.8f};
        s.addObject(std::move(o));
    }

    // Bouncing balls
    const auto makeBall = [&](Pos p, float r, BasicGlossyMaterial mtl) {
        auto o = std::make_unique<SphereObject>(Sphere{p, r}, mtl);
        auto ret = o.get();
        s.addObject(std::move(o));
        return ret;
    };

    makeBall(Pos{-6.3f, 0.0f, 0.0f}, 4.0f, BasicGlossyMaterial{0.05f, Color{0.9f, 0.5f, 0.5f}, Color{0.3f, 0.2f, 0.1f}});
    makeBall(Pos{ 0.0f, 0.0f, 2.0f}, 1.0f, BasicGlossyMaterial{0.05f, Color{0.5f, 0.5f, 0.9f}, Color{0.3f, 0.2f, 0.1f}});
    makeBall(Pos{ 6.3f, 0.0f, 0.0f}, 4.0f, BasicGlossyMaterial{0.05f, Color{0.9f, 0.5f, 0.5f}, Color{0.3f, 0.2f, 0.1f}});

    auto balls = std::vector<SphereObject*>{
        makeBall(Pos{-1.6f, -0.58f, 0.0f}, 0.55f, BasicGlossyMaterial{0.05f, Color{0.6f, 0.4f, 0.5f}, Color{0.4f, 0.4f, 0.6f}}),
        makeBall(Pos{-0.7f, -0.45f, 0.0f}, 0.25f, BasicGlossyMaterial{0.40f, Color{0.2f, 0.8f, 0.7f}, Color{0.2f, 0.3f, 0.0f}}),
        makeBall(Pos{ 0.0f, -0.35f, 0.0f}, 0.30f, BasicGlossyMaterial{0.20f, Color{0.8f, 0.2f, 0.4f}, Color{0.3f, 0.2f, 0.0f}}),
        makeBall(Pos{ 0.6f, -0.40f, 0.0f}, 0.20f, BasicGlossyMaterial{0.01f, Color{0.6f, 0.4f, 0.5f}, Color{0.0f, 0.2f, 0.2f}}),
        makeBall(Pos{ 1.3f, -0.60f, 0.0f}, 0.40f, BasicGlossyMaterial{0.6f, Color{0.6f, 0.4f, 0.5f}, Color{0.0f, 0.2f, 0.2f}})
    };

    auto phases = std::vector<float>{
        0.0476f,
        0.23721f,
        0.89452f,
        0.6372f,
        0.4376458f
    };

    auto frequencies = std::vector<float>{
        1.0f,
        3.0f,
        2.0f,
        4.0f,
        1.0f
    };

    /*{
        auto o = std::make_unique<SphereObject>(Sphere{Pos{0.0f, -0.1f, 0.0f}, 0.25f});
        o->material.diffuseness = 0.2f;
        o->material.emittedRadiance = Color{0.3f, 0.2f, 0.0f};
        o->material.reflectedColor = Color{0.8f, 0.2f, 0.3f};
        s.addObject(std::move(o));
    }*/

    auto c = PerspectiveCamera(Affine{});
    c.setAspectRatio(4.0f / 3.0f);
    c.setFieldOfView(15.0f);
    c.setFocalDistance(10.0f);
    c.setFocalBlurRadius(0.2f);
    // c.setTransform(Affine::Translation(0.0f, 0.0f, -10.0f) * Linear::RotationX(3.141592654f * 0.25f));
    c.transform().linear = Linear::RotationX(3.141592654f * 0.025f);
    c.transform().translation = Vec(0.0f, -1.8f, -10.0f);

    auto r = Renderer{800, 600};
    r.setNumBounces(32);
    r.setSamplesPerPixel(4096);

    std::cout << "Using " << r.numThreads() << " threads\n";

    std::size_t numFrames = 1;
    const auto baseHeight = 2.5f;
    for (std::size_t i = 0; i < numFrames; ++i) {
        const auto t = static_cast<float>(i) / static_cast<float>(numFrames);
        for (std::size_t j = 0; j < balls.size(); ++j) {
            const auto f = frequencies[j];
            const auto x0 = phases[j] + t * f;
            const auto x1 = x0 - std::floor(x0);
            const auto x2 = 4.0f * (-x1 * x1 + x1);
            const auto y = -x2 * baseHeight / f;
            balls[j]->geometry.center.y = y - balls[j]->geometry.radius;
        }

        std::cout << "Rendering frame " << std::to_string(i + 1) << " of " << std::to_string(numFrames) << '\n';

        auto rendered = r.render(s, c);

        // auto tm = ReinhardToneMapper();
        // auto toneMapped = tm(rendered);

        saveImage(rendered, "frame " + std::to_string(i) + ".png");
    }


    return 0;
}