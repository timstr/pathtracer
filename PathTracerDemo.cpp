#include <PathTracer.hpp>

#include <SFML/Graphics.hpp>

int main() {

    auto s = Scene{};

    auto sphere0 = std::make_unique<SphereObject>(Sphere{Pos{ 0.0f,  2.0f, 10.0f}, 2.0f});
    auto sphere1 = std::make_unique<SphereObject>(Sphere{Pos{-0.6f, -0.6f, 5.0f}, 0.4f});
    auto sphere2 = std::make_unique<SphereObject>(Sphere{Pos{-0.6f,  0.6f, 5.0f}, 0.4f});
    auto sphere3 = std::make_unique<SphereObject>(Sphere{Pos{ 0.6f, -0.6f, 5.0f}, 0.4f});
    auto sphere4 = std::make_unique<SphereObject>(Sphere{Pos{ 0.6f,  0.6f, 5.0f}, 0.4f});
    auto sphere5 = std::make_unique<SphereObject>(Sphere{Pos{ 0.0f,  0.0f, 5.0f}, 0.2f});
    auto sphere6 = std::make_unique<SphereObject>(Sphere{Pos{-0.8f,  0.0f, 2.0f}, 0.1f});
    auto sphere7 = std::make_unique<SphereObject>(Sphere{Pos{ 0.8f,  0.0f, 4.0f}, 0.1f});
    auto sphere8 = std::make_unique<SphereObject>(Sphere{Pos{ 0.0f, -0.8f, 6.0f}, 0.1f});
    auto sphere9 = std::make_unique<SphereObject>(Sphere{Pos{ 0.0f,  0.8f, 8.0f}, 0.1f});

    sphere0->material.reflectedColor = Color{1.0f, 1.0f, 1.0f};
    sphere1->material.reflectedColor = Color{1.0f, 0.3f, 0.3f};
    sphere2->material.reflectedColor = Color{0.0f, 1.0f, 0.0f};
    sphere3->material.reflectedColor = Color{0.0f, 0.0f, 1.0f};
    sphere4->material.reflectedColor = Color{1.0f, 1.0f, 0.0f};
    sphere5->material.reflectedColor = Color{1.0f, 1.0f, 1.0f};
    sphere6->material.reflectedColor = Color{1.0f, 1.0f, 1.0f};
    sphere7->material.reflectedColor = Color{1.0f, 1.0f, 1.0f};
    sphere8->material.reflectedColor = Color{1.0f, 1.0f, 1.0f};
    sphere9->material.reflectedColor = Color{1.0f, 1.0f, 1.0f};

    sphere0->material.emittedRadiance = Color{0.0f, 0.0f, 0.0f};
    sphere1->material.emittedRadiance = Color{0.1f, 0.0f, 0.0f};
    sphere2->material.emittedRadiance = Color{0.0f, 0.1f, 0.0f};
    sphere3->material.emittedRadiance = Color{0.0f, 0.0f, 0.1f};
    sphere4->material.emittedRadiance = Color{0.1f, 0.1f, 0.0f};
    sphere5->material.emittedRadiance = Color{5.0f, 5.0f, 5.0f};
    sphere6->material.emittedRadiance = Color{3.0f, 1.0f, 1.0f};
    sphere7->material.emittedRadiance = Color{1.0f, 3.0f, 1.0f};
    sphere8->material.emittedRadiance = Color{1.0f, 1.0f, 3.0f};
    sphere9->material.emittedRadiance = Color{3.0f, 3.0f, 1.0f};

    sphere0->material.diffuseness = 1.0f;
    sphere1->material.diffuseness = 0.05f;
    sphere2->material.diffuseness = 1.0f / 3.0f;
    sphere3->material.diffuseness = 2.0f / 3.0f;
    sphere4->material.diffuseness = 1.0f;
    sphere5->material.diffuseness = 1.0f;
    sphere6->material.diffuseness = 0.5f;
    sphere7->material.diffuseness = 0.5f;
    sphere8->material.diffuseness = 0.5f;
    sphere9->material.diffuseness = 0.5f;

    s.addObject(std::move(sphere0));
    s.addObject(std::move(sphere1));
    s.addObject(std::move(sphere2));
    s.addObject(std::move(sphere3));
    s.addObject(std::move(sphere4));
    s.addObject(std::move(sphere5));
    s.addObject(std::move(sphere6));
    s.addObject(std::move(sphere7));
    s.addObject(std::move(sphere8));
    s.addObject(std::move(sphere9));

    auto c = OrthographicCamera({}, 1.0f);

    auto r = Renderer{};

    const auto w = std::size_t{1024};
    const auto h = std::size_t{1024};
    const auto b = std::size_t{16};
    const auto spp = std::size_t{1024};

    auto rendered = r.render(s, c, h, w, b, spp);

    // auto tm = ReinhardToneMapper();
    // auto toneMapped = tm(rendered);

    const auto& finalImage = rendered; // toneMapped;

    auto img = sf::Image();
    img.create(w, h);
    for (unsigned x = 0; x < w; ++x) {
        for (unsigned y = 0; y < h; ++y) {
            auto color = finalImage(x, y);
            auto px = sf::Color(
                static_cast<std::uint8_t>(std::clamp(color.r, 0.0f, 1.0f) * 255.0f),
                static_cast<std::uint8_t>(std::clamp(color.g, 0.0f, 1.0f) * 255.0f),
                static_cast<std::uint8_t>(std::clamp(color.b, 0.0f, 1.0f) * 255.0f)
            );
            img.setPixel(x, y, px);
        }
    }
    img.saveToFile("output.png");
    return 0;
}