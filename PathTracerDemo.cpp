#include <PathTracer.hpp>

#include <SFML/Graphics.hpp>

int main() {

    auto s = Scene{};

    auto sphere1 = std::make_unique<SphereObject>(Sphere{Pos{0.3f, 0.3f, 5.0f}, 0.2f});
    auto sphere2 = std::make_unique<SphereObject>(Sphere{Pos{0.3f, 0.7f, 5.0f}, 0.2f});
    auto sphere3 = std::make_unique<SphereObject>(Sphere{Pos{0.7f, 0.3f, 5.0f}, 0.2f});
    auto sphere4 = std::make_unique<SphereObject>(Sphere{Pos{0.7f, 0.7f, 5.0f}, 0.2f});
    auto sphere5 = std::make_unique<SphereObject>(Sphere{Pos{0.5f, 0.5f, 5.0f}, 0.1f});

    sphere1->material.diffuseness = 0.0f;
    sphere2->material.diffuseness = 1.0f / 3.0f;
    sphere3->material.diffuseness = 2.0f / 3.0f;
    sphere4->material.diffuseness = 1.0f;
    sphere5->material.diffuseness = 1.0f;

    sphere5->material.emittedRadiance = Color{100.0f, 100.0f, 100.0f};

    s.addObject(std::move(sphere1));
    s.addObject(std::move(sphere2));
    s.addObject(std::move(sphere3));
    s.addObject(std::move(sphere4));
    s.addObject(std::move(sphere5));

    auto c = OrthographicCamera({}, 1.0f);

    auto r = Renderer{};

    const auto w = std::size_t{512};
    const auto h = std::size_t{512};

    const auto spp = std::size_t{4};

    auto rendered = r.render(s, c, h, w, spp);

    auto tm = ReinhardToneMapper();

    auto toneMapped = tm(rendered);

    const auto& finalImage = rendered;

    auto img = sf::Image();
    img.create(w, h);
    for (std::size_t x = 0; x < w; ++x) {
        for (std::size_t y = 0; y < h; ++y) {
            auto color = finalImage(x, y);
            auto px = sf::Color(
                static_cast<std::uint8_t>(color.r * 255.0f),
                static_cast<std::uint8_t>(color.g * 255.0f),
                static_cast<std::uint8_t>(color.b * 255.0f)
            );
            img.setPixel(x, y, px);
        }
    }
    img.saveToFile("output.png");
    return 0;
}