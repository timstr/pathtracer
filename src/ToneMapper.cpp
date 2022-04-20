#include <ToneMapper.hpp>

#include <algorithm>
#include <cmath>

Image ReinhardToneMapper::operator()(const Image& img) const noexcept {
    // https://www.researchgate.net/publication/255682296_Parameter_Estimation_for_Photographic_Tone_Reproduction
    auto acc = 0.0f;
    const auto w = img.width();
    const auto h = img.height();
    auto minLum = std::numeric_limits<float>::max();
    auto maxLum = std::numeric_limits<float>::lowest();
    for (std::size_t i = 0; i < w; ++i) {
        for (std::size_t j = 0; j < h; ++j) {
            const auto& c = img(i, j);
            const auto lum = 0.27f * c.r + 0.67f * c.b + 0.06f * c.b;
            minLum = std::min(minLum, lum);
            maxLum = std::max(maxLum, lum);
            acc += std::log(lum + 1e-6f);
        }
    }

    const auto avgLum = std::exp(acc / static_cast<float>(w * h));
    const auto logAvgLum = std::log2(avgLum);
    const auto logMinLum = std::log2(minLum + 1e-6f);
    const auto logMaxLum = std::log2(maxLum + 1e-6f);

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

Image FilmicToneMapper::operator()(const Image& img) const noexcept {
    // https://knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve/
    const auto a = 2.51f;
    const auto b = 0.03f;
    const auto c = 2.43f;
    const auto d = 0.59f;
    const auto e = 0.14f;
    auto ret = Image(img.width(), img.height());
    for (std::size_t i = 0; i < img.width(); ++i) {
        for (std::size_t j = 0; j < img.height(); ++j) {
            const auto p = img(i, j);
            ret(i, j) = Color{
                std::clamp((p.r * (a * p.r + b)) / (p.r * (c * p.r + d) + e), 0.0f, 1.0f),
                std::clamp((p.g * (a * p.g + b)) / (p.g * (c * p.g + d) + e), 0.0f, 1.0f),
                std::clamp((p.b * (a * p.b + b)) / (p.b * (c * p.b + d) + e), 0.0f, 1.0f)
            };
        }
    }
    return ret;
}
