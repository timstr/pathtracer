#include <SkyBox.hpp>

#include <algorithm>
#include <cassert>

SkyBox::SkyBox(
    const std::string& path_x_positive,
    const std::string& path_x_negative,
    const std::string& path_y_positive,
    const std::string& path_y_negative,
    const std::string& path_z_positive,
    const std::string& path_z_negative
)
    : m_x_positive(Image::load(path_x_positive))
    , m_x_negative(Image::load(path_x_negative))
    , m_y_positive(Image::load(path_y_positive))
    , m_y_negative(Image::load(path_y_negative))
    , m_z_positive(Image::load(path_z_positive))
    , m_z_negative(Image::load(path_z_negative))
{

}

// TODO: this should be in/with Image
Color sample_uv(const Image& image, float u, float v) noexcept {
    auto x = static_cast<size_t>(std::round(
        std::clamp(u, 0.0f, 1.0f) * static_cast<float>(image.width() - 1)
        // + 0.5f
    ));
    auto y = static_cast<size_t>(std::round(
        std::clamp(v, 0.0f, 1.0f) * static_cast<float>(image.height() - 1)
        // + 0.5f
    ));
    assert(x < image.width());
    assert(y < image.height());
    auto c = image(x, y);
    auto avg = (c.r + c.g + c.b) / 3.0f;
    auto k = 4.0f * std::pow(avg, 2.5f);
    return Color(
        k * c.r,
        k * c.g,
        k * c.b
    );
}

std::optional<Pos> SkyBox::hitLocalRay(const Ray& ray) const noexcept {
    return (ray.dir * 1e10f).toPos();
}

ColorBounce SkyBox::deflectLocalRay(const Ray& ray) const noexcept {
    const auto abs_x = std::abs(ray.dir.x);
    const auto abs_y = std::abs(ray.dir.y);
    const auto abs_z = std::abs(ray.dir.z);

    auto bounce = ColorBounce {
        Color(0.0f, 0.0f, 0.0f),
        Color(0.0f, 0.0f, 0.0f),
        ray.dir,
        ray.dir
    };

    if (abs_x > abs_y && abs_x > abs_z) {
        // x is maximum
        if (ray.dir.x > 0.0f) {
            bounce.emitted = sample_uv(
                m_x_positive,
                0.5 - 0.5 * (ray.dir.z / abs_x),
                0.5 + 0.5 * (ray.dir.y / abs_x)
            );
        } else {
            bounce.emitted = sample_uv(
                m_x_negative,
                0.5 + 0.5 * (ray.dir.z / abs_x),
                0.5 + 0.5 * (ray.dir.y / abs_x)
            );
        }
    } else if (abs_y > abs_x && abs_y > abs_z) {
        // y is maximum
        if (ray.dir.y > 0.0f) {
            bounce.emitted = sample_uv(
                m_y_positive,
                0.5 - 0.5 * (ray.dir.x / abs_y),
                0.5 + 0.5 * (ray.dir.z / abs_y)
            );
        } else {
            bounce.emitted = sample_uv(
                m_y_negative,
                0.5 + 0.5 * (ray.dir.x / abs_y),
                0.5 + 0.5 * (ray.dir.z / abs_y)
            );
        }
    } else {
        // z must be maximum
        if (ray.dir.z > 0.0f) {
            bounce.emitted = sample_uv(
                m_z_positive,
                0.5 + 0.5 * (ray.dir.x / abs_z),
                0.5 + 0.5 * (ray.dir.y / abs_z)
            );
        } else {
            bounce.emitted = sample_uv(
                m_z_negative,
                0.5 - 0.5 * (ray.dir.x / abs_z),
                0.5 + 0.5 * (ray.dir.y / abs_z)
            );
        }
    }
    return bounce;
}

AxisAlignedBox SkyBox::getLocalBoundingBox() const noexcept {
    // Uhhhhhhhhhhhhhh my BVH is gonna hate this
    return AxisAlignedBox(Pos(0.0f, 0.0f, 0.0f), Vec(1e10f, 1e10f, 1e10f));
}