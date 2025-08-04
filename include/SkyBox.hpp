#pragma once

#include <Image.hpp>
#include <Object.hpp>
#include <cmath>

class SkyBox : public Object {
public:
    SkyBox(
        const std::string& path_x_positive,
        const std::string& path_x_negative,
        const std::string& path_y_positive,
        const std::string& path_y_negative,
        const std::string& path_z_positive,
        const std::string& path_z_negative
    );

    std::optional<Pos> hitLocalRay(const Ray& ray) const noexcept override;

    ColorBounce deflectLocalRay(const Ray& ray) const noexcept override;

    AxisAlignedBox getLocalBoundingBox() const noexcept override;

private:
    Image m_x_positive;
    Image m_x_negative;
    Image m_y_positive;
    Image m_y_negative;
    Image m_z_positive;
    Image m_z_negative;
};