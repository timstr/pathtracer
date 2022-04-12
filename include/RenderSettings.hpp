#pragma once

#include <cstddef>

class RenderSettings {
public:
    RenderSettings(std::size_t width = 256, std::size_t height = 256) noexcept;

    std::size_t width() const noexcept;
    std::size_t height() const noexcept;
    std::size_t numBounces() const noexcept;
    std::size_t samplesPerPixel() const noexcept;

    void setSize(size_t width, size_t height) noexcept;
    void setNumBounces(std::size_t) noexcept;
    void setSamplesPerPixel(std::size_t) noexcept;

private:
    std::size_t m_width;
    std::size_t m_height;
    std::size_t m_numBounces;
    std::size_t m_samplesPerPixel;
};
