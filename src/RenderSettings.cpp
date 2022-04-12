#include <RenderSettings.hpp>

#include <cassert>

RenderSettings::RenderSettings(std::size_t width, std::size_t height) noexcept
    : m_width(width)
    , m_height(height)
    , m_numBounces(8)
    , m_samplesPerPixel(1) {

    assert(m_width > 0);
    assert(m_height > 0);

}

std::size_t RenderSettings::width() const noexcept {
    return m_width;
}

std::size_t RenderSettings::height() const noexcept {
    return m_height;
}

std::size_t RenderSettings::numBounces() const noexcept {
    return m_numBounces;
}

std::size_t RenderSettings::samplesPerPixel() const noexcept {
    return m_samplesPerPixel;
}

void RenderSettings::setSize(size_t w, size_t h) noexcept {
    assert(w > 0);
    assert(h > 0);
    m_width = w;
    m_height = h;
}

void RenderSettings::setNumBounces(std::size_t n) noexcept {
    assert(n > 0);
    m_numBounces = n;
}

void RenderSettings::setSamplesPerPixel(std::size_t s) noexcept {
    assert(s > 0);
    m_samplesPerPixel = s;
}
