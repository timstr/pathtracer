#include <Image.hpp>

#include <cassert>
#include <climits>
#include <cstdint>
#include <fstream>

Image::Image(std::size_t width, std::size_t height)
    : m_width(width)
    , m_height(height)
    , m_data(width * height, Color{}) {

}

std::size_t Image::width() const noexcept {
    return m_width;
}

std::size_t Image::height() const noexcept {
    return m_height;
}

namespace {
    static_assert(CHAR_BIT == 8);
    static_assert(sizeof(std::uint64_t) == 8);
    static_assert(sizeof(double) == 8);

    void writeUInt64(std::ofstream& f, std::uint64_t x) noexcept {
        f.write(reinterpret_cast<const char*>(&x), sizeof(x));
    }
    void writeFloat32(std::ofstream& f, float x) noexcept {
        f.write(reinterpret_cast<const char*>(&x), sizeof(x));
    }

    std::uint64_t readUInt64(std::ifstream& f) noexcept {
        auto x = std::uint64_t{};
        f.read(reinterpret_cast<char*>(&x), sizeof(x));
        return x;
    }
    float readFloat32(std::ifstream& f) noexcept {
        auto x = float{};
        f.read(reinterpret_cast<char*>(&x), sizeof(x));
        return x;
    }

} // anonymous namespace

void Image::save(const std::string& path) const {
    auto f = std::ofstream{path, std::ios::out | std::ios::binary};
    writeUInt64(f, width());
    writeUInt64(f, height());
    for (std::size_t y = 0, yEnd = height(); y != yEnd; ++y) {
        for (std::size_t x = 0, xEnd = width(); x != xEnd; ++x) {
            auto color = (*this)(x, y);
            writeFloat32(f, color.r);
            writeFloat32(f, color.g);
            writeFloat32(f, color.b);
        }
    }
}

Image Image::load(const std::string& path) {
    auto f = std::ifstream{path, std::ios::in| std::ios::binary};
    if (!f) {
        throw std::runtime_error("Failed to open file: \"" + path + "\"");
    }
    auto w = readUInt64(f);
    auto h = readUInt64(f);
    auto img = Image(w, h);
    for (std::size_t y = 0; y != h; ++y) {
        for (std::size_t x = 0; x != w; ++x) {
            auto color = Color{};
            color.r = readFloat32(f);
            color.g = readFloat32(f);
            color.b = readFloat32(f);
            img(x, y) = color;
        }
    }

    return img;
}

const Color& Image::operator()(std::size_t x, std::size_t y) const noexcept {
    assert(x < m_width);
    assert(y < m_height);
    return m_data[y * m_width + x];
}

Color& Image::operator()(std::size_t x, std::size_t y) noexcept {
    return const_cast<Color&>(const_cast<const Image*>(this)->operator()(x, y));
}

void Image::fill(const Color& c) noexcept {
    for (auto& d : m_data) {
        d = c;
    }
}

Image& Image::operator+=(const Image& other) noexcept {
    assert(m_width == other.m_width);
    assert(m_height == other.m_height);
    for (size_t i = 0; i < m_data.size(); ++i) {
        m_data[i] += other.m_data[i];
    }
    return *this;
}
