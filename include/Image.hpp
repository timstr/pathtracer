#pragma once

#include <Color.hpp>

#include <cstddef>
#include <string>
#include <vector>

class Image {
public:
    Image(std::size_t width, std::size_t height);

    std::size_t width() const noexcept;
    std::size_t height() const noexcept;

    static Image load(const std::string& path);

    void saveUncompressed(const std::string& path) const;
    static Image loadUncompressed(const std::string& path);

    const Color& operator()(std::size_t x, std::size_t y) const noexcept;
    Color& operator()(std::size_t x, std::size_t y) noexcept;

    void fill(const Color&) noexcept;

    Image& operator+=(const Image& other) noexcept;

private:
    size_t m_width;
    size_t m_height;
    std::vector<Color> m_data;
};
