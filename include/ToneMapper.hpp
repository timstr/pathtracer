#pragma once

#include <Image.hpp>

class ReinhardToneMapper{
public:
    Image operator()(const Image&) const noexcept;
};

class FilmicToneMapper {
public:
    Image operator()(const Image&) const noexcept;
};
