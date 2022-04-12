#pragma once

#include <Image.hpp>

class ToneMapper {
public:
    virtual ~ToneMapper() noexcept = default;

    virtual Image operator()(const Image&) const noexcept = 0;
};

class ReinhardToneMapper : public ToneMapper {
public:
    virtual Image operator()(const Image&) const noexcept override final;
};
