#pragma once

#include <random>

// TODO: use something cheaper than a Mersenne twister.
// A small LFSR should be totally fine
std::default_random_engine& randomEngine() noexcept;
