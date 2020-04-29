#include <RandomNumberGenerator.hpp>

std::default_random_engine& randomEngine() noexcept {
    static thread_local auto re = std::default_random_engine{std::random_device{}()};
    return re;
}
