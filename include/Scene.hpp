#pragma once

#include <Color.hpp>
#include <Object.hpp>
#include <ObjectTree.hpp>
#include <LinearAlgebra.hpp>

#include <cstddef>
#include <memory>
#include <vector>

class Scene {
public:
    Scene();

    void addObject(std::unique_ptr<Object> object) noexcept;

    template<typename T, typename... Args>
    T& addObject(Args&&... args) {
        auto up = std::make_unique<T>(std::forward<Args>(args)...);
        auto& r = *up;
        addObject(std::move(up));
        return r;
    }

    Color trace(Ray ray, std::size_t depth) const noexcept;

    void updateGeometry();

private:
    std::vector<std::unique_ptr<Object>> m_objects;

    ObjectTree::Tree m_objectTree;
};
