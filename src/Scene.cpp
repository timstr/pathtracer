#include <Scene.hpp>

#include <cassert>

namespace {
    constexpr float epsilon = 1e-6f;
}

void Scene::addObject(std::unique_ptr<Object> object) noexcept {
    m_objects.push_back(std::move(object));
}

Scene::Scene()
    : m_objectTree(std::vector<const Object*>{}) {

}

Color Scene::trace(Ray ray, std::size_t depth) const noexcept {
    auto color = Color{};
    auto attenuation = Color{1.0f, 1.0f, 1.0f};
    for (std::size_t i = 0; i < depth; ++i) {
        // Brute force search
        // auto closestT = std::numeric_limits<float>::max();
        // auto closestPos = Pos{};
        // const Object* hitObj = nullptr;
        // for (const auto& obj : m_objects) {
        //     auto bb = obj->getBoundingBox();
        //     if (!intersect(ray, bb)) {
        //         continue;
        //     }
        //     if (auto p = obj->hitRay(ray)) {
        //         const auto t = (*p - ray.pos) * ray.dir;
        //         if (t < closestT) {
        //             hitObj = obj.get();
        //             closestT = t;
        //             closestPos = *p;
        //         }
        //     }
        // }
        // if (!hitObj) {
        //     return color;
        // }

        // Tree search
        const auto hit = m_objectTree.hit(ray);
        if (!hit) {
            return color;
        }
        const auto& [closestPos, hitObj] = *hit;

        // TODO: overload +=
        ray.pos = closestPos;
        const auto bounce = hitObj->deflectRay(ray);
        ray.pos += 1e-3f * bounce.rayDirection;
        // TODO: find a way to reliably prevent accidental transmittance
        assert(bounce.emitted.r >= 0.0f);
        assert(bounce.emitted.g >= 0.0f);
        assert(bounce.emitted.b >= 0.0f);
        assert(bounce.attenuation.r >= 0.0f && bounce.attenuation.r <= 1.0f);
        assert(bounce.attenuation.g >= 0.0f && bounce.attenuation.g <= 1.0f);
        assert(bounce.attenuation.b >= 0.0f && bounce.attenuation.b <= 1.0f);
        assert(std::abs(1.0f - bounce.rayDirection.norm()) < epsilon);

        // TODO: overload operators for Color
        color.r += bounce.emitted.r * attenuation.r;
        color.g += bounce.emitted.g * attenuation.g;
        color.b += bounce.emitted.b * attenuation.b;
        attenuation.r *= bounce.attenuation.r;
        attenuation.g *= bounce.attenuation.g;
        attenuation.b *= bounce.attenuation.b;
        ray.dir = bounce.rayDirection;

        if (attenuation.r + attenuation.b + attenuation.g < epsilon) {
            return color;
        }
    }
    return color;
}

void Scene::updateGeometry() {
    std::vector<const Object*> v;
    v.reserve(m_objects.size());
    for (const auto& o : m_objects) {
        v.push_back(o.get());
    }
    m_objectTree = ObjectTree::Tree{v};
}
