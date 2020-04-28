#include <PathTracer.hpp>

std::tuple<Ray, Color, float> BasicGlossyMaterial::deflect(const Vec& inbound, const Vec& normal) const noexcept {
    // TODO
    auto diffuseBounce = randomHemisphereVectorUniform(normal);
    auto specularBounce = a a aa a a a a a a
}


std::optional<float> TriangleObject::hit(const Ray& ray) const noexcept {
    return intersect(ray, geometry);
}

std::tuple<Ray, Color, float> TriangleObject::deflect(const Ray& ray) const noexcept {
    // TODO
    a a a a a a 
}

std::optional<float> SphereObject::hit(const Ray& ray) const noexcept {
    return intersect(ray, geometry);
}

std::tuple<Ray, Color, float> SphereObject::deflect(const Ray& ray) const noexcept {
    // TODO
    a a a a a
}

Color Scene::trace(Ray ray, std::size_t depth) const noexcept {
    auto color = Color{};
    auto attenuation = 1.0f;
    for (std::size_t i = 0; i < depth; ++i) {
        auto closestT = -1.0f;
        const Object* hitObj = nullptr;
        for (const auto& obj : m_objects) {
            if (auto t = obj->hit(ray)) {
                if (*t < closestT) {
                    hitObj = obj.get();
                    closestT = *t;
                }
            }
        }
        if (!hitObj) {
            return color;
        }
        // TODO: overload +=
        ray.pos = ray.pos + closestT * ray.dir;
        const auto [outboundRay, emitted, att] = hitObj->deflect(ray);
        // TODO: overload operators for Color
        color.r += emitted.r * attenuation;
        color.g += emitted.g * attenuation;
        color.b += emitted.b * attenuation;
        attenuation *= att;
        ray = outboundRay;
    }
    return color;
}