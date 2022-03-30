#include <ObjectTree.hpp>

#include <PathTracer.hpp>

#include <cassert>

namespace ObjectTree {


    Node::Node(Box bb) noexcept
        : boundingBox(bb) {

    }

    InternalNode::InternalNode(std::unique_ptr<Node> childA, std::unique_ptr<Node> childB) noexcept
        : Node(boxContaining(childA->boundingBox, childB->boundingBox))
        , m_childA(std::move(childA))
        , m_childB(std::move(childB)) {

    }

    std::optional<std::pair<float, const Object*>> InternalNode::hit(const Ray& ray) const noexcept {
        const auto ta = intersect(ray, m_childA->boundingBox);
        const auto tb = intersect(ray, m_childB->boundingBox);
        if (ta.has_value()) {
            if (tb.has_value()) {
                const auto aIsCloser = *ta < *tb;
                const auto& near = aIsCloser ? m_childA : m_childB;
                const auto& far = aIsCloser ? m_childB : m_childA;

                const auto c = near->hit(ray);
                if (c.has_value()) {
                    const auto p = ray.pos + c->first * ray.dir;
                    if (inside(p, far->boundingBox)) {
                        const auto c2 = far->hit(ray);
                        if (c2.has_value() && c2->first < c->first) {
                            return c2;
                        }
                    }

                    return c;
                }
                return far->hit(ray);
            } else {
                return m_childA->hit(ray);
            }
        } else {
            if (tb.has_value()) {
                return m_childB->hit(ray);
            } else {
                return std::nullopt;
            }
        }
        
    }

    LeafNode::LeafNode(const Object* object) noexcept
        : Node(object->getBoundingBox())
        , m_object(object) {

    }

    std::optional<std::pair<float, const Object*>> LeafNode::hit(const Ray& ray) const noexcept {
        const auto t = m_object->hit(ray);
        if (t.has_value()) {
            return std::pair{ *t, m_object };
        }
        return std::nullopt;
    }


    Tree::Tree(std::vector<const Object*> objects)
        : m_root(makeTree(std::move(objects), 0.5f)) {

    }

    std::optional<std::pair<float, const Object*>> Tree::hit(const Ray& ray) const noexcept {
        assert(m_root);
        return m_root->hit(ray);
    }

    std::unique_ptr<Node> Tree::makeTree(std::vector<const Object*> objects, float volumeVersusSplitCost) {
        // TODO: consider caching bounding boxes i.e. in a std::unordered_map<const Object*, BoundingBox>

        if (objects.size() == 0) {
            return nullptr;
        }

        if (objects.size() == 1) {
            return std::make_unique<LeafNode>(objects[0]);
        }

        if (objects.size() == 2) {
            return std::make_unique<InternalNode>(
                std::make_unique<LeafNode>(objects[0]),
                std::make_unique<LeafNode>(objects[1])
            );
        }

        using Dimension = const float Pos::*;
        auto bestSplitDim = Dimension{nullptr};
        auto bestSplitIndex = std::optional<std::size_t>{};
        auto bestCost = std::numeric_limits<float>::max();

        const auto computeBestSplit = [&](Dimension dim) {
            for (std::size_t i = 1; i < objects.size(); ++i) {
                const auto bb = objects[i]->getBoundingBox();
                const auto m = bb.center.*dim;

                // Bounding boxes containing partitions about current object
                auto bba = std::optional<Box>{};
                auto bbb = std::optional<Box>{};

                for (std::size_t j = 0; j < objects.size(); ++j) {
                    const auto bb2 = objects[j]->getBoundingBox();
                    const auto m2 = bb2.center.*dim;

                    if (m2 < m) {
                        bba = bba.has_value() ? boxContaining(*bba, bb2) : bb2;
                    } else {
                        bbb = bbb.has_value() ? boxContaining(*bbb, bb2) : bb2;
                    }
                }

                if (!(bba.has_value() && bbb.has_value())) {
                    continue;
                }

                const auto vInner = bba->volume() + bbb->volume();
                const auto vTotal = boxContaining(*bba, *bbb).volume();
                const auto volumeCost = vInner / vTotal;

                const auto splitRatio = static_cast<float>(i) / static_cast<float>(objects.size() - i);
                const auto splitCost = std::max(splitRatio, 1.0f / splitRatio);

                const auto cost = volumeVersusSplitCost * volumeCost + (1.0f - volumeVersusSplitCost) * splitCost;

                if (cost < bestCost) {
                    bestCost = cost;
                    bestSplitDim = dim;
                    bestSplitIndex = i;
                }
            }
        };

        computeBestSplit(&Pos::x);
        computeBestSplit(&Pos::z);
        computeBestSplit(&Pos::y);

        assert(bestSplitDim != nullptr);
        assert(bestSplitIndex.has_value());

        auto splitA = std::vector<const Object*>{};
        auto splitB = std::vector<const Object*>{};

        const auto bb = objects[*bestSplitIndex]->getBoundingBox();
        const auto m = bb.center.*bestSplitDim;

        for (const auto& o : objects) {
            const auto bb2 = o->getBoundingBox();
            const auto m2 = bb2.center.*bestSplitDim;
            (m2 < m ? splitA : splitB).push_back(o);
        }

        return std::make_unique<InternalNode>(
            makeTree(splitA, volumeVersusSplitCost),
            makeTree(splitB, volumeVersusSplitCost)
        );
    }


} // namespace ObjectTree
