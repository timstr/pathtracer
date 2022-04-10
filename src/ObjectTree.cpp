#include <ObjectTree.hpp>

#include <PathTracer.hpp>

#include <cassert>

namespace ObjectTree {


    Node::Node(AxisAlignedBox bb) noexcept
        : boundingBox(bb) {

    }

    InternalNode::InternalNode(std::unique_ptr<Node> childA, std::unique_ptr<Node> childB) noexcept
        : Node(boxContaining(childA->boundingBox, childB->boundingBox))
        , m_childA(std::move(childA))
        , m_childB(std::move(childB)) {

    }

    std::optional<std::pair<Pos, const Object*>> InternalNode::hit(const Ray& ray) const noexcept {
        // auto pa = std::optional<Pos>{};
        // if (inside(ray.pos, m_childA->boundingBox)) {
        //     pa = ray.pos;
        // } else {
        //     pa = intersect(ray, m_childA->boundingBox);
        // }
        // 
        // auto pb = std::optional<Pos>{};
        // if (inside(ray.pos, m_childB->boundingBox)) {
        //     pb = ray.pos;
        // } else {
        //     pb = intersect(ray, m_childB->boundingBox);
        // }
        // 
        // if (pa.has_value()) {
        //     if (pb.has_value()) {
                const auto hitA = m_childA->hit(ray);
                const auto hitB = m_childB->hit(ray);
                if (hitA.has_value()) {
                    if (hitB.has_value()) {
                        const auto ta = (hitA->first - ray.pos) * ray.dir;
                        const auto tb = (hitB->first - ray.pos) * ray.dir;
                        if (tb < ta) {
                            return hitB;
                        }
                    }
                    // hitB == nullopt or hitB is further
                    return hitA;
                }
                // hitA == nullopt
                return hitB;
        //     }
        //     // pb == nullopt
        //     return m_childA->hit(ray);
        // }
        // // pa == nullopt
        // if (pb.has_value()) {
        //     return m_childB->hit(ray);
        // }
        // return std::nullopt;
    }

    LeafNode::LeafNode(const Object* object) noexcept
        : Node(object->getBoundingBox())
        , m_object(object) {

    }

    std::optional<std::pair<Pos, const Object*>> LeafNode::hit(const Ray& ray) const noexcept {
        const auto p = m_object->hitRay(ray);
        if (p.has_value()) {
            return std::pair{ *p, m_object };
        }
        return std::nullopt;
    }


    Tree::Tree(std::vector<const Object*> objects)
        : m_root(makeTree(std::move(objects), 0.5f)) {

    }

    std::optional<std::pair<Pos, const Object*>> Tree::hit(const Ray& ray) const noexcept {
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
                auto bba = std::optional<AxisAlignedBox>{};
                auto bbb = std::optional<AxisAlignedBox>{};

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

                const auto va = Rectangle{bba->halfSize}.volume();
                const auto vb = Rectangle{bbb->halfSize}.volume();
                assert(va > 0.0f);
                assert(vb > 0.0f);
                const auto vInner = va + vb;
                const auto vTotal = Rectangle{boxContaining(*bba, *bbb).halfSize}.volume();
                assert(vTotal > 1e-6f);
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

        assert((bestSplitDim == nullptr) != bestSplitIndex.has_value());
        // TODO: fall-back in case objects are co-located
        assert(bestSplitDim != nullptr);

        auto splitA = std::vector<const Object*>{};
        auto splitB = std::vector<const Object*>{};

        const auto bb = objects[*bestSplitIndex]->getBoundingBox();
        const auto m = bb.center.*bestSplitDim;

        for (const auto& o : objects) {
            const auto bb2 = o->getBoundingBox();
            const auto m2 = bb2.center.*bestSplitDim;
            (m2 < m ? splitA : splitB).push_back(o);
        }
        // TODO: fall-back partition in case splitting failed

        return std::make_unique<InternalNode>(
            makeTree(splitA, volumeVersusSplitCost),
            makeTree(splitB, volumeVersusSplitCost)
        );
    }


} // namespace ObjectTree
