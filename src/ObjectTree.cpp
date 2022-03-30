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
        const auto pa = intersect(ray, m_childA->boundingBox);
        const auto pb = intersect(ray, m_childB->boundingBox);
        if (pa.has_value()) {
            // Bounding box A was hit
            const auto ta = ((*pa) - ray.pos) * ray.dir;
            assert(ta >= 0.0f);
            if (pb.has_value()) {
                // Bounding box B was hit
                const auto tb = ((*pb) - ray.pos) * ray.dir;
                assert(tb >= 0.0f);

                const auto aIsCloser = ta < tb;
                const auto& near = aIsCloser ? m_childA : m_childB;
                const auto& far = aIsCloser ? m_childB : m_childA;

                const auto pn = near->hit(ray);
                if (pn.has_value()) {
                    // near object was hit
                    assert(inside(pn->first, near->boundingBox));
                    // NOTE: because both objects' bounding boxes intercept the ray, both need to be fully checked
                    // near object is also inside far bounding box, so far object might also be hit
                    const auto pf = far->hit(ray);
                    if (pf.has_value()) {
                        assert(inside(pf->first, far->boundingBox));
                        // far object was hit
                        const auto cn = (pn->first - ray.pos) * ray.dir;
                        const auto cf = (pf->first - ray.pos) * ray.dir;
                        assert(cn >= 0.0f);
                        assert(cf >= 0.0f);
                        if (cf < cn) {
                            // far object is closer than near object
                            return pf;
                        }
                    }
                    return pn;
                } else {
                    // Near object was not hit
                    return far->hit(ray);
                }
            } else {
                // Bounding box A was hit but bounding box B was not hit
                return m_childA->hit(ray);
            }
        } else {
            // Bounding box A was not hit
            if (pb.has_value()) {
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
                assert(va > 1e-6f);
                assert(vb > 1e-6f);
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
