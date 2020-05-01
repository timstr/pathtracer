#pragma once

#include <Geometry.hpp>

#include <memory>
#include <optional>
#include <utility>
#include <vector>

class Object;

namespace ObjectTree {

    class Node {
    public:
        Node(Box) noexcept;
        virtual ~Node() noexcept = default;

        virtual std::optional<std::pair<float, const Object*>> hit(const Ray&) const noexcept = 0;

        Box const boundingBox;
    };

    class InternalNode : public Node {
    public:
        InternalNode(std::unique_ptr<Node> childA, std::unique_ptr<Node> childB) noexcept;

    private:
        std::optional<std::pair<float, const Object*>> hit(const Ray&) const noexcept override;

        std::unique_ptr<Node> const m_childA;
        std::unique_ptr<Node> const m_childB;
    };

    class LeafNode : public Node {
    public:
        LeafNode(const Object* object) noexcept;

    private:
        std::optional<std::pair<float, const Object*>> hit(const Ray&) const noexcept override;

        const Object* const m_object;
    };

    class Tree {
    public:
        Tree(std::vector<const Object*> objects);

        std::optional<std::pair<float, const Object*>> hit(const Ray&) const noexcept;

    private:
        static std::unique_ptr<Node> makeTree(std::vector<const Object*> objects, float volumeVersusSplitCost);

        std::unique_ptr<Node> m_root;
    };

} // namespace ObjectTree