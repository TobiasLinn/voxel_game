#ifndef TREE_HPP
#define TREE_HPP

#include <vector>

#include "chunk.hpp"

namespace vxg {

template<i32 ... LS>
class tree;
template<i32 LS1, i32 ... LS>
class tree<LS1, LS...> {
public:
    i32 root;
    std::tuple<std::vector<internal_chunk<LS>>...> ichunks;
    std::vector<leaf_chunk<LS1>> lchunks;

    inline tree()
        : root(0) {
    }

    inline i32 get_root() const {
        return root;
    }

    template<i32 LVL>
    inline const auto& get_chunk(i32 chptr, std::integral_constant<i32, LVL>) const {
        return std::get<LVL-1>(ichunks)[chptr];
    }
    inline const auto& get_chunk(i32 chptr, std::integral_constant<i32, 0>) const {
        return lchunks[chptr];
    }
};

}

#endif
