#ifndef VAR_TO_TEMPLATE_HPP
#define VAR_TO_TEMPLATE_HPP

#include <type_traits>

#include "util.hpp"

namespace vxg {

// convert variables to template parameters
template<template<class> class CB, int N, class... TPARAMS>
struct var_to_template {
    template<class R, class... ARGS>
    static inline R impl(bool b, ARGS&&... args) {
        if (b) {
            return var_to_template<CB, N-1, TPARAMS..., std::true_type>::template impl<R>(std::forward<ARGS>(args)...);
        } else {
            return var_to_template<CB, N-1, TPARAMS..., std::false_type>::template impl<R>(std::forward<ARGS>(args)...);
        }
    }
    
    template<class R, class... ARGS>
    static inline R impl(i32 i, ARGS&&... args) {
        switch (i) {
        case -1:
            return var_to_template<CB, N-1, TPARAMS..., std::integral_constant<i32,-1>>::template impl<R>(std::forward<ARGS>(args)...);
        case  0:
            return var_to_template<CB, N-1, TPARAMS..., std::integral_constant<i32, 0>>::template impl<R>(std::forward<ARGS>(args)...);
        case  1:
            return var_to_template<CB, N-1, TPARAMS..., std::integral_constant<i32, 1>>::template impl<R>(std::forward<ARGS>(args)...);
        case  2:
            return var_to_template<CB, N-1, TPARAMS..., std::integral_constant<i32, 2>>::template impl<R>(std::forward<ARGS>(args)...);
        }
        throw "var_to_template: Unsupported i32 value: " + std::to_string(i);
    }
    
    template<class R, class... ARGS>
    static inline R impl(u32 i, ARGS&&... args) {
        switch (i) {
        case  0:
            return var_to_template<CB, N-1, TPARAMS..., std::integral_constant<u32,0>>::template impl<R>(std::forward<ARGS>(args)...);
        case  1:
            return var_to_template<CB, N-1, TPARAMS..., std::integral_constant<u32,1>>::template impl<R>(std::forward<ARGS>(args)...);
        case  2:
            return var_to_template<CB, N-1, TPARAMS..., std::integral_constant<u32,2>>::template impl<R>(std::forward<ARGS>(args)...);
        }
        throw "var_to_template: Unsupported u32 value: " + std::to_string(i);
    }
};
template<template<class> class CB, class... TPARAMS>
struct var_to_template<CB, 0, TPARAMS...> {
    template<class R, class... ARGS>
    static inline R impl(ARGS&&... args) {
        return CB<TPARAMS...>::template impl<R>(std::forward<ARGS>(args)...);
    }
};

}

#endif
