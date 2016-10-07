#ifndef AUDI_BACK_COMPATIBILITY_HPP
#define AUDI_BACK_COMPATIBILITY_HPP

#include<type_traits>

namespace audi{
inline namespace impl
{

template <bool B, typename T = void>
using enable_if_t = typename std::enable_if<B,T>::type;

template< class T >
using decay_t = typename std::decay<T>::type;

}
}
#endif
