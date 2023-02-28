/*****************************************************************************
 *
 * AM utilities
 *
 * released under MIT license
 *
 * 2008-2018 André Müller
 *
 *****************************************************************************/

#ifndef AM_GENERIC_TUPLE_APPLY_H_
#define AM_GENERIC_TUPLE_APPLY_H_

#include <type_traits>
#include <utility>
#include <tuple>


namespace am {

namespace detail {

/*****************************************************************************
 *
 * @brief helper function exploding the tuple arguments into a function call
 *
 *****************************************************************************/
template<class F, class Tuple, std::size_t...ns>
inline auto
apply_helper(F&& f, Tuple&& t, std::integer_sequence<std::size_t,ns...>)
    -> decltype(f(
        std::forward<
            typename std::tuple_element<ns,typename std::decay<Tuple>::type>::type
        >(std::get<ns>(t)) ... ) )
{
    return f(
            std::forward<
            typename std::tuple_element<ns,typename std::decay<Tuple>::type>::type
        >(std::get<ns>(t)) ... );
}

}  // namespace detail



/*****************************************************************************
 *
 * @brief invokes a callable object with a tuple of arguments,
 *        similar to C++17's std::apply
 *
 *****************************************************************************/
template<class F>
inline auto
apply(F&& f, std::tuple<>) -> decltype(f())
{
    return f();
}

//---------------------------------------------------------
template<class F, class...Ts, class = typename
         std::enable_if<(sizeof...(Ts) > 0),F>::type>
inline auto
apply(F&& f, std::tuple<Ts...>& t)
    -> decltype(detail::apply_helper(std::forward<F>(f),t,std::make_index_sequence<sizeof...(Ts)>{}))
{
    return detail::apply_helper(std::forward<F>(f),t,std::make_index_sequence<sizeof...(Ts)>{});
}

//-----------------------------------------------------
template<class F, class...Ts, class = typename
         std::enable_if<(sizeof...(Ts) > 0),F>::type>
inline auto
apply(F&& f, const std::tuple<Ts...>& t)
    -> decltype(detail::apply_helper(std::forward<F>(f),t,std::make_index_sequence<sizeof...(Ts)>{}))
{
    return detail::apply_helper(std::forward<F>(f),t,std::make_index_sequence<sizeof...(Ts)>{});
}

}  // namespace am


#endif
