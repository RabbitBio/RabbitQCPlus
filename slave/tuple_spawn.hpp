#pragma once
#include <tuple>
#ifdef __sw_slave__
#include "dma_funcs.hpp"
namespace cxx11_apply
{
    template <int... Is>
    struct integer_sequence
    {
    };
    template <int N, int... Is>
    struct GenSeq : GenSeq<N - 1, N - 1, Is...>
    {
    };
    template <int... Is>
    struct GenSeq<0, Is...>
    {
        typedef integer_sequence<Is...> type;
    };
    template <typename Tt, typename Tf, int... Is>
    void __apply(Tt t, Tf f, integer_sequence<Is...> seq)
    {
        f(std::get<Is>(t)...);
    }
    template <typename... Ts>
    void apply(void (*F)(Ts...), std::tuple<Ts...> tuple)
    {
        typename GenSeq<sizeof...(Ts)>::type a;
        __apply(tuple, F, a);
    }
}

template <typename ...Ts>
__always_inline void call_tuple(void (*f)(Ts...), Ts ...arg){
  f(arg...);
}
template <typename F, typename ...Ts>
__always_inline void __tuple_spawn_proxy(std::tuple<F, Ts...> *arg) {
  cxx11_apply::apply(call_tuple<Ts...>, *arg);
}
template <typename F, typename ...Ts>
__attribute__((noinline)) void tuple_spawn_proxy(std::tuple<F, Ts...> *arg) {
  std::tuple<F, Ts...> ldm_arg;
  dma_getn(arg, &ldm_arg, 1);
  cxx11_apply::apply(call_tuple<Ts...>, ldm_arg);
}
template <typename ...Ts>
void tmpl_entrance(void (*f)(Ts...)){
  __asm__ __volatile__ (""
    :
    :"r"(tuple_spawn_proxy<void (*)(Ts...), Ts...>)
    : "memory");
}
#define ATHREAD_VISIBLE(x) template void tmpl_entrance(decltype(x)*);

#endif
#ifdef __sw_host__
extern "C"{
#include <athread.h>
}
#include <tuple>

template <typename F, typename ...Ts>
extern void slave_tuple_spawn_proxy(std::tuple<F, Ts...>*);
template <typename ...Ts>
void athread_spawn_tupled(void (*f)(Ts...), Ts ...args) {
  auto arg = std::tuple<decltype(f), Ts...>(f, args...);
  __real_athread_spawn((void*)slave_tuple_spawn_proxy<decltype(f), Ts...>, &arg);
}

#endif