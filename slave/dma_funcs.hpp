#pragma once
#include <sys/cdefs.h>
#ifndef __lsp_clang__
#define __sw_ocn__
#ifdef __sw_ocn__
typedef int desc_t __attribute__ ((__mode__(__V1XI__)));
#elif defined(__sw_thl__)
typedef int desc_t __attribute__ ((__mode__(__V1OI__)));
#endif
#else
typedef int desc_t __attribute__((vector_size(32)));
#endif
extern "C"{
// extern __attribute__((section(".ldm"))) int _MYID;
}
__thread_local desc_t dma_desc;
template <typename TM, typename TL>
__always_inline __attribute__((__artificial__)) void dma_getn(TM *mem, TL *ldm, int count) {
  //long dma_desc[4];
  // if (_MYID == 0)
  // printf("%p %p\n", mem, ldm);
  #ifdef __sw_slave__
//  desc_t dma_desc;
  asm(
      "ble %[SIZE], 2f\n\t"
      "memb\n\t"
      "stl $31, 0(%[RPL])\n\t"
      "vinsw %[RPL], $31, 2, %[DESC]\n\t"
      "ldi %[DESC], 1($31)\n\t"
      "sll %[DESC], 56, %[DESC]\n\t"
      "vinsw %[SIZE], %[DESC], 0, %[DESC]\n\t"
      "dma %[DESC], %[MEM], %[LDM]\n\t"
      "1:\n\t"
      "ldw %[DESC], 0(%[RPL])\n\t"
      "beq %[DESC], 1b\n\t"
      "2:\n\t"
      "memb\n\t"
      : [ DESC ] "=&r"(dma_desc), "+m"(*(TL (*)[count])ldm)
      : [ MEM ] "r"(mem), [ LDM ] "r"(ldm), [ RPL ] "r"(((long*)&dma_desc) + 2), [ SIZE ] "r"(count * sizeof(TL)),
        "m"(*(TL (*)[count])mem));
  #endif
}
template <typename TM, typename TL>
__always_inline __attribute__((__artificial__)) void dma_putn(TM *mem, TL *ldm, int count) {
//  desc_t dma_desc;
  #ifdef __sw_slave__
  asm(
      "ble %[SIZE], 2f\n\t"
      "memb\n\t"
      "stl $31, 0(%[RPL])\n\t"
      "vinsw %[RPL], $31, 2, %[DESC]\n\t"
      "vinsw %[SIZE], %[DESC], 0, %[DESC]\n\t"
      "dma %[DESC], %[MEM], %[LDM]\n\t"
      "1:\n\t"
      "ldw %[DESC], 0(%[RPL])\n\t"
      "beq %[DESC], 1b\n\t"
      "2:\n\t"
      "memb\n\t"
      : [ DESC ] "=&r"(dma_desc), "+m"(*(TL (*)[count])mem)
      : [ MEM ] "r"(mem), [ LDM ] "r"(ldm), [ RPL ] "r"(((long*)&dma_desc) + 2), [ SIZE ] "r"(count * sizeof(TL)),
        "m"(*(TL (*)[count])ldm));
  #endif
}

template <typename TM, typename TL>
__always_inline __attribute__((__artificial__)) void dma_get(TM *mem, TL *ldm, int size) {
  //long dma_desc[4];
  // if (_MYID == 0)
  // printf("%p %p\n", mem, ldm);
  #ifdef __sw_slave__
//  desc_t dma_desc;
  asm volatile(
      "ble %[SIZE], 2f\n\t"
      "memb\n\t"
      "stl $31, 0(%[RPL])\n\t"
      "vinsw %[RPL], $31, 2, %[DESC]\n\t"
      "ldi %[DESC], 1($31)\n\t"
      "sll %[DESC], 56, %[DESC]\n\t"
      "vinsw %[SIZE], %[DESC], 0, %[DESC]\n\t"
      "dma %[DESC], %[MEM], %[LDM]\n\t"
      "1:\n\t"
      "ldw %[DESC], 0(%[RPL])\n\t"
      "beq %[DESC], 1b\n\t"
      "2:\n\t"
      "memb\n\t"
      : [ DESC ] "=&r"(dma_desc)
      : [ MEM ] "r"(mem), [ LDM ] "r"(ldm), [ RPL ] "r"(((long*)&dma_desc) + 2), [ SIZE ] "r"(size)
      : "memory");
  #endif
}
template <typename TM, typename TL>
__always_inline __attribute__((__artificial__)) void dma_put(TM *mem, TL *ldm, int size) {
//  desc_t dma_desc;
  #ifdef __sw_slave__
  asm volatile(
      "ble %[SIZE], 2f\n\t"
      "memb\n\t"
      "stl $31, 0(%[RPL])\n\t"
      "vinsw %[RPL], $31, 2, %[DESC]\n\t"
      "vinsw %[SIZE], %[DESC], 0, %[DESC]\n\t"
      "dma %[DESC], %[MEM], %[LDM]\n\t"
      "1:\n\t"
      "ldw %[DESC], 0(%[RPL])\n\t"
      "beq %[DESC], 1b\n\t"
      "2:\n\t"
      "memb\n\t"
      : [ DESC ] "=&r"(dma_desc)
      : [ MEM ] "r"(mem), [ LDM ] "r"(ldm), [ RPL ] "r"(((long*)&dma_desc) + 2), [ SIZE ] "r"(size)
      : "memory");
  #endif
}
