#ifndef BYTE_H_
#define BYTE_H_

#ifdef __cplusplus
#include <cassert>
#include <cstdint>
#else
#include <assert.h>
#include <stdint.h>
#endif

static int is_big_endian(void) {
  union {
    uint32_t i;
    char c[4];
  } char_rep = { .i = 1 };
  return char_rep.c[3];
}

static uint32_t swap_u32(uint32_t be) {
  return ((be >> 24) & 0xff) | ((be >> 8) & 0xff00) | ((be << 8) & 0xff0000) |
         ((be << 24) & 0xff000000);
}

static uint64_t swap_u64(uint64_t u) {
  uint64_t r = 0;
  r |= (u >> 56) & 0x00000000000000ffull;
  r |= (u >> 40) & 0x000000000000ff00ull;
  r |= (u >> 24) & 0x0000000000ff0000ull;
  r |= (u >> 8) & 0x00000000ff000000ull;
  r |= (u << 8) & 0x000000ff00000000ull;
  r |= (u << 24) & 0x0000ff0000000000ull;
  r |= (u << 40) & 0x00ff000000000000ull;
  r |= (u << 56) & 0xff00000000000000ull;
  return r;
}

uint32_t htole32(uint32_t host) {
  if (!is_big_endian()) {
    return host;
  }
  return swap_u32(host);
}

uint32_t le32toh(uint32_t le) {
  if (is_big_endian()) {
    return swap_u32(le);
  }
  return le;
}

uint64_t htole64(uint64_t host) {
  if (!is_big_endian()) {
    return host;
  }
  return swap_u64(host);
}

uint64_t le64toh(uint64_t le) {
  if (is_big_endian()) {
    return swap_u64(le);
  }
  return le;
}

// BUG: assumes that the float representation is the same
// on systems which call pack and unpack.
static uint64_t pack_double(double d) {
    assert(sizeof(double)==sizeof(uint64_t));
    union {
        uint64_t u;
        double d;
    } rep;
    rep.d = d;
    return htole64(rep.u);
}

static double unpack_double(uint64_t u) {
    assert(sizeof(double)==sizeof(uint64_t));
    union { 
        uint64_t u;
        double d;
    } rep;
    rep.u = le64toh(u);
    return rep.d;
}

#endif

