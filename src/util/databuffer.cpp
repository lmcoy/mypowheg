#include "util/databuffer.h"

#include <cstring>
#include <cassert>

#include "util/byte.h"

using namespace Util;

int DataBuffer::AddUInt8(uint8_t i) {
    if(size == 0) {
        pos += 1;
        return 0;
    }
    union {
        uint8_t i;
        char c;
    } rep;
    if (pos >= size - 1) {
        return -1;
    }
    rep.i = i;
    data[pos] = rep.c;
    pos += 1;
    return 0;
}

uint8_t DataBuffer::GetUInt8() {
    assert(size != 0 && "Is DataBuffer in GETSIZE mode?");
    union {
        uint8_t i;
        char c;
    } rep;
    if (rpos >= size - 1) {
        return 0xff;
    }
    rep.c = data[rpos];
    rpos += 1;
    return rep.i;
}

int DataBuffer::AddUInt32(uint32_t i) {
    if(size == 0) {
        pos += 4;
        return 0;
    }
    union {
        uint32_t i;
        char c[4];
    } rep;
    if (pos >= size - 4) {
        return -1;
    }
    rep.i = htole32(i);
    memcpy(data+pos, rep.c, 4);
    pos += 4;
    return 0;
}

uint32_t DataBuffer::GetUInt32() {
    assert(size != 0 && "Is DataBuffer in GETSIZE mode?");
    union {
        uint32_t i;
        char c[4];
    } rep;
    if (rpos >= size - 4) {
        return 0xffffffff;
    }
    memcpy(rep.c, data+rpos, 4);
    rpos += 4;
    return le32toh(rep.i);
}

int DataBuffer::AddUInt64(uint64_t i) {
    if(size == 0) {
        pos += 8;
        return 0;
    }
    union {
        uint64_t i;
        char c[8];
    } rep;
    if (pos >= size - 8) {
        return -1;
    }
    rep.i = htole64(i);
    memcpy(data+pos, rep.c, 8);
    pos += 8;
    return 0;
}

uint64_t DataBuffer::GetUInt64() {
    assert(size != 0 && "Is DataBuffer in GETSIZE mode?");
    union {
        uint64_t i;
        char c[8];
    } rep;
    if(rpos >= size - 8) {
        return 0xffffffffffffffffull;
    }
    memcpy(rep.c, data+rpos, 8);
    rpos += 8;
    return le64toh(rep.i);
}

int DataBuffer::AddInt32(int32_t i) {
    union {
        uint32_t u;
        int32_t i;
    } rep;
    rep.i = i;
  uint32_t le = htole32(rep.u);
  return AddUInt32(le);
}

int DataBuffer::AddInt64(int64_t i) {
    union {
        uint64_t u;
        int64_t i;
    } rep;
    rep.i = i;
  uint64_t le = htole64(rep.u);
  return AddUInt64(le);
}

int32_t DataBuffer::GetInt32() {
    assert(size != 0 && "Is DataBuffer in GETSIZE mode?");
    union {
        uint32_t u;
        int32_t i;
    } rep;
    rep.u = GetUInt32();
    return rep.i;
}

int64_t DataBuffer::GetInt64() {
    assert(size != 0 && "Is DataBuffer in GETSIZE mode?");
    union {
        uint64_t u;
        int64_t i;
    } rep;
    rep.u = GetUInt64();
    return rep.i;
}

int DataBuffer::AddDouble(double d) {
   uint64_t u = pack_double(d);
   return AddUInt64(u); 
}

double DataBuffer::GetDouble() {
    assert(size != 0 && "Is DataBuffer in GETSIZE mode?");
    uint64_t u = GetUInt64();
    return unpack_double(u);
}

int DataBuffer::AddDoubleArray(int n, double * array) {
    if (size && pos >= size - (int)sizeof(double) * n) {
        return -1;
    }
    for(int i = 0; i < n; i++) {
        int in = AddDouble(array[i]);
        assert(in == 0);
    }
    return 0;
}

void DataBuffer::GetDoubleArray(int n, double * array) {
    assert(size != 0 && "Is DataBuffer in GETSIZE mode?");
    if (rpos >= size - (int)sizeof(double) * n) {
        return;
    }
    for(int i = 0; i < n; i++) {
        array[i] = GetDouble();
    }
}

int DataBuffer::AddIntArray(int n, int *array) {
    if (size && rpos >= size - (int)sizeof(uint64_t) * n) {
        return -1;
    }
    for (int i = 0; i < n; i++) {
        int in = AddInt64((uint64_t)array[i]);
        assert(in == 0);
    }
    return 0;
}

void DataBuffer::GetIntArray(int n, int *array) {
    assert(size != 0 && "Is DataBuffer in GETSIZE mode?");
    if (rpos >= size - (int)sizeof(double) * n) {
        return;
    }
    for(int i = 0; i < n; i++) {
        array[i] = (int)GetInt64();
    }
}
