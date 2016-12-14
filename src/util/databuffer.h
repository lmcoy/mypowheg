#ifndef UTIL_DATABUFFER_H
#define UTIL_DATABUFFER_H

#include <cstdint>
#include <cassert>

namespace Util { 

class DataBuffer {
    public:
        static const int GETSIZE = 0;
        explicit DataBuffer(int s) {
            assert(s >= 0);
            pos = 0;
            rpos = 0;
            if (s != GETSIZE) {
                data = new char[s];
            }
            size = s;
        }
        ~DataBuffer() {
            if(size != 0) {
                delete[] data;
            }
        }
        int AddUInt8(uint8_t i);
        int AddInt32(int32_t i);
        int AddUInt32(uint32_t i);
        int AddUInt64(uint64_t i);
        int AddInt64(int64_t i);
        int AddInt32Array(int size, int32_t *array);
        int AddDouble(double d);
        int AddDoubleArray(int size, double *array);
        int AddIntArray(int size, int * array);

        uint8_t GetUInt8();
        uint32_t GetUInt32();
        uint64_t GetUInt64();
        int32_t GetInt32();
        int64_t GetInt64();

        double GetDouble();
        void GetDoubleArray(int size, double *array);
        void GetIntArray(int size, int *array);

        char *Data() {
            if (size != 0) {
                return data;
            }
            return 0;
        }

        int GetDataSize() const { return pos; }

    private:
        char * data;
        int pos;
        int rpos;
        int size;
};

} // namespace Util

#endif

