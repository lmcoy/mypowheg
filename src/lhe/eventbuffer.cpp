#include <cstdio>
#include <cmath>
#include <cstring>
#include <cstdlib>

#include "lhe/eventbuffer.h"

#include "lhe/event.h"
namespace LHE {

EventBuffer::AllocateStatus EventBuffer::Allocate(size_t cap_in_MB) {
    if (data_ != NULL) {
        return AllocateStatus::EXISTS;
    }
    nmax_ = 0x7fffffff;
    cap_ = cap_in_MB * 1024 * 1024;
    grow_cap_ = cap_ / 2;
    len_ = 0;
    data_ = (char *)malloc(cap_ * sizeof(char));
    if (data_ == NULL) {
        return AllocateStatus::MALLOCERROR;;
    }
    return AllocateStatus::SUCCESS;
}

EventBuffer::AppendStatus EventBuffer::Append(const LHE::Event &event) {
    static const int SIZE = 1024;
    if (n_ >= nmax_) {
        return AppendStatus::FULL;
    }
    int len = 0;
    char buffer[SIZE];
    int NUP = event.N;
    int IDPRUP = event.ID;
    double XWGTUP = event.Weight;
    double SCALUP = event.Scale;
    double AQEDUP = event.Alpha;
    double AQCDUP = event.AlphaS;

    len = snprintf(buffer, SIZE, "<event>\n %d %d %E %E %E %E\n", NUP, IDPRUP,
                   XWGTUP, SCALUP, AQEDUP, AQCDUP);
    if (len >= SIZE) {
        return AppendStatus::INTERNALERROR;
    }
    for (int i = 0; i < NUP; i++) {
        int IDUP = event.Particles[i].PDG;
        int ISTUP = event.Particles[i].Status;
        int MOTHUP1 = event.Particles[i].Mother1;
        int MOTHUP2 = event.Particles[i].Mother2;
        int ICOLUP1 = event.Particles[i].Color1;
        int ICOLUP2 = event.Particles[i].Color2;
        int l = snprintf(buffer + len, SIZE - len, " %3d %2d %2d %d %3d %3d",
                         IDUP, ISTUP, MOTHUP1, MOTHUP2, ICOLUP1, ICOLUP2);
        if (l >= SIZE - len) {
            return AppendStatus::INTERNALERROR;
        }
        len += l;
        double PX = event.Particles[i].Momentum.PX();
        double PY = event.Particles[i].Momentum.PY();
        double PZ = event.Particles[i].Momentum.PZ();
        double E = event.Particles[i].Momentum.E();
        double M2 = E * E - PX * PX - PY * PY - PZ * PZ;
        double M = copysign(sqrt(fabs(M2)), M2);
        double VTIMUP = event.Particles[i].LifeTime;
        double SPINUP = event.Particles[i].Spin;
        l = snprintf(buffer + len, SIZE - len,
                     " %18.10E %18.10E %18.10E %18.10E %18.10E %G %G\n", PX, PY,
                     PZ, E, M, VTIMUP, SPINUP);
        if (l >= SIZE - len) {
            return AppendStatus::INTERNALERROR;
        }
        len += l;
    }
    int l = snprintf(buffer + len, SIZE-len, "</event>\n");
    if (l >= SIZE - len) {
        return AppendStatus::INTERNALERROR;
    }
    len += l;

    if (len_ + len >= cap_) { // not enough space in buffer
        double factor = ceil(static_cast<double>(len_ + len - cap_ + 1) /
                             static_cast<double>(grow_cap_));
        cap_ += factor * grow_cap_;
        char *n = (char *)realloc(data_, sizeof(char) * cap_);
        if (n == NULL) {
            // realloc failed
            return AppendStatus::MALLOCERROR;
        }
        data_ = n;
    }

    memcpy(data_ + len_, buffer, len * sizeof(char));
    n_ += 1;
    len_ += len;

    return AppendStatus::SUCCESS;
}

} // end namespace LHE
