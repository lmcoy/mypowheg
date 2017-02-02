#ifndef FKS_PDF_H_
#define FKS_PDF_H_

#include <cassert>
#include <cmath>
#include <memory>

#include "libconfig.h"

#define LIST

#ifdef MAP
#include <map>
#endif

#include "pdf/pdfinterface.h"

namespace PDF {

class LIB_LOCAL Cache {
  public:
    Cache(std::shared_ptr<PDF::Interface> p) : pdfi(p) {}
#ifdef MAP
    double Get(double x, double muF, int flavour) {
        assert(fabs(muF - _muF) < 1e-9 || _muF < 0.0);
        assert(flavour >= -5);
        assert(flavour <= 5);
        _muF = muF;
        int index = flavour + 5;

        auto it = _cache[index].find(x);
        if (it == _cache[index].end()) {
            double f = pdfi->Xfx(x, muF, flavour) / x;
            _cache[index][x] = f;
            return f;
        }
        return it->second;
    }
#elif defined(LIST)
    double Get(double x, double muF, int flavour) {
        assert(fabs(muF - _muF) < 1e-9 || _muF < 0.0);
        if (flavour == 21) {
            flavour = 0;
        }
        assert(flavour >= -5);
        assert(flavour <= 5);
        _muF = muF;
        int index = flavour + 5;
        size_t i = 0;
        size_t len = _len[index];
        for (i = 0; i < len && i < 20; i++) {
            if (fabs(_cache2[index][i].x - x) < 1e-15) {
                return _cache2[index][i].y;
            }
        }
        double f = pdfi->Xfx(x, muF, flavour) / x;
        if (len < 20) {
            _cache2[index][len].x = x;
            _cache2[index][len].y = f;
            _len[index] += 1;
        }
        return f;
    }

    size_t Len(int flavour) const { return _len[flavour + 5]; }

#else
    double Get(double x, double muF, int flavour) {
        return pdfi->Xfx(x, muF, flavour) / x;
    }
#endif

    std::shared_ptr<PDF::Interface> Interface() { return pdfi; }

  private:
    double _muF = -1.0;
    std::shared_ptr<PDF::Interface> pdfi;

#ifdef MAP
    std::map<double, double> _cache[11];
#elif defined(LIST)
    struct item {
        double x = -1.0;
        double y = -1.0;
    };
    size_t _len[11] = {0};
    item _cache2[11][20];
#endif
};

} // end namespace FKS
#endif
