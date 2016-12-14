#ifndef POWHEG_BTILDE_H_ 
#define POWHEG_BTILDE_H_ 

#include <cstdlib>
#include <utility>
#include <vector>
#include <cassert>

#include "libconfig.h"
#include "powheg/pickelement.h"

namespace UserProcess {
class Data;
}

namespace Phasespace {
class Phasespace;
}

namespace Powheg {

LIB_PUBLIC int Btilde(const Phasespace::Phasespace &, double, double, double, double,
           double *, UserProcess::Data *);

class Btilde_t {
  public:
    typedef std::pair<size_t, size_t> Meta;

    void CalcWOVirtual(const Phasespace::Phasespace &ps, double x1, double x2,
                       double x3, double wgt, UserProcess::Data *params);

    void CalcVirtual(const Phasespace::Phasespace &ps, double x1, double x2,
                     double x3, double wgt, UserProcess::Data *params);

    double TotalWithGuessedVirtual() const {
        return sum_guessed;
    }

    double Total() const {
        assert(has_virtual);
        if (!has_virtual) {
            return 0.0;
        }
        return sum;
    }

    Meta PickElement(double rnd) {
        int n = data.size();
        int index = pick_element(n, data.data(), rnd);
        return meta[index];
    }
            
  private:
    void Append(size_t processIndex, size_t pdfIndex, double btilde, double btildeg) {
        meta.push_back(Meta(processIndex, pdfIndex));
        data.push_back(btilde);
        sum += btilde;
        sum_guessed += btildeg;
    }
    double sum = 0.0;
    double sum_guessed = 0.0;
    std::vector<Meta> meta;
    std::vector<double> data;
    bool has_virtual = false;
};

} // end namespace Powheg

#endif
