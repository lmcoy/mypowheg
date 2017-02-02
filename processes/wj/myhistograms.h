#ifndef MYHISTOGRAMS_H_SDFSR4S4
#define MYHISTOGRAMS_H_SDFSR4S4

#include "process/histograms.h"

class MyHistograms : public Histograms {
  public:
    explicit MyHistograms(int n) : Histograms(n) {
    }
    virtual ~MyHistograms() {
    }

  protected:
    virtual bool Fill(Phasespace::Phasespace const *const ps,
                      double wgt) override {
        const auto &momenta = ps->Momenta;

        if (momenta[5].E() < 1e-12) {
            GetHist(0)->AddValue(momenta[4].PT(), wgt);
        } else {
            double dR = Math::FourMomentum::DeltaR(momenta[4], momenta[5]);
            if (dR < 0.4) {
                // recombine
                auto p = momenta[4].Plus(momenta[5]);
                GetHist(0)->AddValue(p.PT(), wgt);
            } else {
                if (momenta[4].PT() > momenta[5].PT()) {
                    GetHist(0)->AddValue(momenta[4].PT(), wgt);
                } else {
                    GetHist(0)->AddValue(momenta[5].PT(), wgt);
                }
            }
        }
        GetHist(1)->AddValue(momenta[2].PT(), wgt);

        return true;
    }
};

#endif
