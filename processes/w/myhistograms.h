#ifndef MYHISTOGRAMS_H_SDFSR4S4
#define MYHISTOGRAMS_H_SDFSR4S4

#include "process/histograms.h"

class MyHistograms : public Histograms {
  public:
    explicit MyHistograms(int n) : Histograms(n) {}
    virtual ~MyHistograms() {}

  protected:
    virtual bool Fill(Phasespace::Phasespace const *const ps,
                      double wgt) override {
        const auto & k = ps->Momenta;
        Math::FourMomentum p_ll = k[2].Plus(k[3]);
        double mll = sqrt(p_ll.Dot(p_ll));
        GetHist(0)->AddValue(mll, wgt);
        GetHist(1)->AddValue(k[2].PT(), wgt);
        GetHist(2)->AddValue(k[3].Rapidity(), wgt);
        return true;
    }
};

#endif
