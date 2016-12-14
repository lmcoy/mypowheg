#ifndef LHAPDF6_H_I9SC013K
#define LHAPDF6_H_I9SC013K

#include <string>

#include "LHAPDF/LHAPDF.h"

#include "pdf/pdfinterface.h"

namespace PDF {

#ifdef LHAPDF_MAJOR_VERSION

class Lhapdf : public PDF::Interface {
  public:
    virtual ~Lhapdf() {
        if (pdf) {
            delete pdf;
        }
    }

    void InitByLHAID(int id) {
        assert(!pdf);
        pdf = LHAPDF::mkPDF(id);
        mc = pdf->quarkThreshold(4);
        mb = pdf->quarkThreshold(5);
        q2min = pdf->q2Min();
        as_order = pdf->orderQCD();
    }

    virtual double Xfx(double x, double muF, int flavour) const {
        int fl_abs = (flavour >= 0) ? flavour : -flavour;
        if (fl_abs == 5 && muF < mb) {
            return 0.0;
        }

        if (fl_abs == 4 && muF < mc) {
            return 0.0;
        }
        double xfx = pdf->xfxQ(flavour, x, muF);
        if (xfx < 0.0) {
            return 0.0;
        }
        return xfx;
    }

    virtual std::array<double, 13> XfxAll(double x, double muF) const {
        std::vector<double> xfx;
        xfx.reserve(13);
        pdf->xfxQ(x, muF, xfx);
        if (muF < mb) {
            xfx[11] = 0.0; // set b pdf to 0.0
            xfx[1] = 0.0;  // set bbar pdf to 0.0
        }
        if (muF < mc) {
            xfx[10] = 0.0; // set c pdf to 0.0
            xfx[2] = 0.0;  // set cbar pdf to 0.0
        }
        std::array<double, 13> result;
        for (int i = 0; i < 13; i++) {
            if (xfx[i] >= 0.0) {
                result[i] = xfx[i];
            } else {
                result[i] = 0.0;
            }
        }
        return result;
    }

    virtual double LambdaQCD() const {
        assert(0 && "not implemented in LHAPDF 6");
        return 0.2;
    }

    virtual double AlphaS(double Q) const {
        return pdf->alphasQ(Q);
    }

    virtual double ThresholdC() const {
        return mc;
    }

    virtual double ThresholdB() const {
        return mb;
    }

    virtual void SetThresholdC(double thc) { mc = thc; }

    virtual void SetThresholdB(double thb) { mb = thb; }

    virtual int GetOrderAlphaS() const {
        return as_order;
    }

    virtual double MinQ2() const {
        return q2min;
    }
  private:
    LHAPDF::PDF *pdf = 0;
    double mc = 0.0;
    double mb = 0.0;
    double q2min = 0.0;
    int as_order = 0;
};

#else  /* use lhapdf5 interface */

class Lhapdf : public PDF::Interface {
  public:
    virtual ~Lhapdf() {}
    void InitByLHAID(int id) {
        LHAPDF::setVerbosity(LHAPDF::SILENT);
        LHAPDF::initPDFSet(id, 0);
        xmin = LHAPDF::getXmin(0);
        xmax = LHAPDF::getXmax(0);
        q2min = LHAPDF::getQ2min(0);
        mc = LHAPDF::getThreshold(4);
        mb = LHAPDF::getThreshold(5);
        as_order = LHAPDF::getOrderAlphaS();
        LHAPDF::initPDFSet(id, 0);
    }
    virtual double Xfx(double x, double muF, int flavour) const {
        if (std::abs(flavour) == 5 && muF < mb) {
            return 0.0;
        }

        if (std::abs(flavour) == 4 && muF < mc) {
            return 0.0;
        }

        double xfx = LHAPDF::xfx(x, muF, flavour);
        if (xfx < 0.0) {
            return 0.0;
        }
        return xfx;
    }

    virtual std::array<double, 13> XfxAll(double x, double muF) const {
        double xfx[13] = { 0.0 };
        LHAPDF::xfx(x, muF, xfx);
        if (muF < mb) {
            xfx[11] = 0.0; // set b pdf to 0.0
            xfx[1] = 0.0;  // set bbar pdf to 0.0
        }
        if (muF < mc) {
            xfx[10] = 0.0; // set c pdf to 0.0
            xfx[2] = 0.0;  // set cbar pdf to 0.0
        }
        std::array<double, 13> result;
        for (int i = 0; i < 13; i++) {
            if (xfx[i] >= 0.0) {
                result[i] = xfx[i];
            } else {
                result[i] = 0.0;
            }
        }
        return result;
    }

    virtual double AlphaS(double Q) const { return LHAPDF::alphasPDF(Q); }

    virtual double LambdaQCD() const { return LHAPDF::getLam5(0); }
    
    virtual double ThresholdC() const {
        return mc;
    }

    virtual double ThresholdB() const {
        return mb;
    }

    virtual void SetThresholdC(double thc) { mc = thc; }

    virtual void SetThresholdB(double thb) { mb = thb; }

    virtual int GetOrderAlphaS() const {
        return as_order;
    }

    virtual double MinQ2() const {
        return q2min;
    }
  private:
    double mc = 0.0;
    double mb = 0.0;
    double xmin = 0.0;
    double xmax = 1.0;
    double q2min = 0.0;
    int as_order = 0;
};

#endif

} // namespace PDF

#endif /* end of include guard: LHAPDF6_H_I9SC013K */ 
