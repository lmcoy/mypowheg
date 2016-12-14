#ifndef PHYSICS_ALPHAS_H_BQRGG5EY
#define PHYSICS_ALPHAS_H_BQRGG5EY

#include <cmath>
#include <memory>

#include "pdf/pdfinterface.h"

namespace Physics { 

class IAlphaS {
  public:
    virtual double AlphaS(double scale2) const = 0;
    virtual double LambdaQCD() const = 0;
    /**
     * @brief c threshold
     *
     * ThresholdC2 returns the squared flavour activation threshold for c quarks.
     */
    virtual double ThresholdC2() const = 0;
    /**
     * @brief b threshold
     *
     * ThresholdB2 returns the squared flavour activation threshold for b quarks.
     */
    virtual double ThresholdB2() const = 0;
};

class AlphaSfromPDF : public IAlphaS {
    public:
      AlphaSfromPDF(const std::shared_ptr<PDF::Interface> &pdf) : pdf_(pdf) {}
      virtual ~AlphaSfromPDF() {}

    virtual double AlphaS(double scale2) const {
        return pdf_->AlphaS(sqrt(scale2));
    }
    virtual double LambdaQCD() const {
        return pdf_->LambdaQCD();
    }
    virtual double ThresholdC2() const {
        double thc = pdf_->ThresholdC();
        return thc * thc;
    }
    virtual double ThresholdB2() const {
        double thb = pdf_->ThresholdB();
        return thb * thb;
    }

    private:

    std::shared_ptr<PDF::Interface> pdf_;
};

class AlphaSRunning : public IAlphaS {
  public:
    enum class AlphaSOrder {
        LO,
        NLO,
        NNLO
    };
    explicit AlphaSRunning(double alphas_at_mz, AlphaSOrder order, int nf);

    virtual ~AlphaSRunning() {}
    virtual double AlphaS(double scale2) const;
    virtual double LambdaQCD() const {
        return 1.0 / sqrt(invLambda2);
    }

    void SetMc2(double mc2) { M_C2 = mc2; }
    void SetMb2(double mb2) { M_B2 = mb2; }

    virtual double ThresholdC2() const { return M_C2; }
    virtual double ThresholdB2() const { return M_B2; }

  private:
    void LambdaQCD_NL(double as_at_mz, int nf);
    static double alphas_2loop(double scale2, double invLambda2, int nf);
    double invLambda2;
    double M_C2 = 1.5 * 1.5;
    double M_B2 = 5.0 * 5.0;
};

} // end namespace Physics

#endif /* end of include guard: PHYSICS_ALPHAS_H_BQRGG5EY */
