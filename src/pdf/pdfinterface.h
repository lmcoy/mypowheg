#ifndef PDFINTERFACE_H_IW8ST67G
#define PDFINTERFACE_H_IW8ST67G

#include <memory>
#include <array>

namespace PDF {

class Interface {
  public:
    /**
     * @brief pdf value
     *
     * Xfx returns x*f(x, muF, flavour) where muF is the factorization scale.
     *
     * Remarks:
     * - The fitted pdfs can some times return negative values for large values
     *   of x. This function should always return 0.0 for those cases.
     * - This function should return 0.0 if |flavour| = 5 and 
     *   muF < ThresholdB().
     *   (the same for |flavour| = 4 and muF < ThresholdC() )
     *
     */
    virtual double Xfx(double x, double muF, int flavour) const = 0;

    /**
     * @brief pdf values
     *
     * XfxAll returns x*f(x, muF) for all flavours.
     *
     * The numbering is the same as in LHAPDF, i.e.
     *  0..5 : tbar, ..., ubar, dbar
     *     6 : gluon
     * 7..12 : d, u, ..., t
     *
     * The same remarks as for Xfx() apply.
     */
    virtual std::array<double, 13> XfxAll(double x, double muF) const = 0;

    /**
     * @brief alphas from pdf
     *
     * AlphaS returns the alphas value at scale Q from the pdf.
     */
    virtual double AlphaS(double Q) const = 0;

    virtual double LambdaQCD() const = 0;

    /**
     * @brief order of alphas
     *
     * GetOrderAlphaS returns the order of alphas which is returned by AlphaS().
     */
    virtual int GetOrderAlphaS() const = 0;

    /**
     * @brief C quark threshold
     *
     * ThresholdC returns the c quark flavour threshold.
     * This value should be initially set by the pdf and can be changed with
     * SetThresholdC().
     */
    virtual double ThresholdC() const = 0;

    /**
     * @brief B quark threshold
     *
     * ThresholdB returns the b quark flavour threshold.
     * This value should be initially set by the pdf and can be changed with
     * SetThresholdB().
     */
    virtual double ThresholdB() const = 0;

    /**
     * @brief Set c quark threshold
     *
     * See ThresholdC().
     */
    virtual void SetThresholdC(double thc) = 0;

    /**
     * @brief Set b quark threshold
     *
     * See ThresholdB().
     */
    virtual void SetThresholdB(double thb) = 0;
    /**
     * @brief minimal factorization scale
     *
     * MinQ2 returns the squared minimal factorization scale valid for this pdf.
     */
    virtual double MinQ2() const = 0;
};

typedef std::shared_ptr<PDF::Interface> InterfacePtr;

} // end namespace PDF

#endif /* end of include guard: PDFINTERFACE_H_IW8ST67G */ 
