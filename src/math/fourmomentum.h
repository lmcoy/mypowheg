#ifndef _FOURMOMENTUM_H_
#define _FOURMOMENTUM_H_

#include <cassert>
#include <limits>
#include <array>
#include <cmath>

namespace Math {

/**
 \brief 4-momentum
 \ingroup Math
 */
class FourMomentum {
  public:
    /** @name Constructors */
    ///@{
    FourMomentum() : v_{ 0.0, 0.0, 0.0, 0.0 } {
    }
    FourMomentum(double E, double px, double py, double pz)
        : v_{ E, px, py, pz } {
    }
    FourMomentum(const std::array<double, 4> &p)
        : v_{ p[0], p[1], p[2], p[3] } {
    }
    ///@}

    /** @name Getter */
    ///@{
    double E() const {
        return v_[0];
    }
    double PX() const {
        return v_[1];
    }
    double PY() const {
        return v_[2];
    }
    double PZ() const {
        return v_[3];
    }
    double At(const int index) const;
    ///@}

    /** @name Setter */
    ///@{
    void SetE(const double value);
    void SetPX(const double value);
    void SetPY(const double value);
    void SetPZ(const double value);
    void SetP(const double px, const double py, const double pz);
    void Set(const double E, const double px, const double py, const double pz);
    void SetFromPMThetaPhi(const double p, const double m, const double theta,
                           const double phi);
    void SetFromPMCosThetaPhi(const double p, const double m,
                              const double CosTheta, const double phi);
    ///@}

    /** @name 3-momentum functions */
    ///@{
    double MomentumMagnitude() const;
    double MomentumMagnitudeSqr() const { 
        return v_[1] * v_[1] + v_[2] * v_[2] + v_[3] * v_[3];
    }
    double CosTheta() const;
    double Phi() const;
    double Eta() const;
    double Rapidity() const;
    double PT_Sqr() const {
        return v_[1]*v_[1] + v_[2]*v_[2];
    }
    double PT() const { 
        return sqrt(PT_Sqr());
    }
    ///@}

    /** @name math operations */
    ///@{
    double Dot(const FourMomentum &p) const;
    FourMomentum Plus(const FourMomentum &p) const;
    FourMomentum Minus(const FourMomentum &p) const;
    ///@}

    /** @name mathematical modifiers */
    ///@{
    void Scale(const double factor);
    void Add(const FourMomentum &p);
    void Sub(const FourMomentum &p);
    ///@}

    /** @name Comparison */
    ///@{
    bool Equals(const FourMomentum &p) const;
    ///@}

    /** @name string methods */
    ///@{
    std::string ToString(int digits = 6) const;
    ///@}

    /** @name Create new 4-momenta */
    ///@{
    static FourMomentum NewFromPMThetaPhi(const double p, const double m,
                                          const double theta, const double phi);
    static FourMomentum NewMasslessFromPThetaPhi(const double p,
                                                 const double theta,
                                                 const double phi);
    ///@}

    static double DeltaR(const FourMomentum &p1, const FourMomentum &p2);

  private:
    double v_[4]; /**< data */
};

} // namespace Math

#endif
