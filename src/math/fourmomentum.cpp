#include "fourmomentum.h"

#include "math/util.h"
#include "math/math.h"
#include "util/stringutil.h"

#include <cmath>

using namespace Math;

/**
At returns a 4-momentum component.

### Invalid indices
- There is an assert which ensures that only valid indices are used.
- If the code is compiled without assert, At returns a quiet NaN for invalid
indices.

\param index The index of the 4-momentum component (0 <= index <= 3).
*/
double FourMomentum::At(const int index) const {
    assert(index >= 0 && index < 4);
    if (index < 0 || index > 3) {
        return std::numeric_limits<double>::quiet_NaN(); // error: index out of
                                                         // range
    }
    return v_[index];
}

/**
 SetE sets the energy of the 4-momentum.

 \param value The new energy of the 4-momentum.
*/
void FourMomentum::SetE(const double value) {
    v_[0] = value;
}

/**
 SetPX sets the x component of the 3-momentum.

 \param value The new x component of the 3-momentum.
*/
void FourMomentum::SetPX(const double value) {
    v_[1] = value;
}

/**
 SetPY sets the y component of the 3-momentum.

 \param value The new y component of the 3-momentum.
*/
void FourMomentum::SetPY(const double value) {
    v_[2] = value;
}

/**
 SetPZ sets the z component of the 3-momentum.

 \param value The new z component of the 3-momentum.
*/
void FourMomentum::SetPZ(const double value) {
    v_[3] = value;
}

/**
SetP sets the momentum of the 4-momentum.

\param px,py,pz The new momentum components.
*/
void FourMomentum::SetP(const double px, const double py, const double pz) {
    v_[1] = px;
    v_[2] = py;
    v_[3] = pz;
}

/**
Set sets the FourMomentum components to (E,px,py,pz).
*/
void FourMomentum::Set(const double E, const double px, const double py,
                       const double pz) {
    v_[0] = E;
    v_[1] = px;
    v_[2] = py;
    v_[3] = pz;
}

/**
Dot returns the dot product of this momentum vector with another momentum
vector.

The dot product is calculated via
\f[
p_1\cdot p_2 = p_1^0 p_2^0 - \vec{p}_1\cdot \vec{p}_2,
\f]
where \f$\vec{p}_1\f$, \f$\vec{p}_2\f$ are the 3-momenta of the 4-momenta and
\f$p_1^0\f$, \f$p_2^0\f$ are the energy components.
*/
double FourMomentum::Dot(const FourMomentum &p) const {
    return v_[0] * p.v_[0] - v_[1] * p.v_[1] - v_[2] * p.v_[2] -
           v_[3] * p.v_[3];
}

/**
Plus returns a FourMomentum = this + p.
*/
FourMomentum FourMomentum::Plus(const FourMomentum &p) const {
    return FourMomentum(v_[0] + p.v_[0], v_[1] + p.v_[1], v_[2] + p.v_[2],
                        v_[3] + p.v_[3]);
}

/**
Minus returns a FourMomentum = this - p.
*/
FourMomentum FourMomentum::Minus(const FourMomentum &p) const {
    return FourMomentum(v_[0] - p.v_[0], v_[1] - p.v_[1], v_[2] - p.v_[2],
                        v_[3] - p.v_[3]);
}

/**
Scale multiplies every component with factor.
*/
void FourMomentum::Scale(const double factor) {
    v_[0] *= factor;
    v_[1] *= factor;
    v_[2] *= factor;
    v_[3] *= factor;
}

void FourMomentum::Add(const FourMomentum &p) {
    v_[0] += p.v_[0];
    v_[1] += p.v_[1];
    v_[2] += p.v_[2];
    v_[3] += p.v_[3];
}

void FourMomentum::Sub(const FourMomentum &p) {
    assert(false);
}

bool FourMomentum::Equals(const FourMomentum &p) const {
    return AlmostEqualDouble(v_[0], p.v_[0]) &&
           AlmostEqualDouble(v_[1], p.v_[1]) &&
           AlmostEqualDouble(v_[2], p.v_[2]) &&
           AlmostEqualDouble(v_[3], p.v_[3]);
}

void FourMomentum::SetFromPMThetaPhi(const double p, const double m,
                                     const double theta, const double phi) {
    const double E = sqrt(m * m + p * p);
#ifdef _GNU_SOURCE
    double sint, cost;
    sincos(theta, &sint, &cost);
    double sinp, cosp;
    sincos(phi, &sinp, &cosp);
#else
    const double cost = cos(theta);
    const double sint = sin(theta);
    const double sinp = sin(phi);
    const double cosp = cos(phi);
#endif
    Set(E, p * sint * cosp, p * sint * sinp, p * cost);
}

void FourMomentum::SetFromPMCosThetaPhi(const double p, const double m,
                                        const double CosTheta,
                                        const double phi) {
    const double E = sqrt(m * m + p * p);
    const double sint = sqrt(1.0 - CosTheta * CosTheta);
    Set(E, p * sint * cos(phi), p * sint * sin(phi), p * CosTheta);
}

FourMomentum FourMomentum::NewFromPMThetaPhi(const double p, const double m,
                                             const double theta,
                                             const double phi) {
    const double E = sqrt(m * m + p * p);
    const double sint = sin(theta);
    return FourMomentum(E, p * sint * cos(phi), p * sint * sin(phi),
                        p * cos(theta));
}

FourMomentum FourMomentum::NewMasslessFromPThetaPhi(const double p,
                                                    const double theta,
                                                    const double phi) {
    const double P = fabs(p);
    const double sint = sin(theta);
    return FourMomentum(P, P * sint * cos(phi), P * sint * sin(phi),
                        P * cos(theta));
}

double FourMomentum::DeltaR(const FourMomentum &p1, const FourMomentum &p2) {
    const double deta = p2.Eta() - p1.Eta();
    double dphi = fabs(p2.Phi() - p1.Phi());

    if (dphi > Math::Pi) {
        dphi = 2.0 * Math::Pi - dphi;
    }

    return sqrt(deta * deta + dphi * dphi);
}

/**
MomentumMagnitude returns the euclidean magnitude of the 3-momentum.
*/
double FourMomentum::MomentumMagnitude() const {
    return sqrt(v_[1] * v_[1] + v_[2] * v_[2] + v_[3] * v_[3]);
}

/**
CosTheta returns the cosine of the polar angle of the 3-momentum.

### Special Cases
- If the 3-momentum is 0.0, CosTheta returns 1.0. (There is an assert which
ensures that this never happens.)
*/
double FourMomentum::CosTheta() const {
    const double mag = MomentumMagnitude();
    assert(mag != 0.0);
    return (mag == 0.0) ? 1.0 : v_[3] / mag;
}

/**
Phi returns the azimuthal angle in the range [0,2Pi).

### Special Cases
- If the 3-momentum points in z direction, Phi returns 0.0.
*/
double FourMomentum::Phi() const {
    if (v_[1] == 0.0 && v_[2] == 0.0) {
        return 0.0;
    }
    const double p1 = atan2(v_[2], v_[1]);
    if (p1 < 0.0) {
        return p1 + 2.0 * Math::Pi;
    }
    return p1;
}

/**
 Eta returns the pseudorapidity -log(tan(theta/2)) of the 3 momentum.

 ### Special Cases
 - If the 3-momenutm points in z direction, Eta returns
 `std::numeric_limits<double>::max()` or
 `std::numeric_limits<double>::lowest()`.
*/
double FourMomentum::Eta() const {
    const double cosTheta = CosTheta();
    // if the particle is in beam direction, eta = +- infinity => return +- max
    // double
    if (cosTheta >= 1.0) {
        return std::numeric_limits<double>::max();
    }
    if (cosTheta <= -1.0) {
        return std::numeric_limits<double>::lowest();
    }
    return -0.5 * log((1.0 - cosTheta) / (1.0 + cosTheta));
}

/**
ToString returns a string representation of this FourMomentum.
*/
std::string FourMomentum::ToString(int digits) const {
    return Strings::Format("(%8.*g, %8.*g, %8.*g, %8.*g)", digits, v_[0],
                           digits, v_[1], digits, v_[2], digits, v_[3]);
}

/**
 * Rapidity returns the rapdity of the four momentum.
 *
 * If the momentum points in z directions, Rapidity returns 
 * `std::numeric_limits<double>::max()` or
 * `std::numeric_limits<double>::lowest()`.
 */
double FourMomentum::Rapidity() const {
    double num = v_[0] + v_[3];
    double denom = v_[0] - v_[3];
    if (num <= 0.0 && num > -1e-15) {
        return std::numeric_limits<double>::lowest();
    }
    if (denom <= 0.0 && denom > -1e-15) {
        return std::numeric_limits<double>::max();
    }
    assert(num > 0.0);
    assert(denom > 0.0);
    return 0.5 * log(num / denom);
}
