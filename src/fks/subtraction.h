#ifndef FKS_SUBTRACTION_H_
#define FKS_SUBTRACTION_H_

namespace FKS { 

struct Xi;
struct MatrixElement;
struct PartonLuminosity;

/** @brief Subtracted real matrix element */
struct SubtractionTerms {
    double Real = 0.0;
    double Soft = 0.0;
    double Collinear = 0.0;
    double SoftCollinear = 0.0;
};

/**
 * @brief subtracted real matrix element for FSR
 *
 * RealFSR returns the real subtracted part R of the FKS subtraction.
 * It has to be integrated over the born phase-space and three parameters, i.e.
 * \int dPhi_n \int_0^1 x \int_-1^1 y \int_0^2pi phi R
 * (Note: The jacobian for the x,y,phi integration is included. This function
 * uses the jacobian for the massless real phase space from 0709.2092 eq. (5.41).
 *
 * @param lumi        luminosity for the real phase space point
 * @param me          matrix element and its limits
 * @param x           integration parameter in [0,1]
 * @param ximax       maximal energy fraction. It is used to calculate the energy
 *                    fraction xi = x * ximax
 * @param y           usual y integragtion variabel (=cos theta)
 * @param len_kj_born length of born momentum of radiating particle
 * @param s           center of mass energy of the real phase space
 */
SubtractionTerms RealFSR(double lumi, const MatrixElement &me, double x, double ximax,
               double y, double len_kj_born, double s);

/**
 * @brief subtracted real matrix element for ISR
 *
 * @see RealFSR()
 *
 * @param lum     luminosity
 * @param me      matrix element and its limits
 * @param colldir either -1 or 1. choose between 1/(1-y) and 1/(1+y) plus distr.
 * @param x       integration variable in [0,1]
 * @param xi      maximal energy fraction
 * @param y       y integration variable (=cos theta)
 */
SubtractionTerms RealISR2(const PartonLuminosity &lum, const MatrixElement &me,
                int colldir, double x, const Xi &xi, double y);

} // namespace FKS

#endif
