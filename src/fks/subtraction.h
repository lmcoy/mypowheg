#ifndef FKS_SUBTRACTION_H_
#define FKS_SUBTRACTION_H_

namespace FKS { 

class Xi;
class MatrixElement;
class PartonLuminosity;

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
 */
SubtractionTerms RealFSR(double lumi, const MatrixElement &me, double x, double ximax,
               double y, double len_kj_born, double s);

/**
 * @brief subtracted real matrix element for ISR
 *
 */
SubtractionTerms RealISR2(const PartonLuminosity &lum, const MatrixElement &me,
                int colldir, double x, const Xi &xi, double y);

} // namespace FKS

#endif
