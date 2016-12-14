#ifndef FKS_SCALES_H_
#define FKS_SCALES_H_

namespace FKS {

struct Scales {
    double muF; ///< factorization scale
    double mu; ///< renorm scale
    double Q2; ///< Arbitrary QES scale
};

} // end namespace FKS

#endif

