#include <cassert>
#include <iostream>

#include "fks/splitting.h"
#include "physics/pdgcode.h"

namespace {
const double CF = 4.0/3.0;
const double TF = 0.5;
const double CA = 3.0;
}

namespace FKS {

namespace QCD {

/**
 * \internal
 * \brief xi * P(1-xi)
 *
 * splittingTimesXi calculates the QCD splitting function P times xi, i.e.,
 * \f[
 * \xi P_{ab}(1-\xi)
 * \f]
 * where b = real_parton and a = born_parton.
 *
 * \param real_parton pdg of the real parton
 * \param born_parton pdg of the born parton
 * \param xi function parameter
 */
double splittingTimesXi(int real_parton, int born_parton, double xi) {
    if (born_parton == 0 || born_parton == 21) {
        if (real_parton == 0 || real_parton == 21) {
            // gg
            return 2.0 * CA *
                   ((1.0 - xi) + xi * xi / (1.0 - xi) + xi * xi * (1.0 - xi));
        }
        if (real_parton > -6 && real_parton < 6) {
            // qg
            return CF * xi * (1.0 + xi * xi) / (1.0 - xi);
        }
        assert(0 && "real_parton is neither quark nor gluon");
        return -1.0;
    }
    if (born_parton > -6 && born_parton < 6) {
        if (real_parton == 0 || real_parton == 21) {
            // gq
            return TF * xi * (1.0 - 2.0 * xi * (1.0 - xi));
        }
        if (real_parton > -6 && real_parton < 6) {
            assert(born_parton == real_parton &&
                   "no flavour changing splitting");
            // qq
            return CF * (1.0 + (1.0 - xi) * (1.0 - xi));
        }
        assert(0 && "real_parton is neither quark nor gluon");
        return -1.0;
    }
    assert(0 && "born_parton is neither quark nor gluon");
    return -1.0;
}

/**
 * \internal
 * \brief O(eps) of P
 *
 * splittingEps returns the \f$\mathcal{O}(\varepsilon)\f$ term of the
 * QCD splitting function \f$P_{ab}\f$.
 *
 * \param real_parton pdg of the real parton
 * \param born_parton pdg of the born parton
 * \param xi function parameter
 */
double splittingEps(int real_parton, int born_parton, double xi) {
    if (born_parton == 0 || born_parton == 21) {
        if (real_parton == 0 || real_parton == 21) {
            // gg
            return 0.0;
        }
        if (real_parton > -6 && real_parton < 6) {
            // qg
            return -CF * (1.0 - xi);
        }
        assert(0 && "real_parton is neither quark nor gluon");
        return -1.0;
    }
    if (born_parton > -6 && born_parton < 6) {
        if (real_parton == 0 || real_parton == 21) {
            // gq
            return -2.0 * TF * (1.0 - xi) * xi;
        }
        if (real_parton > -6 && real_parton < 6) {
            assert(born_parton == real_parton &&
                   "no flavour changing splitting");
            // qq
            return -CF * xi;
        }
        assert(0 && "real_parton is neither quark nor gluon");
        return -1.0;
    }
    assert(0 && "born_parton is neither quark nor gluon");
    return -1.0;
}

/**
 * \internal
 * \brief xi * P(1-xi) in limit xi = 0
 *
 * splittingTimesXiSoft returns splittinfsTimesXi for xi = 0. 
 *
 * \param real_parton pdg of the real parton
 * \param born_parton pdg of the born parton
 */
double splittingTimesXiSoft(int real_parton, int born_parton) {
    if (born_parton == 0 || born_parton == 21) {
        if (real_parton == 0 || real_parton == 21) {
            return 2.0 * CA;
        }
        if (real_parton > -6 && real_parton < 6) {
            return 0.0;
        }
        assert(0 && "real_parton is neither quark nor gluon");
        return -1.0;
    }
    if (born_parton > -6 && born_parton < 6) {
        if (real_parton == 0 || real_parton == 21) {
            return 0.0;
        }
        if (real_parton > -6 && real_parton < 6) {
            assert(born_parton == real_parton &&
                   "no flavour changing splitting");
            return 2.0 * CF;
        }
        assert(0 && "real_parton is neither quark nor gluon");
        return -1.0;
    }
    assert(0 && "born_parton is neither quark nor gluon");
    return -1.0;
}

} // end namespace QCD

namespace QED {

/**
 * \internal
 * \brief xi * P(1-xi)
 *
 * splittingTimesXi calculates the QED splitting function P times xi, i.e.,
 * \f[
 * \xi P_{ab}(1-\xi)
 * \f]
 * where b = real_pdg and a = born_pdg.
 *
 * \param realpdg pdg of the real particle
 * \param bornpdg pdg of the born particle
 * \param xi function parameter
 */
double splittingTimesXi(int realpdg, int bornpdg, double xi) {
    if (Physics::PDG::IsChargedFermion(realpdg) &&
        Physics::PDG::IsChargedFermion(bornpdg)) {
        assert(realpdg == bornpdg);
        double q = Physics::PDG::Charge(realpdg);

        return q * q * (1.0 + (1.0 - xi) * (1.0 - xi));
    }

    assert(false && "not implemented");
    return -1.0;
}

/**
 * \internal
 * \brief xi * P(1-xi) for xi = 0
 *
 * splittingTimesXi calculates the QED splitting function P times xi, i.e.,
 * \f[
 * \lim_{\xi \rightarrow 0} \xi P_{ab}(1-\xi)
 * \f]
 * where b = real_pdg and a = born_pdg.
 *
 * \param realpdg pdg of the real particle
 * \param bornpdg pdg of the born particle
 */
double splittingTimesXiSoft(int realpdg, int bornpdg) {
    if (Physics::PDG::IsChargedFermion(realpdg) &&
        Physics::PDG::IsChargedFermion(bornpdg)) {
        assert(realpdg == bornpdg);
        double q = Physics::PDG::Charge(bornpdg);
        return 2.0 * q * q;
    }

    assert(false && "not implemented");
    return -1.0;
}

/**
 * \internal
 * \brief O(eps) of P
 *
 * splittingEps returns the \f$\mathcal{O}(\varepsilon)\f$ term of the
 * QED splitting function \f$P_{ab}\f$.
 *
 * \param realpdg pdg of the real particle
 * \param bornpdg pdg of the born particle
 * \param xi function parameter
 */
double splittingEps(int realpdg, int bornpdg, double xi) {
    if (Physics::PDG::IsChargedFermion(realpdg) &&
        Physics::PDG::IsChargedFermion(bornpdg)) {
        assert(realpdg == bornpdg);
        double q = Physics::PDG::Charge(bornpdg);
        return -q * q * xi;
    }

    assert(false && "not implemented");
    return -1.0;
}

} // QED

} // FKS
