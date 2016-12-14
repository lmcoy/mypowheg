#include "virtual_qed.h"

/*******************************************************************
 * fortran functions from VirtualQED/                              *
 *******************************************************************/
#ifdef __cplusplus
extern "C" {
#endif

void set_alphas_(double *as_in, double *Q2ren_in);

/**
 * initializes fortran common blocks. call before msquared_
 *
 * implemented in VirtualQED/squaredME.F
 *
 * @param [in] alfa_in alpha
 * @param [in] mw_in W mass
 * @param [in] gw_in W width
 * @param [in] mz_in Z mass
 * @param [in] gz_in Z width
 * @param [in] me_in e mass
 * @param [in] mm_in mu mass
 * @param [in] ml_in tau mass
 * @param [in] mu_in u mass?
 * @param [in] mc_in c mass?
 * @param [in] mt_in t mass?
 * @param [in] md_in d mass?
 * @param [in] ms_in s mass?
 * @param [in] mb_in b mass
 * @param [in] mh_in higgs mass
 * @param [in] gh_in higgs width
 * @param [in] inveps ?
 */
void init_all_(double *alfa_in, double *mw_in, double *gw_in, double *mz_in,
               double *gz_in, double *me_in, double *mm_in, double *ml_in,
               double *mu_in, double *mc_in, double *mt_in, double *md_in,
               double *ms_in, double *mb_in, double *mh_in, double *gh_in,
               double *inveps);

void abbr_const_();

/**
 * msquared_ calculates the squared matrix element
 *
 * implemented in VirtualQED/squaredME.F
 *
 * @param [out] msqu0 born matrix element squared (2D array: { ubar u, dbar d })
 * @param [out] msqu1 virtual 1-loop corrected matrix element squared (2D array)
 * @param [in] k momenta. format double[4][5] where the first index is the four
 * momentum index and the second enumerates the particles (qbar q lbar l empty).
 */
void msquared_(double *msqu0, double *msqu1, double *k, int *flavour,
               int *virtual_type);

void getdzr_(double *res);

#ifdef __cplusplus
}
#endif

void InitVfin_qxq_lxl(const Parameters_sm &param, AlphaScheme scheme,
                      double inveps, double *counterterm) {
    double alpha = param.alpha;
    double mw = param.MW;
    double gw = param.WidthW;
    double mz = param.MZ;
    double gz = param.WZ;
    double me = param.MElectron;
    double mm = param.MMuon;
    double ml = param.MTau;
    double muq = param.MUQuark;
    double mc = param.MCQuark;
    double mt = param.MTQuark;
    double md = param.MDQuark;
    double ms = param.MSQuark;
    double mb = param.MBquark;
    double mh = param.MH;
    double gh = param.WidthH;
    double inveps_in = inveps;

    init_all_(&alpha, &mw, &gw, &mz, &gz, &me, &mm, &ml, &muq, &mc, &mt, &md,
              &ms, &mb, &mh, &gh, &inveps_in);

    abbr_const_();

    
    double coupling = 2.0 * M_PI / param.alpha;

    double cterm = 0.0;
    switch (scheme) {
    case AlphaScheme::Alpha0:
        cterm = 0.0;
        break;
    case AlphaScheme::Gmu:
        double dzr = -1.0;
        getdzr_(&dzr);
        cterm = 2.0 * dzr * coupling;
        break;
    }
    *counterterm = cterm;
}

void Vfin_qxq_lxl(double *born, double *virt, const Math::FourMomentum *momenta,
                  const int *perm, CorrectionType type, int quark_flavour,
                  double alpha, double alphas, double mu2, double counterterm) {

    assert(quark_flavour == 1 || quark_flavour == 2);
    assert(born != NULL);
    assert(virt != NULL);
    
    set_alphas_(&alphas, &mu2);

    int flavour = quark_flavour;
    // copy momenta
    double k[4][5] = { { 0.0 } };
    for (int i = 0; i < 4; ++i) {
        int index = perm[i];
        assert(index >= 0 && index < 5);
        k[0][i] = momenta[index].E();
        k[1][i] = momenta[index].PX();
        k[2][i] = momenta[index].PY();
        k[3][i] = momenta[index].PZ();
    }
    // last momentum is empty!

    double me_born[2] = { 0.0 };
    double me_virtual[2] = { 0.0 };

    int t = 0;
    double charge[2] = {0.0};
    double coupling = 0.0;
    switch (type) {
    case CorrectionType::None:
        t = 0;
        break;
    case CorrectionType::QCD:
        t = 1;
        charge[0] = 2.0 * 4.0 / 3.0; // 2 CF
        charge[1] = charge[0];
        coupling = 2.0 * M_PI / alphas;
        break;
    case CorrectionType::EW:
        charge[1] = 2.0 + 2.0 / 9.0;
        charge[0] = 2.0 + 8.0 / 9.0;
        coupling = 2.0 * M_PI / alpha;
        t = 2;
        break;
    }
    

    msquared_(&me_born[0], &me_virtual[0], &k[0][0], &flavour, &t);
    // avergage over color and spin = 1/36
    if (quark_flavour == 1) {
        // the result is normalized to (4pi)^eps Gamma(1+eps) while
        // my calculation is normalized to (4pi)^eps/Gamma(1-eps).
        // Therefore, there is an additional contribution M_b^2 * pi^2/6 * (sum
        // of squared charges).
        double eps2pole = me_born[1] * charge[1] * M_PI * M_PI / 6.0;
        *born = me_born[1] / 36.0;
        *virt =
            (me_virtual[1] * coupling - eps2pole - me_born[1] * counterterm) /
            36.0;
        return;
    }

    double eps2pole = me_born[0] * charge[0] * M_PI * M_PI / 6.0;
    *born = me_born[0] / 36.0;
    *virt =
        (me_virtual[0] * coupling - eps2pole - me_born[0] * counterterm) / 36.0;
}

