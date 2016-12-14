#include "virtual.h"

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
 * @param [in] inveps
 */
void init_all_(double *alfa_in, double *mw_in, double *gw_in, double *mz_in,
               double *gz_in, double *me_in, double *mm_in, double *ml_in,
               double *mu_in, double *mc_in, double *mt_in, double *md_in,
               double *ms_in, double *mb_in, double *mh_in, double *gh_in,
               double *inveps, double *lambdain, int *drin);

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
void msquared_(double *msqu0, double *msqu1, double *k, int *virtual_type);

void getdzr_(double *res);

/**
 * getreg_ returns whether init was called with mass or with dim. reg.
 */
void getreg_(int *reg);

#ifdef __cplusplus
}
#endif

static void init(const Parameters_sm &param, double inveps, double lambda,
                 int dr) {
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
              &ms, &mb, &mh, &gh, &inveps_in, &lambda, &dr);

    abbr_const_();
}

void W_InitMassReg(const Parameters_sm &param, double lambda) {
    init(param, 0.0, lambda, 0);
}

void W_InitDimReg(const Parameters_sm &param, double inveps) {
    init(param, inveps, 1.0, 1);
}

void W_MatrixElement(double *born, double *virt,
                     const Math::FourMomentum *momenta, const int *perm,
                     double alpha, double alphas, double mu, Type type, Scheme scheme) {

    assert(born != NULL);
    assert(virt != NULL);
    double coupling = 1.0;

    double mu2 = mu * mu;
    set_alphas_(&alphas, &mu2);

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

    double me_born = 0.0;
    double me_virtual = 0.0;

    int t = 0; // EW
    double sum_charge2 = 1.0 + 5.0 / 9.0;
    switch (type) {
    case Type::EW:
        t = 2;
        coupling = 2.0 * M_PI / alpha;
        sum_charge2 = 1.0 + 5.0 / 9.0;
        break;
    case Type::QCD:
        t = 1;
        coupling = 2.0 * M_PI / alphas;
        sum_charge2 = 2.0 * 4.0 / 3.0; // 2 CF
        break;
    case Type::None:
	t = 0;
	break;
    default:
        assert(0);
    }

    double counterterm = 0.0;
    if (scheme == Scheme::Gmu) {
        getdzr_(&counterterm);
        counterterm *= 2.0;
    }

    msquared_(&me_born, &me_virtual, &k[0][0], &t);

    // avergage over color and spin = 1/36
    *born = me_born / 36.0;
    *virt = coupling * (me_virtual - counterterm * me_born) / 36.0;

    int dr = 0;
    getreg_(&dr);

    if (dr == 1) { // dim. reg.
        // the result is normalized to (4pi)^eps Gamma(1+eps) while
        // my calculation is normalized to (4pi)^eps/Gamma(1-eps).
        // Therefore, there is an additional contribution M_b^2 * pi^2/6 * (sum
        // of squared charges).
        double eps2pole = me_born * sum_charge2 * M_PI * M_PI / 6.0;
        *virt -= eps2pole / 36.0;
    }
    if (t == 0) {
      *virt = 0.0;
    }
}

