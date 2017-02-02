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
 * implemented in virtual/squaredME.F
 *
 * @param [out] msqu0 born matrix element squared (5D array)
 * @param [out] msqu1 virtual 1-loop corrected matrix element squared (5D array)
 * @param [in] k momenta. format double[4][6] where the first index is the four
 * momentum index and the second enumerates the particles.
 * @param [in] virtual_type 0: no virtual, 1: QCD, 2: EW
 * @param [in] proc choose process: from 1 to 5 for QCD and 1, 3, 5 for EW.
 */
void msquared_(double *msqu0, double *msqu1, double *k, int *virtual_type,
               int *proc);

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

void Wj_InitDimReg(const Parameters_sm &param, double inveps) {
    init(param, inveps, 1.0, 1);
}

void Wj_InitMassReg(const Parameters_sm &param, double lambda) {
    init(param, 0.0, lambda, 0);
}

void Wj_MatrixElement(int proc, double *born, double *virt,
                      const Math::FourMomentum *momenta, const int *perm,
                      double alpha, double alphas, double mu, Type type,
                      Scheme scheme) {

    assert(born != NULL);
    assert(virt != NULL);
    assert(proc >= 0 && proc <= 5);
    double coupling = 1.0;

    double mu2 = mu * mu;
    set_alphas_(&alphas, &mu2);

    // copy momenta
    double k[4][6] = {{0.0}};
    for (int i = 0; i < 5; ++i) {
        int index = perm[i];
        assert(index >= 0 && index < 6);
        k[0][i] = momenta[index].E();
        k[1][i] = momenta[index].PX();
        k[2][i] = momenta[index].PY();
        k[3][i] = momenta[index].PZ();
    }
    k[0][5] = 0.0;
    k[1][5] = 0.0;
    k[2][5] = 0.0;
    k[3][5] = 0.0;
    // last momentum is empty!

    // 5 sub processes:
    // 0: u d~ -> nu mu+ g
    // 1: u a -> nu mu+ d
    // 2: u g -> nu mu+ d
    // 3: d~ a -> nu mu+ u~
    // 4: d~ g -> nu mu+ u~
    double me_born[31] = {0.0};
    double me_virtual[31] = {0.0};

    // average over spins and colors
    double denominators[5] = {36.0, 96.0, 96.0, 96.0, 96.0};

    int t = 0;
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
        sum_charge2 = 2.0 * 4.0 / 3.0 + 3.0; // 2 CF + CA
        break;
    case Type::None:
        t = 0;
        break;
    default:
        assert(0);
    }

    double counterterm = 0.0;
    if (type == Type::EW && scheme == Scheme::Gmu) {
        getdzr_(&counterterm);
        counterterm *= 2.0;
    }
    int fproc = proc + 1;
    msquared_(&me_born[0], &me_virtual[0], &k[0][0], &t, &fproc);

    // average over color and spin
    *born = me_born[proc] / denominators[proc];
    *virt = coupling * (me_virtual[proc] - counterterm * me_born[proc]) /
            denominators[proc];

    int dr = 0;
    getreg_(&dr);

    if (dr == 1) { // dim. reg.
        // the result is normalized to (4pi)^eps Gamma(1+eps) while
        // my calculation is normalized to (4pi)^eps/Gamma(1-eps).
        // Therefore, there is an additional contribution M_b^2 * pi^2/6 * (sum
        // of squared charges).
        double eps2pole = me_born[proc] * sum_charge2 * M_PI * M_PI / 6.0;
        *virt -= eps2pole / denominators[proc];
    }
    if (t == 0) {
        *virt = 0.0;
    }
}
