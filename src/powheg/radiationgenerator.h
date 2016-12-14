#include "util/staticmatrix.h"
#include "phasespace/phasespace.h"
#include "random/rnd.h"
#include "process/data.h"
#include "powheg/generateradiation.h"
#include "powheg/roverb.h"
#include "powheg/pickelement.h"
#include "phasespace/realphasespace.h"
#include "fks/param.h"

namespace {

template <typename T> class RadiationGenerator {
  public:
    explicit RadiationGenerator(double pT2min) : pT2min_(pT2min) {}

    Powheg::RadiationType Generate(const T &functions, int pdf, double B,
                                   const FKS::RadiationRegion &radreg,
                                   const Phasespace::Phasespace &ps,
                                   UserProcess::Data *userdata,
                                   Random::RNG *rng, Powheg::Radiation *rad) {

        Util::StaticMatrix32 roverb(0, 0, 0.0);
        double xi_max = T::XiMax(ps, radreg.Region.J);
        double sb = ps.X1 * ps.X2 * ps.S;
        double kT2tilde = T::PT2Max(xi_max, sb);
        double p2 = kT2tilde;
        double N = radreg.Norm[pdf];
        double minQ = sqrt(userdata->pdf->MinQ2());
        while (true) {
            double pT2 =
                functions.GenPT2(pT2min_, p2, kT2tilde, xi_max, sb, N, rng);
            if (pT2 == 0.0) {
                // born like event
                return Powheg::RadiationType::BORN;
            }
            double xi = 0.0;
            double y = 0.0;
            double phi = 0.0;
            bool suc =
                T::GenRadiationVariables(pT2, sb, xi_max, rng, &xi, &y, &phi);
            if (!suc) {
                fprintf(stderr, "warning: radiation variable out of bounds: xi = %g, y = %g\n", xi, y);
                return Powheg::RadiationType::ERADVAR;
            }
            double alpha = functions.GetUpperboundingCoupling(pT2, userdata);
            double U = T::UpperBounding(N, xi, y, alpha);
            double pdfscale = T::PDFScale(pT2, userdata);
            if (pdfscale < minQ) {
                pdfscale = minQ;
            }

            int jfks = radreg.Region.J;
            assert(jfks >= 0);

            Phasespace::Phasespace ps_real;
            Phasespace::GenRealPhasespace(&ps_real, &ps, jfks, xi, y, phi);

            Util::StaticMatrix32 lumi_ratio(0, 0, 0);
            Powheg::LumiRatio(radreg, ps, ps_real, userdata, pdfscale, pdf,
                              &lumi_ratio);

            bool pdfzero = true;
            for (int i = 0; i < lumi_ratio.Cols(); i++) {
                double L = lumi_ratio.Get(pdf, i);
                if (L < 0.0) {
                    // The born pdfs in the denominator are zero. Therefore, the
                    // propability to radiate would be infinity. This is a weird
                    // situation. We generate a born event and let the
                    // parton shower generate the radiation.
                    return Powheg::RadiationType::BORN;
                }
                if (L != 0.0) {
                    pdfzero = false;
                    break;
                }
            }
            if (pdfzero) {
                // if all pdf values are 0, we can veto the event before
                // evaluating the matrix elements.
                p2 = pT2;
                continue;
            }

            double rad_alpha = functions.GetCoupling(pT2, userdata);
            auto roverb = Powheg::RoverB(B, radreg, ps, ps_real, rad_alpha, userdata);

            double sum = 0.0;
            for (int i = 0; i < lumi_ratio.Cols(); i++) {
                double L = lumi_ratio.Get(pdf, i);
                sum += L * roverb[i];
            }
            if (sum / U > 1.0) {
                fprintf(stderr, "warning: wrong norm for upper bounding "
                                "function: radreg: i = %d j = %d, ",
                        radreg.Region.I, radreg.Region.J);
                fprintf(stderr, "Norm = %g, ", N);
                fprintf(stderr, "(would need Norm >= %g)\n", N * sum / U);
                return Powheg::RadiationType::ENORM;
            }
            if (rng->Random() < sum / U) {
                rad->xi = xi;
                rad->y = y;
                rad->phi = phi;
                rad->kT2 = pT2;
                double prob[8] = { 0.0 };
                for (int i = 0; i < lumi_ratio.Cols(); i++) {
                    prob[i] = lumi_ratio.Get(pdf, i) * roverb[i];
                }
                int rflv = pick_element(lumi_ratio.Cols(), prob, rng->Random());
                rad->j = radreg.Region.J;
                rad->Real = radreg.RealFlavour[rflv];

                return Powheg::RadiationType::REAL;
            }
            p2 = pT2;
        }
        assert(0 && "unreachable");
        return Powheg::RadiationType::ERROR;
    }

  private:
    double pT2min_;
};

} // end namespace
