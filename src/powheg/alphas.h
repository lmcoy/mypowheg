#ifndef ALPHAS_H_WH6TICDD
#define ALPHAS_H_WH6TICDD

#include <cmath>

#include "physics/alphas.h"

namespace {

class RadiationAlphaS {
    public:
      RadiationAlphaS(std::shared_ptr<Physics::IAlphaS> alphas)
          : alphas_(alphas) {
          double LambdaQCD = alphas->LambdaQCD();
          double LambdaQCD2x4 = 4.0 * LambdaQCD * LambdaQCD;
          invLambda2_ = exp(1.0 / b0 / AlphaSNLL(LambdaQCD2x4)) / LambdaQCD2x4;
      }

      double SimpleAlphaS(double scale2) const {
          return 1.0 / (b0 * log(scale2 * invLambda2_));
      }

      double AlphaSNLL(double scale2) const {
          double as = alphas_->AlphaS(scale2);
          double lnf = 5.0;
          double M_C2 = alphas_->ThresholdC2();
          double M_B2 = alphas_->ThresholdB2();
          if (scale2 < M_B2 && scale2 > M_C2) {
              lnf = 4.0;
          }
          if (scale2 < M_C2) {
              lnf = 3.0;
          }
          return as * (1.0 + as / (2.0 * M_PI) *
                                 ((67.0 / 18.0 - M_PI * M_PI / 6.0) * 3.0 -
                                  5.0 / 9.0 * lnf));
      }

      double SimpleInvLambda2() const { return invLambda2_; }

      static constexpr double nf = 5.0;
      static constexpr double b0 = (33.0 - 2.0 * nf) / (12.0 * M_PI);

    private:
      std::shared_ptr<Physics::IAlphaS> alphas_;
      double invLambda2_;
};

} // namespace

#endif /* end of include guard: ALPHAS_H_WH6TICDD */ 
