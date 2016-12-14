#ifndef FKS_PARAM_H_
#define FKS_PARAM_H_

namespace FKS {

class Param {
  public:
    virtual ~Param() {}
    virtual double Mass(int pdg) const = 0;

    virtual void SetAlphaS(double as) = 0;

    virtual double alphaS() const = 0;

    double alpha;
};

} // namespace FKS

#endif
