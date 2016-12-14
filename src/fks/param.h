#ifndef FKS_PARAM_H_
#define FKS_PARAM_H_

namespace FKS {

class Param {
  public:
    virtual ~Param() {}
    virtual double Mass(int pdg) const = 0;
    double alpha;
};

class Param_as {
  public:
    virtual ~Param_as() {}
    virtual void Set(double as) = 0;
    double aS;
};

} // namespace FKS

#endif
