#ifndef PROCESS_CUTS_H_
#define PROCESS_CUTS_H_

#include <memory>

namespace Math {
class FourMomentum;
}

namespace UserProcess {

class ICuts {
  public:
    virtual bool ApplyCuts(int n, const int *pdgs,
                           const Math::FourMomentum *momenta) const = 0;
};

using ICutsPtr = std::shared_ptr<ICuts>;

} // end namespace UserProcess

#endif

