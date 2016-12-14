#ifndef FKS_SFUNCTIONS_H_
#define FKS_SFUNCTIONS_H_


namespace Phasespace { class Phasespace; }

namespace FKS { 

class Real_t;
class Region;

double SFunction(const Phasespace::Phasespace &ps, const Region &region,
                 const Real_t &real, int resonance);

} // namespace FKS
#endif
