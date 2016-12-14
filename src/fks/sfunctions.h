#ifndef FKS_SFUNCTIONS_H_
#define FKS_SFUNCTIONS_H_


namespace Phasespace { class Phasespace; }

namespace FKS { 

struct Real_t;
struct Region;

/**
 * @brief FKS S function
 *
 * SFunctions returns the value of a FKS S function for the phase space
 * point ps and the singular region "region". It uses all regions in
 * real to compute the S function.
 * If resonance > 0, resonance enhanced regions are used.
 *
 * @param ps        phase space point
 * @param region    singular region
 * @param real      information about the whole process
 * @param resonance If resonance > 0, use enhanced regions. 
 *                  If resonance > 1 and EW rad., use also the charge
 */
double SFunction(const Phasespace::Phasespace &ps, const Region &region,
                 const Real_t &real, int resonance);

} // namespace FKS
#endif
