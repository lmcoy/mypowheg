/**\mainpage

All information used in the POWHEG event generation are stored in a  UserProcess::Data object.
UserProcess::Data is an abstract base class. Every process has to derive its own verison.
This makes it possible to make process specific adjustments.
In UserProcess::Data::ProcessInit() you have to set 
- cuts (UserProcess::ICuts)
- matrix elements (UserProcess::IMatrixElement)
- scales (UserProcess::IScales)
- ISR and FSR resonances for the LHE file (UserProcess::Data::ResonanceISR, UserProcess::Data::ResonanceFSR)
- a list of flavour configs (FKS::FlavourConfig)

In principle you have to set the PDFs and AlphaS too, but this is usually done by the POWHEG event generator.
Only when you don't want to use the `powheggen` event generator, you have to set them.

The function FKS::XSecFull() returns the FKS subtracted real cross section. When you integrate over this function
you get the NLO cross section for one specific born sub process (UserProcess::Data::ProcessID).

The POWHEG event generation is implemented in the function Powheg::GenerateEvents(). 
Before you can generate events you have to use Powheg::FindNorm() to find a norm for the upper bounding function
in POWHEG. The event generation is implemented in EventGenerator in `powheggen`.
 

*/

/**
@namespace Powheg @brief POWHEG method

*/