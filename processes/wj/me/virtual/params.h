      integer np,nproc,nproctree,nprocbrems,nres,npole,ncuts,nhel
      integer nheltree,ncol,maxnminvc,maxnmtransc,maxv,maxvphys,maxch
      integer maxg,maxo,maxh,ndipolereal,ndipole,nndipole,nhelpairs

* Parameters for the Process Definition:

* maximal Number of partonic subprocesses per generator
      parameter(nproc = 31)
* Number of partonic subprocesses at tree level
      parameter(nproctree = 5)
* Number of partonic subprocesses for bremsstrahlung
      parameter(nprocbrems = 31)
* Maximal number of particles for all considered processes
      parameter(np = 6)
* Maximal number of resonances for all considered processes
      parameter(nres = 0)
* Maximal number of resonant propagators
      parameter(npole = 1)
* Maximal number of helicities for all considered processes
      parameter(nhel = 11)
* number of helicities for tree processes
      parameter(nheltree = 5)
* Maximal number of color structures for all considered processes
      parameter(ncol = 9)
* Number of explicit dipole contributions for real subtraction
      parameter(ndipolereal = 112)
* Number of explicit dipole contributions for convolution
      parameter(ndipole = 60)
* Number of dipoles for specific tree process, emitter and spectator
      parameter(nndipole = 9)
* Number of pairs of Helicity combinations with opposite helicity gluons
      parameter(nhelpairs = 2)

* Parameters for the MC Integration:

* maximal number of particles in the considered model
      parameter(maxv = 24)
* maximal number of physical particles in the considered model
      parameter(maxvphys = 20)
* maximal number of channels
      parameter(maxch = 2000)
* number of generators
      parameter(maxg = 2)
* maximal number of optimization steps
      parameter(maxo = 20)
* maximal number of histograms
      parameter(maxh = 25)
* Maximal number of cuts
      parameter(ncuts = 30)
