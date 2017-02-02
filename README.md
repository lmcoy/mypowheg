MyPOWHEG
========

POWHEG event generator for QCD and EW corrections.

Requirements
------------

The program requires an MPI implementation. It is only tested with openmpi.
LHAPDF 6 is used for pdf evaluation. There is also an interface to
LHAPDF 5 available but it is not tested extensively. `lhapdf-config` is
used to find your LHAPDF installation. Another requirement is the boost C++
library.

requirement list:
* MPI (openmpi)
* LHAPDF
* Collier (downloaded during compilation)
* boost


Code structure
--------------

The source code consists of five parts:

* _mypowheg_ contains a FKS and POWHEG implemenation. The code is in the src directory.
* _vegas_ contains a VEGAS implementation written in C. It uses MPI to parallelize the integration.
* _powheggen_ is a POWHEG event generator. It uses _vegas_ and _mypowheg_.
* _external_ contains external code for matrix elements etc.
* _processes_ contains POWHEG event generators for processes. Each individual process uses _powheggen_.

Currently three processes are implemented:

* _w_: charged current Drell-Yan: p p -> W+ -> mu+ nu
* _z_: neutral current Drell-Yan: p p -> Z/a -> mu+ mu-
* _wj_: p p -> W+ j -> mu+ nu j **(under development, do not use!)**

Earlier versions of _w_ and _z_ were used to make the plots in
[arXiv:1612.04292](https://inspirehep.net/search?p=find+eprint+1612.04292).

All files with a name ending with "\_test" contain unit tests implemented with
googletest.

Building
--------

CMake is used to generate build files.
To create build files, you have to create a separate build directory.
Run

    mkdir build
    cd build
    cmake PATH_TO_SOURCE

By default the code is compiled optimized and with debugging symbols.
If you use the release mode in cmake (`-DCMAKE_BUILD_TYPE=Release`), the most
aggressive optimizations are used and the program is tuned to the CPU
architecture. This mode is not tested. The debug mode in
cmake (`-DCMAKE_BUILD_TYPE=Debug`) disables all optimizations.

After calling cmake the build directory contains all required build files.
We assume that Makefiles are used. Other options can be chosen with the -G
option of cmake.
You can now use

    make Collier
    make

to compile the code. If you run `make test` all tests are called. The different

Running the code
----------------

POWHEG generator executables are located in the processes directory. You can
run them with

    cd processes/w
    mpirun -np 2 ./w

An example config file `config.txt` is copied to each process directory.

Each process generates the files `vegas_grid.dat` and `run.dat`. These files
contain the vegas state after the Btilde integration.
If the program is called with the `-load` flag, `vegas_grid.dat` and `run.dat`
is loaded and the program can generate events.
Furthermore, the file `norm.txt`, `maxima.txt`, and `vegas_max.txt` are created.
`norm.txt` contains cumulative histograms for the upperbounding function norms.
They can be plotted with `gnuplot`, e.g.
```gnuplot
gnuplot> set terminal pdfcairo
gnuplot> set output 'norm.pdf'
gnuplot> load 'norm.txt'
```
`maxima.txt` contains a histogram for the `Maximum` found during the B~ integration.
It can be used to modify the maximum in `run.dat` in order to make the event generation
more efficient. `vegas_max.txt` contains the maxima found during the vegas integration all
VEGAS bins. The idea is to use those maxima for event generation to make it more efficient.
This feature is not fully implemented so far. One would determine the current VEGAS bin
in every dimension and the maximum. The maximum of all maxima could then be used as
unweighting maximum.


Config files
------------

Every process expects a configuration file `config.txt` in its working directory.
The configuration files contains global options for the event generator and can
also contain process specific options (e.g. cuts, scales etc.).
The config format is inspired by the [ini-format](https://en.wikipedia.org/wiki/INI_file).
The global flags are listed below:

```ini
[recombination]
  dR = 0.00000001 # recombination parameter for NLO calculations. Should be set
                  # tiny values for event generation.

[config]
  SqrtS = 13000.0   # center of mass energy in GeV
  Verbose = 0       # either 0 or 1 (optional, default: 0)
  CalculateXSec = 0 # calculate the cross section separately for born, virtual,
                    # real etc. This calculation is only needed for a precise
                    # cross section determination. (optional, default: 0)

# grid setup for Btilde integration
[integration.setup]
  iterations = 5    # number of VEGAS iterations for the grid setup
  nevents = 1000000 # number of points in each VEGAS iteration
  IgnoreVirtual = 0 # Ignore virtual corrections during grid setup
                    # (optional, default:0)

# Btilde integration
[integration]
  iterations = 5    # number of VEGAS iterations
  nevents = 1000000 # number of points in each VEGAS iteration

# grid setup for cross section calculation. Only used if config.CalculateXSec = 1.
[xsec.setup]
  iterations = 5    # number of VEGAS iterations
  nevents = 1000000 # number of points in each VEGAS iteration

# cross section calculation. Only used if config.CalculateXSec = 1.
[xsec]
  iterations = 5    # number of VEGAS iterations
  nevents = 1000000 # number of points in each VEGAS iteration

# random seeds for the VEGAS integrations
[random]
  seed1 = 1 # random seed for grid setups (0 = random seed)
  seed2 = 2 # random seed for the integration (0 = random seed)

# folding of real part in Btilde (reduces number of negative events, increases
# the run time).
[RadiationSampling]
  RealXi = 1    # folding of xi
  RealY = 1     # folding of y
  RealPhi = 1   # folding of phi
  RemnantXi = 1 # folding of pdf remnant xi

[RadiationType]
  QCD = 1              # enable/disable QCD corrections
  EW = 1               # enable/disable EW corrections
  RadiationQCD = 1     # generate QCD radiation (optional, default: 1)
  RadiationEW = 1      # generate EW radiation (optional, default: 1)
  SplitME = 0          # NC Drell-Yan (z specific): split matrix element in ISR/FSR
  ModS = 0             # use modified FKS S functions. (0: default,
                       # 1: Breit-Wigner improved, 2: Charge & Breit-Wigner improved)
  NoInterference = 0   # NC Drell-Yan (z specific): disable interference terms
  QCDPDFScheme = MSbar # PDF renormalization scheme
  EWPDFScheme = DIS    # PDF renormalization scheme
  # Approximations = NoPhotonRadiation # different approximations. deprecated!

[EventGeneration]
  GenerateEvents = 1    # enable POWHEG event generation (optional, default: 1)
  nevents = 100000      # total number of events
  nperiteration = 10000 # number of events per iteration
  seed = 0              # random seed (0 = random seed)
  NegativeEvents = 0    # use negative events or discard them. (optional, default: 0)
  BornUnweighting = 0   # born unweighting. Do not use! (optional, default: 0)
  BornOnly = 0          # use born instead of B~ (optional, default: 0)

[PDF]
  LHAID = 244600 # select pdfs

[cuts]
  CutReal = 0 # cut on real phase space. Enables real NLO calculation.
              # Do not use for event generation. (optional, default:0)

[Unweighting]
  GuessVirtual = 1 # guess virtual in unweighting, only calculate if necessary
  MaxMultiple = 3  # write events with Btilde > Max multiple times to event file

[Radiation]
  kT2min = 0.8 # minimal radiation kT2
```

NLO Distributions
-----------------

It is possible to generate fixed-order NLO distribution with the code. You have to
enable the option `CutReal` (see above) and you have to implement a process-specific
histogramming class.
The histogramming class has to inherit from the `Histograms` class in `src/process/histograms.h`.
The only method that has to be reimplemented is the `Fill()` function.
A simple example is

```c++
class MyHistograms : public Histograms {
  public:
    explicit MyHistograms(int n) : Histograms(n) {}
    virtual ~MyHistograms() {}

  protected:
    virtual bool Fill(Phasespace::Phasespace const *const ps,
                      double wgt) override {
        const auto &momenta = ps->Momenta;

        // add to histograms
        GetHist(0)->AddValue(momenta[4].PT(), wgt);

        return true;
    }
};
```

An instance of `MyHistograms` has to be used as `hists` variable in `UserProcess::Data`, i.e.

```c++
// use one histogram
hists = new MyHistograms(1);
hists->Init1D(0, 1000, 0.0, 1000.0);
hists->SetName(0, "some description");
```

The histograms are written to `histogram_%d.hist` files, where `%d` is the number of the histogram.

If `CutReal` is not set, only a cut on the underlying born phase space is performed.

Code Documentation
------------------

Doxygen can be used to generate code documentation. Go to the doc directory
and run

    doxygen

to generate `html/index.html`.
