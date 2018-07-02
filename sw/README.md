This directory contains the software to generate the FFT.  It compiles into a
program called `fftgen`, which you can then call to generate the FFT you are
interested in.

Components of this coregen include:

- [fftgen.cpp](fftgen.cpp) - This is the top level or 'main' FFT generation program.
- [bldstage.cpp](bldstage.cpp) - Generates the code for a single FFT stage,
  called [fftstage.v](../rtl/fftstage.v) in the RTL directory.
- [softmpy.cpp](softmpy.cpp) - Generates a soft multiply.
- [bitreverse.cpp](bitreverse.cpp) - Generates a bit reverse module
