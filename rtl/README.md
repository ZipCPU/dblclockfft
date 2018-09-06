This directory contains a demonstration FFT design.

Should you wish to use this core, I would recommend you run `fftgen` from the
[sw](../sw) directory to create an FFT tailored to your own needs.

In sum, from the top down, the modules are:

- [fftmain](fftmain.v) is the top level FFT file.

 - [fftstage](fftstage.v) calculates one FFT stage

  - [hwbfly](hwbfly.v) implements a butterfly that uses the `*` operator
    for its multiply
  - [butterfly](butterfly.v) implements a butterfly that uses a logic
    multiply at the cost of more logic and a greater delay.

   - [longbimpy](longbimpy.v) is the logic binary multiply.

   - [bimpy](bimpy.v) multiplies a small set of bits together.  It is a
     component of [longbimpy](longbimpy.v)

 - [qtrstage](qtrstage.v) is the 4-pt stage of the FFT

 - [laststage](laststage.v) is the 2-pt stage of the FFT

 - [bitreverse](bitreverse.v), the final step in the multiply, bit-reverses
   the outgoing data.


