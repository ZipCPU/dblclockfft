This directory contains several SymbiYosys scripts useful for
formally verifying parts and pieces of the design.  Admittedly,
the entire design has yet to be formally verfified, however many
components have been verified successfully.  These include:

- The butterflies, both the hardware enabled butterflies and the
  soft multiplies.

- The penultimate (4-pt) stage of the FFT

- The final stage (2-pt) of the FFT

- The bitreverse

My intention is not to place formal properties into the repository.
Within the [defaults.h](../../sw/defaults.h) there's a
``formal_property_flag`` used for controlling whether or not the
formal properties are included into the RTL files.
