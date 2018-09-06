This directory contains several SymbiYosys scripts useful for formally verifying parts and pieces
of the design.  While the top level of the design has yet to be formally verified, scripts are
present for verifying all of the other components--save the [long binary multiply,
longbimpy](../../rtl/longbimpy.v).  The properties within the [long binary
multiply](../../rtl/longbimpy.v) only proves certain properties of the binary multiply.  The [full
proof](../cpp/mpy_tb.cpp) is done by exhaustion via Verilator in the [bench/cpp](../cpp) directory.

Within the [defaults.h](../../sw/defaults.h) there's a ``formal_property_flag`` used for
controlling whether or not the formal properties are included into the RTL files.
