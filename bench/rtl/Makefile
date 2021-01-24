all: ifft_tb
################################################################################
##
## Filename: 	Makefile
## {{{
## Project:	A General Purpose Pipelined FFT Implementation
##
## Purpose:	This programs the build process for part of the ifft_tb test
##		bench associated with the double clocked FFT project.
##
##		This is only part of the test bench program.  This one is
##		different from the others, in that the ifft_tb includes
##		both verilator code in the test bench as well as the
##		C++ code.  The C++ code will depend upon the verilog
##		code found in this directory and built here.
##
## Creator:	Dan Gisselquist, Ph.D.
##		Gisselquist Technology, LLC
##
################################################################################
## }}}
## Copyright (C) 2015-2021, Gisselquist Technology, LLC
## {{{
## This program is free software (firmware): you can redistribute it and/or
## modify it under the terms of  the GNU General Public License as published
## by the Free Software Foundation, either version 3 of the License, or (at
## your option) any later version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
## for more details.
##
## You should have received a copy of the GNU General Public License along
## with this program.  (It's in the $(ROOT)/doc directory.  Run make with no
## target there if the PDF file isn't present.)  If not, see
## <http://www.gnu.org/licenses/> for a copy.
## }}}
## License:	GPL, v3, as defined and found on www.gnu.org,
## {{{
##		http://www.gnu.org/licenses/gpl.html
##
##
################################################################################
##
## }}}
# This is really simple ...
all: ifft_tb
.DELETE_ON_ERROR:
ifneq ($(VERILATOR_ROOT),)
VERILATOR:=$(VERILATOR_ROOT)/bin/verilator
else
VERILATOR:=verilator
VERILATOR_ROOT ?= $(shell bash -c 'verilator -V|grep VERILATOR_ROOT | head -1 | sed -e " s/^.*=\s*//"')
endif
export	$(VERILATOR)
VROOT   := $(VERILATOR_ROOT)
VSRCD   := ../../rtl
VFLAGS  := -Wall -MMD --trace -y $(VSRCD) -cc

LCLDR:= obj_dir
IFTLB:= $(LCLDR)/Vifft_tb__ALL.a

.PHONY: ifft_tb
## {{{
ifft_tb: $(IFTLB)

$(IFTLB): $(LCLDR)/Vifft_tb.cpp
	cd $(LCLDR); make -f Vifft_tb.mk
$(LCLDR)/Vifft_tb.cpp: ifft_tb.v $(VSRCD)/fftmain.v $(VSRCD)/ifftmain.v
	$(VERILATOR) $(VFLAGS) ifft_tb.v
## }}}

.PHONY: clean
## {{{
clean:
	rm -rf $(LCLDR)
## }}}

-include $(VERILATOR_ROOT)/include/verilated.mk
