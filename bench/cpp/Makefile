################################################################################
##
## Filename: 	Makefile
## {{{
## Project:	A General Purpose Pipelined FFT Implementation
##
## Purpose:	This programs the build process for the test benches
##		associated with the double clocked FFT project.  These
##	test benches are designed for the size and arguments of the FFT as
##	given by the Makefile in the trunk/sw directory, although they shouldn't
##	be too difficult to modify for other FFT parameters.
##
##	Please note that running these test benches requires access to the
##	*cmem_*.hex files found in trunk/rtl.  They will be linked into this
##	directory as part of the tests.
##
## Creator:	Dan Gisselquist, Ph.D.
##		Gisselquist Technology, LLC
##
################################################################################
## }}}
## Copyright (C) 2015-2024 Gisselquist Technology, LLC
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
################################################################################
## }}}
all: mpy_tb bitreverse_tb hwbfly_tb butterfly_tb fftstage_tb fft_tb
all: qtrstage_tb laststage_tb

OBJDR:= ../../rtl/obj_dir
VSRCD = ../../rtl
TBODR:= ../rtl/obj_dir
ifneq ($(VERILATOR_ROOT),)
VERILATOR:=$(VERILATOR_ROOT)/bin/verilator
else
VERILATOR:=verilator
VERILATOR_ROOT ?= $(shell bash -c 'verilator -V|grep VERILATOR_ROOT | head -1 | sed -e " s/^.*=\s*//"')
endif
export	$(VERILATOR)
VROOT   := $(VERILATOR_ROOT)
VDEFS:= $(shell ./vversion.sh)
VINC := -I$(VROOT)/include -I$(OBJDR)/ -I$(TBODR)/
# MPYLB:= $(OBJDR)/Vshiftaddmpy__ALL.a
MPYLB:= $(OBJDR)/Vlongbimpy__ALL.a
BTREV:= $(OBJDR)/Vbitreverse__ALL.a
STAGE:= $(OBJDR)/Vfftstage__ALL.a
QTRSG:= $(OBJDR)/Vqtrstage__ALL.a
LSTSG:= $(OBJDR)/Vlaststage__ALL.a
BFLYL:= $(OBJDR)/Vbutterfly__ALL.a
HWBFY:= $(OBJDR)/Vhwbfly__ALL.a
FFTLB:= $(OBJDR)/Vfftmain__ALL.a
IFTLB:= $(TBODR)/Vifft_tb__ALL.a
STGLB:= $(OBJDR)/Vfftstage__ALL.a
VSRCS:= $(VROOT)/include/verilated.cpp $(VROOT)/include/verilated_vcd_c.cpp $(VROOT)/include/verilated_threads.cpp

mpy_tb: mpy_tb.cpp fftsize.h twoc.h $(MPYLB)
	g++ -g $(VINC) $(VDEFS) $< twoc.cpp $(MPYLB) $(VSRCS) -lpthread -o $@

bitreverse_tb: bitreverse_tb.cpp twoc.cpp twoc.h fftsize.h $(BTREV)
	g++ -g $(VINC) $(VDEFS) $< twoc.cpp $(BTREV) $(VSRCS) -lpthread -o $@

laststage_tb: laststage_tb.cpp twoc.cpp twoc.h $(LSTSG)
	g++ -g $(VINC) $(VDEFS) $< twoc.cpp $(LSTSG) $(VSRCS) -lpthread -o $@

qtrstage_tb: qtrstage_tb.cpp twoc.cpp twoc.h $(QTRSG)
	g++ -g $(VINC) $(VDEFS) $< twoc.cpp $(QTRSG) $(VSRCS) -lpthread -o $@

butterfly_tb: butterfly_tb.cpp twoc.cpp twoc.h fftsize.h $(BFLYL)
	g++ -g $(VINC) $(VDEFS) $< twoc.cpp $(BFLYL) $(VSRCS) -lpthread -o $@

hwbfly_tb: hwbfly_tb.cpp twoc.cpp twoc.h $(HWBFY)
	g++ -g $(VINC) $(VDEFS) $< twoc.cpp $(HWBFY) $(VSRCS) -lpthread -o $@

fftstage_tb: fftstage_tb.cpp twoc.cpp twoc.h $(STGLB)
	g++ -g $(VINC) $(VDEFS) $< twoc.cpp $(STGLB) $(VSRCS) -lpthread -o $@

fft_tb: fft_tb.cpp twoc.cpp twoc.h fftsize.h $(FFTLB)
	g++ -g $(VINC) $(VDEFS) $< twoc.cpp $(FFTLB) $(VSRCS) -lpthread -lfftw3 -o $@

ifft_tb: ifft_tb.cpp twoc.cpp twoc.h fftsize.h $(IFTLB)
	g++ -g $(VINC) $(VDEFS) $< twoc.cpp $(IFTLB) $(VSRCS) -lpthread -lfftw3 -o $@

.PHONY: HEX
HEX:
	ln -sf $(VSRCD)/*.hex .

.PHONY: test
test: mpy_tb.pass bitreverse_tb.pass fftstage_tb.pass qtrstage_tb.pass
test: laststage_tb.pass butterfly_tb.pass hwbfly_tb.pass
test: fft_tb HEX # ifft_tb
	./fft_tb

mpy_tb.pass: mpy_tb
	./mpy_tb
	touch mpy_tb.pass

fftstage_tb.pass: fftstage_tb HEX
	./fftstage_tb
	touch fftstage_tb.pass

butterfly_tb.pass: butterfly_tb
	./butterfly_tb
	touch butterfly_tb.pass

hwbfly_tb.pass: hwbfly_tb
	./hwbfly_tb
	touch hwbfly_tb.pass

qtrstage_tb.pass: qtrstage_tb
	./qtrstage_tb
	touch qtrstage_tb.pass

laststage_tb.pass: laststage_tb
	./laststage_tb
	touch laststage_tb.pass

bitreverse_tb.pass: bitreverse_tb
	./bitreverse_tb
	touch bitreverse_tb.pass

.PHONY: clean
clean:
	rm -f mpy_tb bitreverse_tb fftstage_tb qtrstage_tb butterfly_tb
	rm -f fftstage_tb fft_tb ifft_tb hwbfly_tb laststage_tb
	rm -rf fft_tb.dbl ifft_tb.dbl
	rm -rf *cmem_*.hex
	rm -rf *.pass *.vcd

include $(VERILATOR_ROOT)/include/verilated.mk
