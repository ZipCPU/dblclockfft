############################################################################/
##
## Filename: 	Makefile
##
## Project:	A Doubletime Pipelined FFT
##
## Purpose:	This is the main Makefile for the FFT core generator.
##		It is very simple in its construction, the most complicated
##		parts being the building of the Verilator simulation--a
##		step that may not be required for your project.
##
##		To build the FFT generator, just type 'make' on a line
##		by itself.  For a quick tutorial in how to run the
##		generator, just type './fftgen -h' to read the usage()
##		statement.
##
## Creator:	Dan Gisselquist, Ph.D.
##		Gisselquist Technology, LLC
##
##########################################################################/
##
## Copyright (C) 2015, Gisselquist Technology, LLC
##
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
## with this program.  (It's in the $(ROOT)/doc directory, run make with no
## target there if the PDF file isn't present.)  If not, see
## <http:##www.gnu.org/licenses/> for a copy.
##
## License:	GPL, v3, as defined and found on www.gnu.org,
##		http:##www.gnu.org/licenses/gpl.html
##
##
##########################################################################/
##
##
# This is really simple ...
all: fftgen
CORED := fft-core
OBJDR := $(CORED)/obj_dir
TESTSZ := 2048
BENCHD := ../bench/cpp
ifneq ($(VERILATOR_ROOT),)
VERILATOR:=$(VERILATOR_ROOT)/bin/verilator
else
VERILATOR:=verilator
VERILATOR_ROOT ?= $(shell bash -c 'verilator -V|grep VERILATOR_ROOT | head -1 | sed -e " s/^.*=\s*//"')
endif
export	$(VERILATOR)
VROOT   := $(VERILATOR_ROOT)
VFLAGS  := -Wall -MMD --trace -cc

fftgen: fftgen.o
	$(CXX) $< -o $@

%.o: %.cpp
	$(CXX) -c $< -o $@

.PHONY: test
test: fft ifft butterfly dblreverse qtrstage dblstage fftstage_o2048
test: hwbfly shiftaddmpy longbimpy

#
# Although these parameters, a 2048 point FFT of 16 bits input, aren't
# the only parameters the FFT can do, they are the ones that the test
# benches depend upon.  If you change these, and you are welcome to do so,
# you may need to adjust the test benches if you wish to prove that your
# changes work.
#
.PHONY: fft
fft: fftgen
	./fftgen -f $(TESTSZ) -n 16 -p 6 -a $(BENCHD)/fftsize.h
	cd $(CORED)/; $(VERILATOR) $(VFLAGS) fftmain.v
	cd $(OBJDR); make -f Vfftmain.mk

.PHONY: ifft
ifft: fftgen
	./fftgen -f $(TESTSZ) -i -n 22  -p 6 -a $(BENCHD)/ifftsize.h
	cd $(CORED)/; $(VERILATOR) $(VFLAGS) ifftmain.v
	cd $(OBJDR); make -f Vifftmain.mk

.PHONY: shiftaddmpy
shiftaddmpy: $(OBJDR)/Vshiftaddmpy__ALL.a

$(CORED)/shiftaddmpy.v: fft
$(OBJDR)/Vshiftaddmpy.cpp $(OBJDR)/Vshiftaddmpy.h: $(CORED)/shiftaddmpy.v
	cd $(CORED)/; $(VERILATOR) $(VFLAGS) shiftaddmpy.v
$(OBJDR)/Vshiftaddmpy__ALL.a: $(OBJDR)/Vshiftaddmpy.h
$(OBJDR)/Vshiftaddmpy__ALL.a: $(OBJDR)/Vshiftaddmpy.cpp
	cd $(OBJDR)/; make -f Vshiftaddmpy.mk

.PHONY: longbimpy
longbimpy: $(OBJDR)/Vlongbimpy__ALL.a

$(CORED)/longbimpy.v: fft
$(OBJDR)/Vlongbimpy.cpp $(OBJDR)/Vlongbimpy.h: $(CORED)/longbimpy.v
	cd $(CORED)/; $(VERILATOR) $(VFLAGS) longbimpy.v
$(OBJDR)/Vlongbimpy__ALL.a: $(OBJDR)/Vlongbimpy.h
$(OBJDR)/Vlongbimpy__ALL.a: $(OBJDR)/Vlongbimpy.cpp
	cd $(OBJDR)/; make -f Vlongbimpy.mk

.PHONY: butterfly
butterfly: $(OBJDR)/Vbutterfly__ALL.a

$(CORED)/butterfly.v: fft
$(OBJDR)/Vbutterfly.cpp $(OBJDR)/Vbutterfly.h: $(CORED)/butterfly.v
	cd $(CORED)/; $(VERILATOR) $(VFLAGS) butterfly.v
$(OBJDR)/Vbutterfly__ALL.a: $(OBJDR)/Vbutterfly.h
$(OBJDR)/Vbutterfly__ALL.a: $(OBJDR)/Vbutterfly.cpp
	cd $(OBJDR)/; make -f Vbutterfly.mk

.PHONY: hwbfly
hwbfly: $(OBJDR)/Vhwbfly__ALL.a

$(CORED)/hwbfly.v: fft
$(OBJDR)/Vhwbfly.cpp $(OBJDR)/Vhwbfly.h: $(CORED)/hwbfly.v
	cd $(CORED)/; $(VERILATOR) $(VFLAGS) hwbfly.v
$(OBJDR)/Vhwbfly__ALL.a: $(OBJDR)/Vhwbfly.h
$(OBJDR)/Vhwbfly__ALL.a: $(OBJDR)/Vhwbfly.cpp
	cd $(OBJDR)/; make -f Vhwbfly.mk

.PHONY: dblreverse
dblreverse: $(OBJDR)/Vdblreverse__ALL.a

$(CORED)/dblreverse.v: fft
$(OBJDR)/Vdblreverse.cpp $(OBJDR)/Vdblreverse.h: $(CORED)/dblreverse.v
	cd $(CORED)/; $(VERILATOR) $(VFLAGS) dblreverse.v
$(OBJDR)/Vdblreverse__ALL.a: $(OBJDR)/Vdblreverse.h
$(OBJDR)/Vdblreverse__ALL.a: $(OBJDR)/Vdblreverse.cpp
	cd $(OBJDR)/; make -f Vdblreverse.mk

.PHONY: qtrstage
qtrstage: $(OBJDR)/Vqtrstage__ALL.a

$(CORED)/qtrstage.v: fft
$(OBJDR)/Vqtrstage.cpp $(OBJDR)/Vqtrstage.h: $(CORED)/qtrstage.v
	cd $(CORED)/; $(VERILATOR) $(VFLAGS) qtrstage.v
$(OBJDR)/Vqtrstage__ALL.a: $(OBJDR)/Vqtrstage.h
$(OBJDR)/Vqtrstage__ALL.a: $(OBJDR)/Vqtrstage.cpp
	cd $(OBJDR)/; make -f Vqtrstage.mk

.PHONY: dblstage
dblstage: $(OBJDR)/Vdblstage__ALL.a

$(CORED)/dblstage.v: fft
$(OBJDR)/Vdblstage.cpp $(OBJDR)/Vdblstage.h: $(CORED)/dblstage.v
	cd $(CORED)/; $(VERILATOR) $(VFLAGS) dblstage.v
$(OBJDR)/Vdblstage__ALL.a: $(OBJDR)/Vdblstage.h
$(OBJDR)/Vdblstage__ALL.a: $(OBJDR)/Vdblstage.cpp
	cd $(OBJDR)/; make -f Vdblstage.mk

.PHONY: fftstage_o2048
dblstage: $(OBJDR)/Vfftstage_o2048__ALL.a

$(CORED)/fftstage_o2048.v: fft
$(OBJDR)/Vfftstage_o2048.cpp $(OBJDR)/Vfftstage_o2048.h: $(CORED)/fftstage_o2048.v
	cd $(CORED)/; $(VERILATOR) $(VFLAGS) fftstage_o2048.v
$(OBJDR)/Vfftstage_o2048__ALL.a: $(OBJDR)/Vfftstage_o2048.h
$(OBJDR)/Vfftstage_o2048__ALL.a: $(OBJDR)/Vfftstage_o2048.cpp
	cd $(OBJDR)/; make -f Vfftstage_o2048.mk

.PHONY: clean
clean:
	rm fftgen fftgen.o
	rm -rf $(CORED)


