############################################################################/
##
## Filename: 	Makefile
##
## Project:	A Generic Pipelined FFT Implementation
##
## Purpose:	This is the main Makefile for the FFT core generator.
##		It is very simple in its construction, the most complicated
##	parts being the building of the Verilator simulation--a step that may
##	not be required for your project.
##
##	To build the FFT generator, just type 'make' on a line by itself.  For
##	a quick tutorial in how to run the generator, just type './fftgen -h'
##	to read the usage() statement.
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
VFLAGS  := -Wall -O3 -MMD --trace -cc
CFLAGS  := -g -Wall

fftgen: fftgen.o
	$(CXX) $(CFLAGS) $< -o $@

%.o: %.cpp
	$(CXX) $(CFLAGS) -c $< -o $@

.PHONY: test
test: fft ifft butterfly dblreverse qtrstage dblstage fftstage
test: hwbfly shiftaddmpy longbimpy

#
# Although these parameters, a 2048 point FFT of 16 bits input, aren't
# the only parameters the FFT can do, they are the ones that the test
# benches depend upon.  If you change these, and you are welcome to do so,
# you may need to adjust the test benches if you wish to prove that your
# changes work.
#
.PHONY: fft
fft: $(OBJDR)/Vfftmain__ALL.so
$(CORED)/fftmain.v: fftgen
	./fftgen -v -1 -f $(TESTSZ) -n 16 -p 6 -a $(BENCHD)/fftsize.h
$(OBJDR)/Vfftmain.h: $(CORED)/fftmain.v
	cd $(CORED)/; $(VERILATOR) $(VFLAGS) fftmain.v
$(OBJDR)/Vfftmain__ALL.so: $(OBJDR)/Vfftmain.h
	cd $(OBJDR); make -f Vfftmain.mk

.PHONY: dblfft
dblfft: $(OBJDR)/Vdblfftmain__ALL.so
$(CORED)/dblfftmain.v: fftgen
	./fftgen -v -f $(TESTSZ) -n 16 -p 6 -a $(BENCHD)/dblfftsize.h
$(OBJDR)/Vdblfftmain.h: $(CORED)/dblfftmain.v
	cd $(CORED)/; $(VERILATOR) $(VFLAGS) dblfftmain.v
$(OBJDR)/Vdblfftmain__ALL.so: $(OBJDR)/Vdblfftmain.h
	cd $(OBJDR); make -f Vdblfftmain.mk

.PHONY: idblfft
idblfft: $(OBJDR)/Vidblfftmain__ALL.so
$(CORED)/idblfftmain.v: fftgen
	./fftgen -f $(TESTSZ) -i -n 22  -p 6 -a $(BENCHD)/idblfftsize.h
$(OBJDR)/Vidblfftmain.h: $(CORED)/idblfftmain.v
	cd $(CORED)/; $(VERILATOR) $(VFLAGS) idblfftmain.v
$(OBJDR)/Vidblfftmain__ALL.so: $(OBJDR)/Vidblfftmain.h
	cd $(OBJDR); make -f Vidblfftmain.mk

.PHONY: ifft
ifft: $(OBJDR)/Vifftmain__ALL.so
$(CORED)/ifftmain.v: fftgen
	./fftgen -1 -f $(TESTSZ) -i -n 22  -p 6 -a $(BENCHD)/ifftsize.h
$(OBJDR)/Vifftmain.h: $(CORED)/ifftmain.v
	cd $(CORED)/; $(VERILATOR) $(VFLAGS) ifftmain.v
$(OBJDR)/Vifftmain__ALL.so: $(OBJDR)/Vifftmain.h
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

.PHONY: fftstage
dblstage: $(OBJDR)/Vfftstage__ALL.a

$(CORED)/fftstage.v: fft
$(OBJDR)/Vfftstage.cpp $(OBJDR)/Vfftstage.h: $(CORED)/fftstage.v
	cd $(CORED)/; $(VERILATOR) $(VFLAGS) fftstage.v
$(OBJDR)/Vfftstage__ALL.a: $(OBJDR)/Vfftstage.h
$(OBJDR)/Vfftstage__ALL.a: $(OBJDR)/Vfftstage.cpp
	cd $(OBJDR)/; make -f Vfftstage.mk

.PHONY: clean
clean:
	rm fftgen fftgen.o
	rm -rf $(CORED)


