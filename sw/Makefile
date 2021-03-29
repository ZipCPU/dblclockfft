################################################################################
##
## Filename: 	Makefile
## {{{
## Project:	A General Purpose Pipelined FFT Implementation
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
all: fftgen
CORED := ../rtl
VOBJDR  := $(CORED)/obj_dir
OBJDIR   := obj-pc
BENCHD  := ../bench/cpp
SOURCES := bitreverse.cpp bldstage.cpp butterfly.cpp fftgen.cpp fftlib.cpp \
		legal.cpp rounding.cpp softmpy.cpp
TESTSZ  := -f 2048
# CKPCE   := -1 -k 1
CKPCE   := -2
MPYS    := -p 0
IWID    := -n 15
FFTPARAMS := -d $(CORED) $(TESTSZ) $(CKPCE) $(MPYS) $(IWID)
OBJECTS := $(addprefix $(OBJDIR)/,$(subst .cpp,.o,$(SOURCES)))
HEADERS := $(wildcard *.h)
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

$(OBJDIR)/%.o: %.cpp
	$(mk-objdir)
	$(CXX) -c $(CFLAGS) $< -o $@

fftgen: $(OBJECTS)
	$(CXX) $(CFLAGS) $^ -o $@

.PHONY: test
test: fft ifft butterfly fftstage hwbfly longbimpy qtrstage
test: bitreverse laststage

#
# Although these parameters, a 2048 point FFT of 16 bits input, aren't
# the only parameters the FFT can do, they are the ones that the test
# benches depend upon.  If you change these, and you are welcome to do so,
# you may need to adjust the test benches if you wish to prove that your
# changes work.
#
.PHONY: fft forcedfft
## {{{
fft: $(VOBJDR)/Vfftmain__ALL.so
$(CORED)/fftmain.v: fftgen
	./fftgen -v $(FFTPARAMS) -a $(BENCHD)/fftsize.h
forcedfft: fftgen
	./fftgen -v $(FFTPARAMS) -a $(BENCHD)/fftsize.h
$(VOBJDR)/Vfftmain.h: $(CORED)/fftmain.v
	cd $(CORED)/; $(VERILATOR) $(VFLAGS) fftmain.v
$(VOBJDR)/Vfftmain__ALL.so: $(VOBJDR)/Vfftmain.h
	cd $(VOBJDR); make -f Vfftmain.mk
## }}}

.PHONY: dblfft
## {{{
dblfft: $(VOBJDR)/Vdblfftmain__ALL.so
$(CORED)/dblfftmain.v: fftgen
	./fftgen -v $(FFTPARAMS) -a $(BENCHD)/fftsize.h
$(VOBJDR)/Vdblfftmain.h: $(CORED)/dblfftmain.v
	cd $(CORED)/; $(VERILATOR) $(VFLAGS) dblfftmain.v
$(VOBJDR)/Vdblfftmain__ALL.so: $(VOBJDR)/Vdblfftmain.h
	cd $(VOBJDR); make -f Vdblfftmain.mk
## }}}

.PHONY: idblfft
## {{{
idblfft: $(VOBJDR)/Vidblfftmain__ALL.so
$(CORED)/idblfftmain.v: fftgen
	./fftgen -i $(FFTPARAMS) -a $(BENCHD)/ifftsize.h
$(VOBJDR)/Vidblfftmain.h: $(CORED)/idblfftmain.v
	cd $(CORED)/; $(VERILATOR) $(VFLAGS) idblfftmain.v
$(VOBJDR)/Vidblfftmain__ALL.so: $(VOBJDR)/Vidblfftmain.h
	cd $(VOBJDR); make -f Vidblfftmain.mk
## }}}

.PHONY: ifft
## {{{
ifft: $(VOBJDR)/Vifftmain__ALL.so
$(CORED)/ifftmain.v: fftgen
	./fftgen -i $(FFTPARAMS) -a $(BENCHD)/ifftsize.h
$(VOBJDR)/Vifftmain.h: $(CORED)/ifftmain.v
	cd $(CORED)/; $(VERILATOR) $(VFLAGS) ifftmain.v
$(VOBJDR)/Vifftmain__ALL.so: $(VOBJDR)/Vifftmain.h
	cd $(VOBJDR); make -f Vifftmain.mk
## }}}

.PHONY: longbimpy
## {{{
longbimpy: $(VOBJDR)/Vlongbimpy__ALL.a

$(CORED)/longbimpy.v: fft
$(VOBJDR)/Vlongbimpy.cpp $(VOBJDR)/Vlongbimpy.h: $(CORED)/longbimpy.v
	cd $(CORED)/; $(VERILATOR) $(VFLAGS) longbimpy.v
$(VOBJDR)/Vlongbimpy__ALL.a: $(VOBJDR)/Vlongbimpy.h
$(VOBJDR)/Vlongbimpy__ALL.a: $(VOBJDR)/Vlongbimpy.cpp
	cd $(VOBJDR)/; make -f Vlongbimpy.mk
## }}}

.PHONY: butterfly
## {{{
butterfly: $(VOBJDR)/Vbutterfly__ALL.a

$(CORED)/butterfly.v: fft
$(VOBJDR)/Vbutterfly.cpp $(VOBJDR)/Vbutterfly.h: $(CORED)/butterfly.v
	cd $(CORED)/; $(VERILATOR) $(VFLAGS) butterfly.v
$(VOBJDR)/Vbutterfly__ALL.a: $(VOBJDR)/Vbutterfly.h
$(VOBJDR)/Vbutterfly__ALL.a: $(VOBJDR)/Vbutterfly.cpp
	cd $(VOBJDR)/; make -f Vbutterfly.mk
## }}}

.PHONY: hwbfly
## {{{
hwbfly: $(VOBJDR)/Vhwbfly__ALL.a

$(CORED)/hwbfly.v: fft
$(VOBJDR)/Vhwbfly.cpp $(VOBJDR)/Vhwbfly.h: $(CORED)/hwbfly.v
	cd $(CORED)/; $(VERILATOR) $(VFLAGS) hwbfly.v
$(VOBJDR)/Vhwbfly__ALL.a: $(VOBJDR)/Vhwbfly.h
$(VOBJDR)/Vhwbfly__ALL.a: $(VOBJDR)/Vhwbfly.cpp
	cd $(VOBJDR)/; make -f Vhwbfly.mk
## }}}

.PHONY: bitreverse
## {{{
bitreverse: $(VOBJDR)/Vbitreverse__ALL.a

$(CORED)/bitreverse.v: fft
$(VOBJDR)/Vbitreverse.cpp $(VOBJDR)/Vbitreverse.h: $(CORED)/bitreverse.v
	cd $(CORED)/; $(VERILATOR) $(VFLAGS) bitreverse.v
$(VOBJDR)/Vbitreverse__ALL.a: $(VOBJDR)/Vbitreverse.h
$(VOBJDR)/Vbitreverse__ALL.a: $(VOBJDR)/Vbitreverse.cpp
	cd $(VOBJDR)/; make -f Vbitreverse.mk
## }}}

.PHONY: qtrstage
## {{{
qtrstage: $(VOBJDR)/Vqtrstage__ALL.a

$(CORED)/qtrstage.v: fft
$(VOBJDR)/Vqtrstage.cpp $(VOBJDR)/Vqtrstage.h: $(CORED)/qtrstage.v
	cd $(CORED)/; $(VERILATOR) $(VFLAGS) qtrstage.v
$(VOBJDR)/Vqtrstage__ALL.a: $(VOBJDR)/Vqtrstage.h
$(VOBJDR)/Vqtrstage__ALL.a: $(VOBJDR)/Vqtrstage.cpp
	cd $(VOBJDR)/; make -f Vqtrstage.mk
## }}}

.PHONY: dblstage
## {{{
dblstage: $(VOBJDR)/Vdblstage__ALL.a

$(CORED)/dblstage.v: fft
$(VOBJDR)/Vdblstage.cpp $(VOBJDR)/Vdblstage.h: $(CORED)/dblstage.v
	cd $(CORED)/; $(VERILATOR) $(VFLAGS) dblstage.v
$(VOBJDR)/Vdblstage__ALL.a: $(VOBJDR)/Vdblstage.h
$(VOBJDR)/Vdblstage__ALL.a: $(VOBJDR)/Vdblstage.cpp
	cd $(VOBJDR)/; make -f Vdblstage.mk
## }}}

.PHONY: laststage
## {{{
laststage: $(VOBJDR)/Vlaststage__ALL.a

$(CORED)/laststage.v: fft
$(VOBJDR)/Vlaststage.cpp $(VOBJDR)/Vlaststage.h: $(CORED)/laststage.v
	cd $(CORED)/; $(VERILATOR) $(VFLAGS) laststage.v
$(VOBJDR)/Vlaststage__ALL.a: $(VOBJDR)/Vlaststage.h
$(VOBJDR)/Vlaststage__ALL.a: $(VOBJDR)/Vlaststage.cpp
	cd $(VOBJDR)/; make -f Vlaststage.mk
## }}}

.PHONY: fftstage
## {{{
fftstage: $(VOBJDR)/Vfftstage__ALL.a

$(CORED)/fftstage.v: fft
$(VOBJDR)/Vfftstage.cpp $(VOBJDR)/Vfftstage.h: $(CORED)/fftstage.v
	cd $(CORED)/; $(VERILATOR) $(VFLAGS) fftstage.v
$(VOBJDR)/Vfftstage__ALL.a: $(VOBJDR)/Vfftstage.h
$(VOBJDR)/Vfftstage__ALL.a: $(VOBJDR)/Vfftstage.cpp
	cd $(VOBJDR)/; make -f Vfftstage.mk
## }}}


.PHONY: clean
## {{{
clean:
	rm -rf fft-core/
	rm -rf fftgen $(OBJDIR)/
	rm -rf $(CORED)/obj_dir
	rm -rf $(CORED)/fftmain.v $(CORED)/fftstage.v
	rm -rf $(CORED)/ifftmain.v $(CORED)/ifftstage.v
	rm -rf $(CORED)/qtrstage.v $(CORED)/laststage.v $(CORED)/bitreverse.v
	rm -rf $(CORED)/butterfly.v $(CORED)/hwbfly.v
	rm -rf $(CORED)/longbimpy.v $(CORED)/bimpy.v
	rm -rf $(CORED)/convround.v
	rm -rf $(CORED)/*cmem_*.hex
## }}}

## Automatic dependency handling
## {{{
# The "depends" target, to know what files things depend upon.  The depends
# file itself is kept in $(OBJDIR)/depends.txt
#
define	build-depends
	$(mk-objdir)
	@echo "Building dependency file"
	@$(CXX) $(CFLAGS) $(INCS) -MM $(SOURCES) > $(OBJDIR)/xdepends.txt
@sed -e 's/^.*.o: /$(OBJDIR)\/&/' < $(OBJDIR)/xdepends.txt > $(OBJDIR)/depends.txt
	@rm $(OBJDIR)/xdepends.txt
endef

.PHONY: depends
depends: tags
	$(build-depends)

$(OBJDIR)/depends.txt: depends

#
# Make a directory to hold all of the FFT-gen (i.e. the C++) build products
# (object files)
#
define	mk-objdir
	@bash -c "if [ ! -e $(OBJDIR) ]; then mkdir -p $(OBJDIR); fi"
endef

#
# The "tags" target
#
tags:	$(SOURCES) $(HEADERS)
	@echo "Generating tags"
	@ctags $(SOURCES) $(HEADERS)

ifneq ($(MAKECMDGOALS),clean)
-include $(OBJDIR)/depends.txt
endif
## }}}
