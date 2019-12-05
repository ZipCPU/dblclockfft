################################################################################
##
## Filename: 	bench/formal/Makefile
##
## Project:	A General Purpose Pipelined FFT Implementation
##
## Purpose:	Runs all of the SymbiYosys enabled formal proofs.  Requires
##		the formal code within the various modules to be enabled.
##
## Creator:	Dan Gisselquist, Ph.D.
##		Gisselquist Technology, LLC
##
################################################################################
##
## Copyright (C) 2018, Gisselquist Technology, LLC
##
## This file is part of the general purpose pipelined FFT project.
##
## The pipelined FFT project is free software (firmware): you can redistribute
## it and/or modify it under the terms of the GNU Lesser General Public License
## as published by the Free Software Foundation, either version 3 of the
## License, or (at your option) any later version.
##
## The pipelined FFT project is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTIBILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
## General Public License for more details.
##
## You should have received a copy of the GNU Lesser General Public License
## along with this program.  (It's in the $(ROOT)/doc directory.  Run make
## with no target there if the PDF file isn't present.)  If not, see
## <http://www.gnu.org/licenses/> for a copy.
##
## License:	LGPL, v3, as defined and found on www.gnu.org,
##		http://www.gnu.org/licenses/lgpl.html
##
##
################################################################################
##
##
TARGETS := bimpy longbimpy fftstage hwbfly butterfly qtrstage laststage bitreverse
.PHONY: $(TARGETS)
all: bimpy longbimpy hwbfly fftstage qtrstage laststage bitreverse butterfly
RTL := ../../rtl

.PHONY: bimpy
bimpy: bimpy/PASS
bimpy/PASS: bimpy.sby ../../rtl/bimpy.v
	sby -f bimpy.sby

longbimpy: longbimpy/PASS
longbimpy/PASS: longbimpy.sby ../../rtl/longbimpy.v
	sby -f longbimpy.sby

fftstage: fftstage_s2/PASS fftstage_s3/PASS fftstage_s4/PASS fftstage_s5/PASS fftstage_s6/PASS
fftstage_s2/PASS: fftstage.sby $(RTL)/fftstage.v
	sby -f fftstage.sby s2
fftstage_s3/PASS: fftstage.sby $(RTL)/fftstage.v
	sby -f fftstage.sby s3
fftstage_s4/PASS: fftstage.sby $(RTL)/fftstage.v
	sby -f fftstage.sby s4
fftstage_s5/PASS: fftstage.sby $(RTL)/fftstage.v
	sby -f fftstage.sby s5
fftstage_s6/PASS: fftstage.sby $(RTL)/fftstage.v
	sby -f fftstage.sby s6

hwbfly: hwbfly_one/PASS hwbfly_two/PASS hwbfly_three/PASS
hwbfly_one/PASS: hwbfly.sby $(RTL)/hwbfly.v
	sby -f hwbfly.sby one
hwbfly_two/PASS: hwbfly.sby $(RTL)/hwbfly.v
	sby -f hwbfly.sby two
hwbfly_three/PASS: hwbfly.sby $(RTL)/hwbfly.v
	sby -f hwbfly.sby three

.PHONY: butterfly_ck1 butterfly_ck2 butterfly_ck3
butterfly: butterfly_ck1 butterfly_ck2 butterfly_ck3

butterfly_ck1: butterfly_ck1/PASS
butterfly_ck2: butterfly_ck2_r0/PASS butterfly_ck2_r1/PASS
butterfly_ck3: butterfly_ck3_r0/PASS butterfly_ck3_r1/PASS butterfly_ck3_r2/PASS

SOFTMPY := $(RTL)/longbimpy.v $(RTL)/bimpy.v $(RTL)/convround.v

butterfly_ck1/PASS: butterfly.sby    $(SOFTMPY)
	sby -f butterfly.sby ck1
butterfly_ck2_r0/PASS: butterfly.sby $(SOFTMPY)
	sby -f butterfly.sby ck2_r0
butterfly_ck2_r1/PASS: butterfly.sby $(SOFTMPY)
	sby -f butterfly.sby ck2_r1
butterfly_ck3_r0/PASS: butterfly.sby $(SOFTMPY)
	sby -f butterfly.sby ck3_r0
butterfly_ck3_r1/PASS: butterfly.sby $(SOFTMPY)
	sby -f butterfly.sby ck3_r1
butterfly_ck3_r2/PASS: butterfly.sby $(SOFTMPY)
	sby -f butterfly.sby ck3_r2

qtrstage: qtrstage/PASS
qtrstage/PASS: qtrstage.sby $(RTL)/qtrstage.v
	sby -f qtrstage.sby

laststage: laststage/PASS
laststage/PASS: laststage.sby $(RTL)/laststage.v
	sby -f laststage.sby

bitreverse: bitreverse/PASS
bitreverse/PASS: bitreverse.sby $(RTL)/bitreverse.v
	sby -f bitreverse.sby

windowfn: windowfn_cvr/PASS windowfn_proof/PASS
windowfn_cvr/PASS: windowfn.sby $(RTL)/windowfn.v
	sby -f windowfn.sby cvr
windowfn_proof/PASS: windowfn.sby $(RTL)/windowfn.v
	sby -f windowfn.sby proof

clean:
	rm -rf bimpy/ longbimpy/
	rm -rf bitreverse/
	rm -rf butterfly_ck1/
	rm -rf butterfly_ck2_r0/
	rm -rf butterfly_ck2_r1/
	rm -rf butterfly_ck3_r0/
	rm -rf butterfly_ck3_r1/
	rm -rf butterfly_ck3_r2/
	rm -rf fftstage_s2/ fftstage_s3/ fftstage_s4/ fftstage_s5/ fftstage_s6/
	rm -rf hwbfly_one/ hwbfly_two/ hwbfly_three/
	rm -rf laststage/ qtrstage/
	rm -rf windowfn_cover/ windowfn_proof/
