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
all: bimpy longbimpy hwbfly fftstage qtrstage laststage bitreverse # butterfly

bimpy:
	sby -f bimpy.sby

longbimpy:
	sby -f longbimpy.sby

fftstage:
	sby -f fftstage.sby

hwbfly:
	sby -f hwbfly.sby

butterfly:
	sby -f butterfly.sby

qtrstage:
	sby -f qtrstage.sby

laststage:
	sby -f laststage.sby

bitreverse:
	sby -f bitreverse.sby

windowfn:
	sby -f windowfn.sby

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
