################################################################################
##
## Filename: 	doc/Makefile
##
## Project:	A General Purpose Pipelined FFT Implementation
##
## Purpose:	This simple Makefile builds the documentation for the Pipelined
##		FFT Generator.
##
## Targets:
##
##	The "all" target builds both a copy of the LGPL license as well as the
##	specification (and directions for) the FFT core generator.
##
## Creator:	Dan Gisselquist, Ph.D.
##		Gisselquist Technology, LLC
##
################################################################################
##
## Copyright (C) 2018-2020, Gisselquist Technology, LLC
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
.PHONY: all license spec pdf gpl lgpl
all:	lgpl spec
license: lgpl
gpl: lgpl
lgpl: lgpl-3.0.pdf
spec:	spec.pdf
pdf: lgpl spec

DSRC := src

.DELETE_ON_ERROR:

LICENSE := lgpl-3.0
$(LICENSE).pdf: $(DSRC)/$(LICENSE).tex
	latex $(DSRC)/$(LICENSE).tex
	latex $(DSRC)/$(LICENSE).tex
	dvips -q -z -t letter -P pdf -o $(LICENSE).ps $(LICENSE).dvi
	ps2pdf -dAutoRotatePages=/All $(LICENSE).ps $(LICENSE).pdf
	rm -f $(LICENSE).dvi $(LICENSE).log $(LICENSE).aux $(LICENSE).ps

spec.pdf: $(DSRC)/spec.tex $(DSRC)/gqtekspec.cls
	cd $(DSRC)/; latex spec.tex
	cd $(DSRC)/; latex spec.tex
	dvips -q -z -t letter -P pdf -o spec.ps $(DSRC)/spec.dvi
	ps2pdf -dAutoRotatePages=/All spec.ps spec.pdf
	rm -f $(DSRC)/spec.dvi $(DSRC)/spec.log
	rm -f $(DSRC)/spec.aux $(DSRC)/spec.toc
	rm -f $(DSRC)/spec.lof $(DSRC)/spec.lot
	rm -f spec.ps

