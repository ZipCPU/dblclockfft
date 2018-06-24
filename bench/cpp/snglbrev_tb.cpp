////////////////////////////////////////////////////////////////////////////
//
// Filename: 	snglbrev_tb.cpp
//
// Project:	A General Purpose Pipelined FFT Implementation
//
// Purpose:	A test-bench for the snglbrev.v subfile of the pipelined
//		FFT.  This file may be run autonomously.  If so, the last line
//	output will either read "SUCCESS" on success, or some other failure
//	message otherwise.
//
//	This file depends upon verilator to both compile, run, and therefore
//	test snglbrev.v
//
// Creator:	Dan Gisselquist, Ph.D.
//		Gisselquist Technology, LLC
//
///////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2015,2018, Gisselquist Technology, LLC
//
// This program is free software (firmware): you can redistribute it and/or
// modify it under the terms of  the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or (at
// your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
//
// You should have received a copy of the GNU General Public License along
// with this program.  (It's in the $(ROOT)/doc directory, run make with no
// target there if the PDF file isn't present.)  If not, see
// <http://www.gnu.org/licenses/> for a copy.
//
// License:	GPL, v3, as defined and found on www.gnu.org,
//		http://www.gnu.org/licenses/gpl.html
//
//
///////////////////////////////////////////////////////////////////////////
#include "Vsnglbrev.h"
#include "verilated.h"

#include "fftsize.h"

#define	FFTBITS	TST_DBLREVERSE_LGSIZE
#define	FFTSIZE	(1<<(FFTBITS))
#define	FFTMASK	(FFTSIZE-1)
#define	DATALEN	(1<<(FFTBITS+1))
#define	DATAMSK	(DATALEN-1)
#define	PAGEMSK	(FFTSIZE)

#ifdef	NEW_VERILATOR
#define	VVAR(A)	snglbrev__DOT_ ## A
#else
#define	VVAR(A)	v__DOT_ ## A
#endif

#define	iaddr		VVAR(_wraddr)
#define	in_reset	VVAR(_in_reset)

void	tick(Vsnglbrev *snglbrev) {
	snglbrev->i_clk = 0;
	snglbrev->eval();
	snglbrev->i_clk = 1;
	snglbrev->eval();

	snglbrev->i_ce = 0;
}

void	reset(Vsnglbrev *snglbrev) {
	snglbrev->i_ce  = 0;
	snglbrev->i_reset = 1;
	tick(snglbrev);
	snglbrev->i_ce  = 0;
	snglbrev->i_reset = 0;
	tick(snglbrev);
}

unsigned long	bitrev(const int nbits, const unsigned long vl) {
	unsigned long	r = 0;
	unsigned long	val = vl;

	for(int k=0; k<nbits; k++) {
		r <<= 1;
		r |= (val & 1);
		val >>= 1;
	}

	return r;
}

int	main(int argc, char **argv, char **envp) {
	Verilated::commandArgs(argc, argv);
	Vsnglbrev	*snglbrev = new Vsnglbrev;
	int syncd = 0;
	unsigned long	datastore[DATALEN], dataidx=0;
	const	int	BREV_OFFSET = 0;

	reset(snglbrev);

	printf("FFTSIZE = %08x\n", FFTSIZE);
	printf("FFTMASK = %08x\n", FFTMASK);
	printf("DATALEN = %08x\n", DATALEN);
	printf("DATAMSK = %08x\n", DATAMSK);

	for(int k=0; k<4*(FFTSIZE); k++) {
		snglbrev->i_ce = 1;
		snglbrev->i_in = k;
		datastore[(dataidx++)&(DATAMSK)] = snglbrev->i_in;
		tick(snglbrev);

		printf("k=%3d: IN = %6lx, OUT = %6lx, SYNC = %d\t(%2x) %d\n",
			k, snglbrev->i_in, snglbrev->o_out, snglbrev->o_sync,
			snglbrev->iaddr, snglbrev->in_reset);

		if ((k>BREV_OFFSET)&&((BREV_OFFSET==(k&FFTMASK))?1:0) != snglbrev->o_sync) {
			fprintf(stdout, "FAIL, BAD SYNC (k = %d > %d)\n", k, BREV_OFFSET);
			exit(EXIT_FAILURE);
		} else if (snglbrev->o_sync) {
			syncd = 1;
		}
		if ((syncd)&&((snglbrev->o_out&FFTMASK) != bitrev(FFTBITS, k-BREV_OFFSET))) {
			fprintf(stdout, "FAIL: BITREV.0 of k (%2x) = %2lx, not %2lx\n",
				k, snglbrev->o_out, bitrev(FFTBITS, (k-BREV_OFFSET)));
			exit(EXIT_FAILURE);
		}
	}

	for(int k=0; k<4*(FFTSIZE); k++) {
		snglbrev->i_ce = 1;
		snglbrev->i_in = rand() & 0x0ffffff;
		datastore[(dataidx++)&(DATAMSK)] = snglbrev->i_in;
		tick(snglbrev);

		printf("k=%3d: IN = %6lx, OUT = %6lx, SYNC = %d\n",
			k, snglbrev->i_in, snglbrev->o_out, snglbrev->o_sync);

		/*
		if ((k>BREV_OFFSET)&&(((BREV_OFFSET==(k&(FFTMASK>>1)))?1:0) != snglbrev->o_sync)) {
			fprintf(stdout, "FAIL, BAD SYNC\n");
			// exit(EXIT_FAILURE);
		} else
		*/
		if (snglbrev->o_sync)
			syncd = 1;
		if ((syncd)&&(snglbrev->o_out != datastore[
				(((dataidx-1-FFTSIZE)&PAGEMSK)
				+ bitrev(FFTBITS,
					(dataidx-FFTSIZE-1)&FFTMASK))])) {
			fprintf(stdout, "FAIL: BITREV.0 of k (%2x) = %2lx, not %2lx (expected %lx -> %lx)\n",
				k, snglbrev->o_out,
				datastore[(((dataidx-1-FFTSIZE)&PAGEMSK)
					+ bitrev(FFTBITS, (dataidx-FFTSIZE-1)&FFTMASK))],
				(dataidx-2)&DATAMSK,
				(((dataidx-2)&PAGEMSK)
					+ bitrev(FFTBITS, (dataidx-FFTSIZE-1)&FFTMASK)));
			exit(EXIT_FAILURE);
		}
	}

	delete	snglbrev;

	printf("SUCCESS!\n");
	exit(0);
}
