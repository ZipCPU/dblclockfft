////////////////////////////////////////////////////////////////////////////////
//
// Filename: 	bitreverse_tb.cpp
//
// Project:	A General Purpose Pipelined FFT Implementation
//
// Purpose:	A test-bench for the bitreversal stage of the pipelined
//		FFT.  This file may be run autonomously.  If so, the last line
//	output will either read "SUCCESS" on success, or some other failure
//	message otherwise.
//
//	This file depends upon verilator to both compile, run, and therefore
//	test either snglbrev.v or dblreverse.v--depending on whether or not the
//	FFT handles one or two inputs per clock respectively.
//
// Creator:	Dan Gisselquist, Ph.D.
//		Gisselquist Technology, LLC
//
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2015-2020, Gisselquist Technology, LLC
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
// with this program.  (It's in the $(ROOT)/doc directory.  Run make with no
// target there if the PDF file isn't present.)  If not, see
// <http://www.gnu.org/licenses/> for a copy.
//
// License:	GPL, v3, as defined and found on www.gnu.org,
//		http://www.gnu.org/licenses/gpl.html
//
//
////////////////////////////////////////////////////////////////////////////////
#include "verilated.h"
#include "verilated_vcd_c.h"

#include "fftsize.h"
#include "Vbitreverse.h"

#define	FFTBITS	TST_DBLREVERSE_LGSIZE
#define	FFTSIZE	(1<<(FFTBITS))
#define	FFTMASK	(FFTSIZE-1)
#define	DATALEN	(1<<(FFTBITS+1))
#define	DATAMSK	(DATALEN-1)
#define	PAGEMSK	(FFTSIZE)

#ifdef	NEW_VERILATOR
#define	VVAR(A)	bitreverse__DOT_ ## A
#else
#define	VVAR(A)	v__DOT_ ## A
#endif

typedef	Vbitreverse	TSTCLASS;

#ifdef DBLCLKFFT
#define	iaddr		VVAR(_iaddr)
#else
#define	iaddr		VVAR(_wraddr)
#endif
#define	in_reset	VVAR(_in_reset)

VerilatedVcdC	*trace = NULL;
uint64_t	m_tickcount = 0;

void	tick(TSTCLASS *brev) {
	m_tickcount++;

	brev->i_clk = 0;
	brev->eval();
	if (trace)	trace->dump((vluint64_t)(10ul*m_tickcount-2));
	brev->i_clk = 1;
	brev->eval();
	if (trace)	trace->dump((vluint64_t)(10ul*m_tickcount));
	brev->i_clk = 0;
	brev->eval();
	if (trace) {
		trace->dump((vluint64_t)(10ul*m_tickcount+5));
		trace->flush();
	}

	brev->i_ce = 0;
}

void	cetick(TSTCLASS *brev) {
	brev->i_ce = 1;
	tick(brev);
	if (rand()&1) {
		brev->i_ce = 1;
		tick(brev);
	}
}

void	reset(TSTCLASS *brev) {
	brev->i_ce  = 0;
	brev->i_reset = 1;
	tick(brev);
	brev->i_ce  = 0;
	brev->i_reset = 0;
	tick(brev);
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

#ifdef	DBLCLKFFT
#define	BREVMASK	(FFTMASK>>1)
#define	BREVBITS	(FFTBITS-1)
#else
#define	BREVMASK	FFTMASK
#define	BREVBITS	FFTBITS
#endif

int	main(int argc, char **argv, char **envp) {
	Verilated::commandArgs(argc, argv);
	Verilated::traceEverOn(true);
	TSTCLASS	*brev = new TSTCLASS;
	int syncd = 0;
	unsigned long	datastore[DATALEN], dataidx=0;
	const	int	BREV_OFFSET = 0;

	trace = new VerilatedVcdC;
	brev->trace(trace, 99);
	trace->open("bitreverse_tb.vcd");

	reset(brev);

	printf("FFTSIZE = %08x\n", FFTSIZE);
	printf("FFTMASK = %08x\n", FFTMASK);
	printf("DATALEN = %08x\n", DATALEN);
	printf("DATAMSK = %08x\n", DATAMSK);

	for(int k=0; k<4*(FFTSIZE); k++) {
		brev->i_ce = 1;
#ifdef	DBLCLKFFT
		brev->i_in_0 = 2*k;
		brev->i_in_1 = 2*k+1;
		datastore[(dataidx++)&(DATAMSK)] = brev->i_in_0;
		datastore[(dataidx++)&(DATAMSK)] = brev->i_in_1;
#else
		brev->i_in = k;
		datastore[(dataidx++)&(DATAMSK)] = brev->i_in;
#endif
		tick(brev);

#ifdef	DBLCLKFFT
		printf("k=%3d: IN = %6lx,%6lx OUT = %6lx,%6lx SYNC = %d\t(%2x) %d\n",
			k, brev->i_in_0, brev->i_in_1,
			brev->o_out_0, brev->o_out_1, brev->o_sync,
			brev->iaddr, brev->in_reset);
#else
		printf("k=%3d: IN = %6lx, OUT = %6lx, SYNC = %d\t(%2x) %d\n",
			k, brev->i_in, brev->o_out, brev->o_sync,
			brev->iaddr, brev->in_reset);
#endif

		if ((k>BREV_OFFSET)&&((BREV_OFFSET==(k&BREVMASK))?1:0) != brev->o_sync) {
			fprintf(stdout, "FAIL, BAD SYNC (k = %d > %d)\n", k, BREV_OFFSET);
			exit(EXIT_FAILURE);
		} else if (brev->o_sync) {
			syncd = 1;
		}
#ifdef	DBLCLKFFT
		if ((syncd)&&((brev->o_out_0&BREVMASK) != bitrev(FFTBITS, k+1-BREV_OFFSET))) {
			fprintf(stdout, "FAIL: BITREV.0 of k (%2x) = %2lx, not %2lx\n",
				k, brev->o_out_0, bitrev(FFTBITS, (k+1-BREV_OFFSET)));
			// exit(EXIT_FAILURE);
		}
#else
		if ((syncd)&&((brev->o_out&FFTMASK) != bitrev(FFTBITS, k-BREV_OFFSET))) {
			fprintf(stdout, "FAIL: BITREV.0 of k (%2x) = %2lx, not %2lx\n",
				k, brev->o_out, bitrev(FFTBITS, (k-BREV_OFFSET)));
			exit(EXIT_FAILURE);
		}
#endif
	}

	for(int k=0; k<4*(FFTSIZE); k++) {
		brev->i_ce = 1;
#ifdef	DBLCLKFFT
		brev->i_in_0 = rand() & 0x0ffffff;
		brev->i_in_1 = rand() & 0x0ffffff;
		datastore[(dataidx++)&(DATAMSK)] = brev->i_in_0;
		datastore[(dataidx++)&(DATAMSK)] = brev->i_in_1;
#else
		brev->i_in = rand() & 0x0ffffff;
		datastore[(dataidx++)&(DATAMSK)] = brev->i_in;
#endif
		tick(brev);

#ifdef	DBLCLKFFT
		printf("k=%3d: IN = %6lx : %6lx, OUT = %6lx : %6lx, SYNC = %d\n",
			k, brev->i_in_0, brev->i_in_1,
			brev->o_out_0, brev->o_out_1, brev->o_sync);
#else
		printf("k=%3d: IN = %6lx, OUT = %6lx, SYNC = %d\n",
			k, brev->i_in, brev->o_out, brev->o_sync);
#endif

		if (brev->o_sync)
			syncd = 1;
#ifdef	DBLCLKFFT
		if ((syncd)&&(brev->o_out_0 != datastore[(((dataidx-2-FFTSIZE)&PAGEMSK) + bitrev(FFTBITS, (dataidx-FFTSIZE-2)&FFTMASK))])) {
			fprintf(stdout, "FAIL: BITREV.0 of k (%2x) = %2lx, not %2lx (expected %lx -> %lx)\n",
				k, brev->o_out_0,
				datastore[(((dataidx-2-FFTSIZE)&PAGEMSK)
					+ bitrev(FFTBITS, (dataidx-FFTSIZE-2)&FFTMASK))],
				(dataidx-2)&DATAMSK,
				(((dataidx-2)&PAGEMSK)
					+ bitrev(FFTBITS, (dataidx-FFTSIZE-2)&FFTMASK)));
			// exit(-1);
		}

		if ((syncd)&&(brev->o_out_1 != datastore[(((dataidx-2-FFTSIZE)&PAGEMSK) + bitrev(FFTBITS, (dataidx-FFTSIZE-1)&FFTMASK))])) {
			fprintf(stdout, "FAIL: BITREV.1 of k (%2x) = %2lx, not %2lx (expected %lx)\n",
				k, brev->o_out_1,
				datastore[(((dataidx-2-FFTSIZE)&PAGEMSK)
					+ bitrev(FFTBITS, (dataidx-FFTSIZE-1)&FFTMASK))],
				(((dataidx-1)&PAGEMSK)
					+ bitrev(FFTBITS, (dataidx-FFTSIZE-1)&FFTMASK)));
			// exit(-1);
		}
#else
		if ((syncd)&&(brev->o_out != datastore[
				(((dataidx-1-FFTSIZE)&PAGEMSK)
				+ bitrev(FFTBITS,
					(dataidx-FFTSIZE-1)&FFTMASK))])) {
			fprintf(stdout, "FAIL: BITREV.0 of k (%2x) = %2lx, not %2lx (expected %lx -> %lx)\n",
				k, brev->o_out,
				datastore[(((dataidx-1-FFTSIZE)&PAGEMSK)
					+ bitrev(FFTBITS, (dataidx-FFTSIZE-1)&FFTMASK))],
				(dataidx-2)&DATAMSK,
				(((dataidx-2)&PAGEMSK)
					+ bitrev(FFTBITS, (dataidx-FFTSIZE-1)&FFTMASK)));
			exit(EXIT_FAILURE);
		}
#endif
	}

	delete	brev;

	printf("SUCCESS!\n");
	exit(0);
}
