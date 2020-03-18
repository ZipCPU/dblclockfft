////////////////////////////////////////////////////////////////////////////////
//
// Filename: 	fftstage_tb.cpp
//
// Project:	A General Purpose Pipelined FFT Implementation
//
// Purpose:	A test-bench for a generic FFT stage which has been
//		instantiated by fftgen.  Without loss of (much) generality,
//	we'll examine the 2048 fftstage.v.  This file may be run autonomously.
//	If so, the last line output will either read "SUCCESS" on success, or
//	some other failure message otherwise.  Likewise the exit code will
//	also indicate success (exit(0)) or failure (anything else).
//
//	This file depends upon verilator to both compile, run, and therefore
//	test fftstage.v.  Also, you'll need to place a copy of the cmem_*2048
//	hex file into the directory where you run this test bench.
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
//
//
#include "Vfftstage.h"
#include "verilated.h"
#include "verilated_vcd_c.h"
#include "twoc.h"
#include "fftsize.h"


#ifdef	NEW_VERILATOR
#define	VVAR(A)	fftstage__DOT_ ## A
#else
#define	VVAR(A)	v__DOT_ ## A
#endif

#define	cmem	VVAR(_cmem)
#define	iaddr	VVAR(_iaddr)
#define	ib_sync	VVAR(_ib_sync)
#define	ib_a	VVAR(_ib_a)
#define	ib_b	VVAR(_ib_b)
#define	ib_c	VVAR(_ib_c)

#define	ob_sync	VVAR(_ob_sync)
#define	ob_a	VVAR(_ob_a)
#define	ob_b	VVAR(_ob_b)

#define	FFTBITS	(FFT_LGWIDTH)
#define	FFTLEN	(1<<FFTBITS)
#define	FFTSIZE	FFTLEN
#define	FFTMASK	(FFTLEN-1)
#define	IWIDTH	FFT_IWIDTH
#define	CWIDTH	20
#define	OWIDTH	(FFT_IWIDTH+1)
#define	BFLYSHIFT	0
#define	LGWIDTH	(FFT_LGWIDTH)
#ifdef	DBLCLKFFT
#define	LGSPAN	(LGWIDTH-2)
#else
#define	LGSPAN	(LGWIDTH-1)
#endif
#define	ROUND	true

#define	SPANLEN		(1<<LGSPAN)
#define	SPANMASK	(SPANLEN-1)
#define	DBLSPANLEN	(1<<(LGSPAN+4))
#define	DBLSPANMASK	(DBLSPANLEN-1)

const	bool	gbl_debug = false;

class	FFTSTAGE_TB {
public:
	Vfftstage	*m_ftstage;
	VerilatedVcdC	*m_trace;
	long		m_iaddr;
	long		m_vals[SPANLEN], m_out[DBLSPANLEN];
	bool		m_syncd, m_ib_syncd, m_ob_syncd, m_input_sync;
	int		m_offset, m_ib_offset, m_ob_offset;
	uint64_t	m_tickcount;

	FFTSTAGE_TB(void) {
		Verilated::traceEverOn(true);
		m_ftstage = new Vfftstage;
		m_syncd = false;
		m_input_sync = false;
		m_ib_syncd   = false;
		m_ob_syncd   = false;
		m_iaddr = 0;
		m_offset = 0;
		m_tickcount = 0;
assert(OWIDTH == IWIDTH+1);
	}

	void	opentrace(const char *vcdname) {
		if (!m_trace) {
			m_trace = new VerilatedVcdC;
			m_ftstage->trace(m_trace, 99);
			m_trace->open(vcdname);
		}
	}

	void	closetrace(void) {
		if (m_trace) {
			m_trace->close();
			delete	m_trace;
			m_trace = NULL;
		}
	}

	void	tick(void) {
		m_tickcount++;

		m_ftstage->i_clk = 0;
		m_ftstage->eval();
		if (m_trace)	m_trace->dump((vluint64_t)(10ul*m_tickcount-2));
		m_ftstage->i_clk = 1;
		m_ftstage->eval();
		if (m_trace)	m_trace->dump((vluint64_t)(10ul*m_tickcount));
		m_ftstage->i_clk = 0;
		m_ftstage->eval();
		if (m_trace) {
			m_trace->dump((vluint64_t)(10ul*m_tickcount+5));
			m_trace->flush();
		}
	}

	void	cetick(void) {
		int	ce = m_ftstage->i_ce, nkce;

		tick();
		nkce = 0; // (rand()&1);
#ifdef	FFT_CKPCE
		nkce += FFT_CKPCE;
#endif
		if ((ce)&&(nkce > 0)) {
			m_ftstage->i_ce = 0;
			for(int kce = 1; kce < nkce; kce++)
				tick();
		}

		m_ftstage->i_ce = ce;
	}

	void	reset(void) {
		m_ftstage->i_ce  = 0;
		m_ftstage->i_reset = 1;
		m_ftstage->i_sync  = 0;
		tick();

		// Let's give it several ticks with no sync
		m_ftstage->i_ce = 0;
		m_ftstage->i_reset = 0;
		m_ftstage->i_sync  = 0;
		/*
		for(int i=0; i<8192; i++) {
			m_ftstage->i_data = rand();
			m_ftstage->i_sync = 0;
			m_ftstage->i_ce = 1;

			cetick();

			assert(m_ftstage->o_sync == 0);
		}
		*/

		m_iaddr = 0;
		m_offset = 0;
		m_syncd = false;
		m_ib_syncd = false;
		m_ob_syncd = false;
	}

	void	butterfly(const long cv, const long lft, const long rht,
				long &o_lft, long &o_rht) {
		long	cv_r, cv_i;
		long	lft_r, lft_i, rht_r, rht_i;
		long	o_lft_r, o_lft_i, o_rht_r, o_rht_i;

		cv_r = sbits(cv>>CWIDTH, CWIDTH);
		cv_i = sbits(cv, CWIDTH);

		lft_r = sbits(lft>>IWIDTH, IWIDTH);
		lft_i = sbits(lft, IWIDTH);

		rht_r = sbits(rht>>IWIDTH, IWIDTH);
		rht_i = sbits(rht, IWIDTH);

		o_lft_r = lft_r + rht_r;
		o_lft_i = lft_i + rht_i;
		o_lft_r = ubits(o_lft_r, OWIDTH);
		o_lft_i = ubits(o_lft_i, OWIDTH);

		// o_lft_r >>= 1;
		// o_lft_i >>= 1;
		o_lft = (o_lft_r << OWIDTH) | (o_lft_i);

		o_rht_r = (cv_r * (lft_r-rht_r)) - (cv_i * (lft_i-rht_i));
		o_rht_i = (cv_r * (lft_i-rht_i)) + (cv_i * (lft_r-rht_r));

		if (ROUND) {
			/* Non-convergent rounding
			if (o_rht_r & (1<<(CWIDTH-3)))
				o_rht_r += (1<<(CWIDTH-3))-1;
			if (o_rht_i & (1<<(CWIDTH-3)))
				o_rht_i += (1<<(CWIDTH-3))-1;
			*/

			// Convergent rounding
			if (o_rht_r & (1<<(CWIDTH-2)))
				o_rht_r += (1<<(CWIDTH-3));
			else
				o_rht_r += (1<<(CWIDTH-3))-1;

			if (o_rht_i & (1<<(CWIDTH-2)))
				o_rht_i += (1<<(CWIDTH-3));
			else
				o_rht_i += (1<<(CWIDTH-3))-1;
		}

		o_rht_r >>= (CWIDTH-2);
		o_rht_i >>= (CWIDTH-2);

		o_rht_r = ubits(o_rht_r, OWIDTH);
		o_rht_i = ubits(o_rht_i, OWIDTH);
		o_rht = (o_rht_r << OWIDTH) | (o_rht_i);

		/*
		printf("BUTTERFLY: %10lx %10lx %10lx -> %10lx %10lx\n",
			ubits(cv, 2*CWIDTH),
			ubits(lft,   2*IWIDTH), ubits(rht,   2*IWIDTH),
			ubits(o_lft, 2*OWIDTH), ubits(o_rht, 2*OWIDTH));
		*/
	}

	void	test(bool i_sync, long i_data) {
		long	cv;
		bool	bc;
		int	raddr;
		bool	failed = false;

		m_ftstage->i_reset = 0;
		m_ftstage->i_ce   = 1;
		m_ftstage->i_sync = i_sync;
		i_data = ubits(i_data, 2*IWIDTH);
		m_ftstage->i_data = i_data;

		if (!m_input_sync) {
			if (i_sync)
				m_input_sync = true;
			m_iaddr = 0;
		}

		cv = m_ftstage->cmem[m_iaddr & SPANMASK];
		bc = m_iaddr & (1<<LGSPAN);
		if (!bc)
			m_vals[m_iaddr & (SPANMASK)] = i_data;
		else {
			int	waddr = m_iaddr ^ (1<<LGSPAN);
			waddr &= (DBLSPANMASK);
			if (m_iaddr & (1<<(LGSPAN+1)))
				waddr |= (1<<(LGSPAN));
			butterfly(cv,
				m_vals[m_iaddr & (SPANMASK)],
				i_data,
				m_out[(m_iaddr-SPANLEN) & (DBLSPANMASK)],
				m_out[m_iaddr & (DBLSPANMASK)]);
		}

		cetick();

		unsigned	ib_addr = (m_iaddr - m_ib_offset)
					& ((SPANMASK<<1)|1);
		if ((!m_ib_syncd)&&(m_ftstage->ib_sync)) {
			m_ib_syncd = true;
			m_ib_offset = m_iaddr;

			if (gbl_debug) printf("IB-SYNC!!!!  Offset = %d\n", m_ib_offset);
		} else if ((m_ib_syncd)&&(ib_addr <(1<<(LGSPAN)))) {
			unsigned long	ib_a, ib_b, ib_c;
			ib_a = m_vals[ib_addr&SPANMASK];
			ib_b = i_data; // m_vals[(m_iaddr-m_ib_offset+(1<<(LGSPAN-1)))&(SPANMASK)];
			ib_c = m_ftstage->cmem[(m_iaddr-m_ib_offset)&(SPANMASK)];

			assert(m_ftstage->ib_a == ib_a);
			assert(m_ftstage->ib_b == ib_b);
			// assert(m_ftstage->ib_a == ib_a);
		}

		if ((!m_ob_syncd)&&(m_ftstage->ob_sync)) {
			m_ob_syncd = true;
			m_ob_offset = m_iaddr;

			if (gbl_debug) printf("OB-SYNC!!!!  Offset = %d\n", m_ob_offset);
		}

		if ((!m_syncd)&&(m_ftstage->o_sync)) {
			m_syncd = true;
			m_offset = m_iaddr;

			if (gbl_debug) printf("SYNC!!!!\n");
		}

		raddr = (m_iaddr-m_offset) & DBLSPANMASK;

		if (gbl_debug)
			printf("%4ld, %4ld: %d %9lx -> %d %9lx ... "
				"%4x %15lx (%10lx)",
			(long)m_iaddr, (long)raddr,
			// Input values
			i_sync, ubits(i_data, 2*IWIDTH),
			// Output values
			m_ftstage->o_sync,
			(long)m_ftstage->o_data,

			m_ftstage->iaddr&FFTMASK,
			(long)(ubits(m_ftstage->cmem[
				m_ftstage->iaddr&(SPANMASK)],
				2*CWIDTH)),
			(long)m_out[raddr]);

		unsigned long	oba, obb;

		oba = ubits(
		m_ftstage->VVAR(_HWBFLY__DOT__bfly__DOT__rnd_left_r),OWIDTH)
				<<OWIDTH;
		oba |= ubits(
		m_ftstage->VVAR(_HWBFLY__DOT__bfly__DOT__rnd_left_i), OWIDTH);

		obb = ubits(
		m_ftstage->VVAR(_HWBFLY__DOT__bfly__DOT__rnd_right_r), OWIDTH)
				<<OWIDTH;
		obb |= ubits(
		m_ftstage->VVAR(_HWBFLY__DOT__bfly__DOT__rnd_right_i), OWIDTH);

		if ((gbl_debug)&&(m_ob_syncd))
			printf(" [%d %10lx %10lx]",
				m_ftstage->ob_sync,
				ubits(oba, 2*OWIDTH), ubits(obb, 2*OWIDTH));

		if ((m_syncd)&&(m_ftstage->o_sync != ((((m_iaddr-m_offset)&((1<<(LGSPAN+1))-1))==0)?1:0))) {
			fprintf(stderr, "Bad output sync (m_iaddr = %lx, m_offset = %x)\n",
				(m_iaddr-m_offset) & SPANMASK, m_offset);
			failed = true;
		} if (gbl_debug) printf("\n");

		if (m_syncd) {
			if (m_out[raddr] != m_ftstage->o_data) {
				printf("Bad output data, ([%lx - %x = %x (%d)] "
					"%lx(exp) != %lx(sut))\n",
					m_iaddr, m_offset, raddr, raddr,
					m_out[raddr], (long)m_ftstage->o_data);
				failed = true;
			} // printf("Checked #%d\n", raddr);
		} else if (m_iaddr > FFTSIZE/2+128) {
			printf("NO OUTPUT SYNC!\n");
			failed = true;
		}
		m_iaddr++;

		if (failed)
			exit(-1);
	}
};
		


int	main(int argc, char **argv, char **envp) {
	Verilated::commandArgs(argc, argv);
	FFTSTAGE_TB	*ftstage = new FFTSTAGE_TB;
#ifdef	DBLCLKFFT
#define	STEP	2
#else
#define	STEP	1
#endif

	// ftstage->opentrace("fftstage.vcd");
	ftstage->reset();

	// Medium real (constant) value ... just for starters
	for(int k=0; k<FFTSIZE; k+=STEP)
		ftstage->test((k==0), 0x00200000l);
	// Medium imaginary (constant) value ... just for starters
	for(int k=0; k<FFTSIZE; k+=STEP)
		ftstage->test((k==0), 0x00000020l);
	// Medium sine wave, real
	for(int k=0; k<FFTSIZE; k+=STEP) {
		long vl;
		vl= (long)(cos(2.0 * M_PI * 1.0 / FFTSIZE * k)*(1l<<30) + 0.5);
		vl &= (-1l << 16); // Turn off the imaginary bit portion
		vl = ubits(vl,2*IWIDTH); // Turn off unused high order bits
		ftstage->test((k==1), vl);
	}
	// Smallest real value
	for(int k=0; k<FFTSIZE; k+=STEP)
		ftstage->test((k==0), 0x00080000l);
	// Smallest imaginary value
	for(int k=0; k<FFTSIZE; k+=STEP)
		ftstage->test((k==0), 0x00000001l);
	// Largest real value
	for(int k=0; k<FFTSIZE; k+=STEP)
		ftstage->test((k==0), 0x200000000l);
	// Largest negative imaginary value
	for(int k=0; k<FFTSIZE; k+=STEP)
		ftstage->test((k==0), 0x000010000l);
	// Let's try an impulse
	for(int k=0; k<FFTSIZE; k+=STEP)
		ftstage->test((k==0), (k==0)?0x020000000l:0l);
	// Now, let's clear out the result
	for(int k=0; k<FFTSIZE; k+=STEP)
		ftstage->test((k==0), 0x000000000l);
	for(int k=0; k<FFTSIZE; k+=STEP)
		ftstage->test((k==0), 0x000000000l);
	for(int k=0; k<FFTSIZE; k+=STEP)
		ftstage->test((k==0), 0x000000000l);
	for(int k=0; k<FFTSIZE; k+=STEP)
		ftstage->test((k==0), 0x000000000l);

	printf("SUCCESS! (Offset = %d)\n", ftstage->m_offset);
	delete	ftstage;

	exit(EXIT_SUCCESS);
}
