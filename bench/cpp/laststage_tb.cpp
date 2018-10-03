////////////////////////////////////////////////////////////////////////////////
//
// Filename: 	laststage_tb.cpp
//
// Project:	A General Purpose Pipelined FFT Implementation
//
// Purpose:	A test-bench for the laststage.v subfile of the general purpose
//		pipelined FFT.  This file may be run autonomously.  If so,
//	the last line output will either read "SUCCESS" on success, or some
//	other failure message otherwise.
//
//	This file depends upon verilator to both compile, run, and therefore
//	test laststage.v
//
// Creator:	Dan Gisselquist, Ph.D.
//		Gisselquist Technology, LLC
//
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2015,2018 Gisselquist Technology, LLC
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
#include <stdio.h>
#include <stdint.h>

#include "verilated.h"
#include "verilated_vcd_c.h"
#include "Vlaststage.h"
#include "fftsize.h"
#include "twoc.h"

#define	IWIDTH	16
#define	OWIDTH	(IWIDTH+1)
#define	SHIFT	0
#define	ROUND	1

#define	ASIZ	32
#define	AMSK	(ASIZ-1)

class	LASTSTAGE_TB {
public:
	Vlaststage	*m_last;
	VerilatedVcdC	*m_trace;
#ifdef	DBLCLKFFT
	unsigned long	m_left[ASIZ], m_right[ASIZ];
#else
	unsigned long	m_data[ASIZ];
#endif
	bool		m_syncd;
	int		m_addr, m_offset;
	unsigned long	m_tickcount;

	LASTSTAGE_TB(void) {
		Verilated::traceEverOn(true);
		m_last = new Vlaststage;
		m_tickcount = 0;
		m_syncd = false; m_addr = 0, m_offset = 0;
	}

	void	opentrace(const char *vcdname) {
		if (!m_trace) {
			m_trace = new VerilatedVcdC;
			m_last->trace(m_trace, 99);
			m_trace->open(vcdname);
		}
	}

	void	closetrace(void) {
		if (m_trace) {
			m_trace->close();
			delete m_trace;
			m_trace = NULL;
		}
	}

	void	tick(void) {
		m_tickcount++;

		m_last->i_clk = 0;
		m_last->eval();
		if (m_trace) m_trace->dump((vluint64_t)(10ul * m_tickcount - 2));
		m_last->i_clk = 1;
		m_last->eval();
		if (m_trace) m_trace->dump((vluint64_t)(10ul * m_tickcount));
		m_last->i_clk = 0;
		m_last->eval();
		if (m_trace) {
			m_trace->dump((vluint64_t)(10ul * m_tickcount + 5));
			m_trace->flush();
		}
		m_last->i_reset  = 0;
		m_last->i_sync = 0;
	}

	void	cetick(void) {
		int	nkce;

		tick();
		nkce = (rand()&1);
#ifdef	FFT_CKPCE
		nkce += FFT_CKPCE;
#endif
		if ((m_last->i_ce)&&(nkce > 0)) {
			m_last->i_ce = 0;
			for(int kce = 1; kce < nkce; kce++)
				tick();
			m_last->i_ce = 1;
		}
	}

	void	reset(void) {
		m_last->i_reset = 1;
		tick();

		m_syncd = false; m_addr = 0, m_offset = 0;
	}

	void	check_results(void) {
		bool	failed = false;

		if ((!m_syncd)&&(m_last->o_sync)) {
			m_syncd = true;
			m_offset = m_addr;
			printf("SYNCD at %d\n", m_addr);
		}

#ifdef	DBLCLKFFT
		int	ir0, ir1, ii0, ii1, or0, oi0, or1, oi1;

		ir0 = sbits(m_left[ (m_addr-m_offset)&AMSK]>>IWIDTH, IWIDTH);
		ir1 = sbits(m_right[(m_addr-m_offset)&AMSK]>>IWIDTH, IWIDTH);
		ii0 = sbits(m_left[ (m_addr-m_offset)&AMSK], IWIDTH);
		ii1 = sbits(m_right[(m_addr-m_offset)&AMSK], IWIDTH);


		or0 = sbits(m_last->o_left  >> OWIDTH, OWIDTH);
		oi0 = sbits(m_last->o_left           , OWIDTH);
		or1 = sbits(m_last->o_right >> OWIDTH, OWIDTH);
		oi1 = sbits(m_last->o_right          , OWIDTH);


		// Sign extensions
		printf("k=%3d: IN = %08x:%08x, OUT =%09lx:%09lx, S=%d\n",
			m_addr, m_last->i_left, m_last->i_right,
			m_last->o_left, m_last->o_right,
			m_last->o_sync);

		/*
		printf("\tI0 = { %x : %x }, I1 = { %x : %x }, O0 = { %x : %x }, O1 = { %x : %x }\n",
			ir0, ii0, ir1, ii1, or0, oi0, or1, oi1);
		*/

		if (m_syncd) {
			if (or0 != (ir0 + ir1))	{
				printf("FAIL 1: or0 != (ir0+ir1), or %x(exp) != %x(sut)\n", (ir0+ir1), or0);
				failed=true;}
			if (oi0 != (ii0 + ii1))	{printf("FAIL 2\n"); failed=true;}
			if (or1 != (ir0 - ir1))	{printf("FAIL 3\n"); failed=true;}
			if (oi1 != (ii0 - ii1))	{printf("FAIL 4\n"); failed=true;}
		} else if (m_addr > 20) {
			printf("NO SYNC!\n");
			failed = true;
		}
#else
		int	or0, oi0;
		int	sumr, sumi, difr, difi;
		int	ir0, ii0, ir1, ii1, ir2, ii2, ir3, ii3, irn, iin;

		irn = sbits(m_data[(m_addr-m_offset+2)&AMSK]>>IWIDTH, IWIDTH);
		iin = sbits(m_data[(m_addr-m_offset+2)&AMSK], IWIDTH);
		ir0 = sbits(m_data[(m_addr-m_offset+1)&AMSK]>>IWIDTH, IWIDTH);
		ii0 = sbits(m_data[(m_addr-m_offset+1)&AMSK], IWIDTH);
		ir1 = sbits(m_data[(m_addr-m_offset  )&AMSK]>>IWIDTH, IWIDTH);
		ii1 = sbits(m_data[(m_addr-m_offset  )&AMSK], IWIDTH);
		ir2 = sbits(m_data[(m_addr-m_offset-1)&AMSK]>>IWIDTH, IWIDTH);
		ii2 = sbits(m_data[(m_addr-m_offset-1)&AMSK], IWIDTH);
		ir3 = sbits(m_data[(m_addr-m_offset-2)&AMSK]>>IWIDTH, IWIDTH);
		ii3 = sbits(m_data[(m_addr-m_offset-2)&AMSK], IWIDTH);

		sumr = ir1 + ir0;
		sumi = ii1 + ii0;

		difr = ir2 - ir1;
		difi = ii2 - ii1;

		or0 = sbits(m_last->o_val  >> OWIDTH, OWIDTH);
		oi0 = sbits(m_last->o_val           , OWIDTH);

		printf("IR0 = %08x, IR1 = %08x, IR2 = %08x, ",
				ir0, ir1, ir2);
		printf("II0 = %08x, II1 = %08x, II2 = %08x, ",
				ii0, ii1, ii2);
		// Sign extensions
		printf("k=%3d: IN = %08x, %c, OUT =%09lx, S=%d\n",
			m_addr, m_last->i_val,
			m_last->i_sync ? 'S':' ',
			m_last->o_val, m_last->o_sync);


		if ((m_syncd)&&(0 == ((m_addr-m_offset)&1))) {
			if (or0 != sumr)	{
				printf("FAIL 1: or0 != (ir0+ir1), or %x(exp) != %x(sut)\n", sumr, or0);
				failed=true;
			} if (oi0 != sumi)	{
				printf("FAIL 2\n");
				failed=true;
			}
		} else if ((m_syncd)&&(1 == ((m_addr-m_offset)&1))) {
			if (or0 != difr) {
				printf("FAIL 3: or0 != (ir1-ir0), or %x(exp) != %x(sut)\n", difr, or0);
				failed=true;
			} if (oi0 != difi) {
				printf("FAIL 4: oi0 != (ii1-ii0), or %x(exp) != %x(sut)\n", difi, oi0);
				failed=true;
			}
		} else if (m_addr > 20) {
			printf("NO SYNC!\n");
			failed = true;
		}
#endif
		if (failed)
			exit(-2);
	}

	void	sync(void) {
		m_last->i_sync = 1;
	}

	void	test(unsigned long left, unsigned long right) {
		m_last->i_ce    = 1;
		if (m_last->i_sync)
			m_addr = 0;
#ifdef	DBLCLKFFT
		m_last->i_left  = left;
		m_last->i_right = right;

		m_left[ m_addr&AMSK] = m_last->i_left;
		m_right[m_addr&AMSK] = m_last->i_right;
		m_addr++;

		cetick();
#else
		m_last->i_val = left;
		m_data[ m_addr&AMSK] = m_last->i_val;
		m_addr = (m_addr+1);
		cetick();

		check_results();

		m_last->i_val = right;
		m_data[m_addr&AMSK] = m_last->i_val;
		m_addr = (m_addr+1)&AMSK;
		cetick();
#endif

		check_results();
	}

	void	test(int ir0, int ii0, int ir1, int ii1) {
		unsigned long	left, right, mask = (1<<IWIDTH)-1;

		left  = ((ir0&mask) << IWIDTH) | (ii0 & mask);
		right = ((ir1&mask) << IWIDTH) | (ii1 & mask);
		test(left, right);
	}
};

int	main(int argc, char **argv, char **envp) {
	Verilated::commandArgs(argc, argv);
	LASTSTAGE_TB	*tb = new LASTSTAGE_TB;

	// tb->opentrace("laststage.vcd");
	tb->reset();

	tb->sync();

	tb->test( 1, 0,0,0);
	tb->test( 0, 2,0,0);
	tb->test( 0, 0,4,0);
	tb->test( 0, 0,0,8);

	tb->test( 0, 0,0,0);

	tb->test(16,16,0,0);
	tb->test(0,0,16,16);
	tb->test(16,-16,0,0);
	tb->test(0,0,16,-16);
	tb->test(16,16,0,0);
	tb->test(0,0,16,16);

	for(int k=0; k<64; k++) {
		int16_t	ir0, ii0, ir1, ii1;

		// Let's pick some random values, ...
		ir0 = rand(); if (ir0&4) ir0 = -ir0;
		ii0 = rand(); if (ii0&2) ii0 = -ii0;
		ir1 = rand(); if (ir1&1) ir1 = -ir1;
		ii1 = rand(); if (ii1&8) ii1 = -ii1;

		tb->test(ir0, ii0, ir1, ii1);

	}

	delete	tb;

	printf("SUCCESS!\n");
	exit(0);
}






