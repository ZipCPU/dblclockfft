////////////////////////////////////////////////////////////////////////////////
//
// Filename: 	qtrstage_tb.cpp
//
// Project:	A General Purpose Pipelined FFT Implementation
//
// Purpose:	A test-bench for the qtrstage.v subfile of the double
//		clocked FFT.  This file may be run autonomously.  If so,
//	the last line output will either read "SUCCESS" on success, or some
//	other failure message otherwise.
//
//	This file depends upon verilator to both compile, run, and therefore
//	test qtrstage.v
//
// Creator:	Dan Gisselquist, Ph.D.
//		Gisselquist Technology, LLC
//
////////////////////////////////////////////////////////////////////////////////
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
#include "Vqtrstage.h"
#include "twoc.h"
#include "fftsize.h"

#define	IWIDTH	TST_QTRSTAGE_IWIDTH
#define	OWIDTH	(IWIDTH+1)
#define	LGWIDTH	TST_QTRSTAGE_LGWIDTH

#define	ASIZ	32
#define	AMSK	(ASIZ-1)

#ifdef	NEW_VERILATOR
#define	VVAR(A)	qtrstage__DOT_ ## A
#else
#define	VVAR(A)	v__DOT_ ## A
#endif

#define	sum_r		VVAR(_sum_r)
#define	sum_i		VVAR(_sum_i)
#define	diff_r		VVAR(_diff_r)
#define	diff_i		VVAR(_diff_i)
#define	pipeline	VVAR(_pipeline)
#define	iaddr		VVAR(_iaddr)
#define	imem		VVAR(_imem)
#define	wait_for_sync	VVAR(_wait_for_sync)

class	QTRTEST_TB {
public:
	Vqtrstage	*m_qstage;
	VerilatedVcdC	*m_trace;
	unsigned long	m_data[ASIZ], m_tickcount;
	int		m_addr, m_offset;
	bool		m_syncd;

	QTRTEST_TB(void) {
		Verilated::traceEverOn(true);
		m_trace = NULL;
		m_qstage = new Vqtrstage;
		m_addr = 0;
		m_offset = 6;
		m_syncd = false;
		m_tickcount = 0;
	}

	void	opentrace(const char *vcdname) {
		if (!m_trace) {
			m_trace = new VerilatedVcdC;
			m_qstage->trace(m_trace, 99);
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

		m_qstage->i_clk = 0;
		m_qstage->eval();
		if (m_trace)	m_trace->dump((vluint64_t)(10ul*m_tickcount-2));
		m_qstage->i_clk = 1;
		m_qstage->eval();
		if (m_trace)	m_trace->dump((vluint64_t)(10ul*m_tickcount));
		m_qstage->i_clk = 0;
		m_qstage->eval();
		if (m_trace) {
			m_trace->dump((vluint64_t)(10ul*m_tickcount+5));
			m_trace->flush();
		}

		m_qstage->i_sync = 0;
	}

	void	cetick(void) {
		int	nkce;

		tick();
		nkce = (rand()&1);
#ifdef	FFT_CKPCE
		nkce += FFT_CKPCE;
#endif
		if ((m_qstage->i_ce)&&(nkce>0)) {
			m_qstage->i_ce = 0;
			for(int kce = 1; kce < nkce; kce++)
				tick();
			m_qstage->i_ce = 1;
		}
	}

	void	reset(void) {
		m_qstage->i_ce  = 0;
		m_qstage->i_reset = 1;
		tick();
		m_qstage->i_ce  = 0;
		m_qstage->i_reset = 0;
		tick();

		m_addr = 0; m_offset = 6; m_syncd = false;
	}

	void	check_results(void) {
		int	sumr, sumi, difr, difi, or0, oi0;
		bool	fail = false;

		if ((!m_syncd)&&(m_qstage->o_sync)) {
			m_syncd = true;
			assert(m_addr == m_offset);
			m_offset = m_addr;
			printf("VALID-SYNC!!\n");
		}

		if (!m_syncd)
			return;

#ifdef	DBLCLKFFT
		int	ir0, ii0, ir1, ii1, ir2, ii2;

		ir0 = sbits(m_data[(m_addr-m_offset-1)&AMSK]>>IWIDTH, IWIDTH);
		ii0 = sbits(m_data[(m_addr-m_offset-1)&AMSK], IWIDTH);
		ir1 = sbits(m_data[(m_addr-m_offset  )&AMSK]>>IWIDTH, IWIDTH);
		ii1 = sbits(m_data[(m_addr-m_offset  )&AMSK], IWIDTH);
		ir2 = sbits(m_data[(m_addr-m_offset+1)&AMSK]>>IWIDTH, IWIDTH);
		ii2 = sbits(m_data[(m_addr-m_offset+1)&AMSK], IWIDTH);

		sumr = ir1 + ir2;
		sumi = ii1 + ii2;
		difr = ir0 - ir1;
		difi = ii0 - ii1;

		or0 = sbits(m_qstage->o_data >> OWIDTH, OWIDTH);
		oi0 = sbits(m_qstage->o_data, OWIDTH);

		if (0==((m_addr-m_offset)&1)) {
			if (or0 != sumr)	{
				printf("FAIL 1: or0 != sumr (%x(exp) != %x(sut))\n", sumr, or0); fail = true;}
			if (oi0 != sumi)	{
				printf("FAIL 2: oi0 != sumi (%x(exp) != %x(sut))\n", sumi, oi0); fail = true;}
		} else if (1==((m_addr-m_offset)&1)) {
			if (or0 != difr)	{
				printf("FAIL 3: or0 != difr (%x(exp) != %x(sut))\n", difr, or0); fail = true;}
			if (oi0 != difi)	{
				printf("FAIL 4: oi0 != difi (%x(exp) != %x(sut))\n", difi, oi0); fail = true;}
		}
#else
		int	locn = (m_addr-m_offset)&AMSK;
		int	ir1, ii1, ir3, ii3, ir5, ii5;

		ir5 = sbits(m_data[(m_addr-m_offset-2)&AMSK]>>IWIDTH, IWIDTH);
		ii5 = sbits(m_data[(m_addr-m_offset-2)&AMSK], IWIDTH);
		ir3 = sbits(m_data[(m_addr-m_offset  )&AMSK]>>IWIDTH, IWIDTH);
		ii3 = sbits(m_data[(m_addr-m_offset  )&AMSK], IWIDTH);
		ir1 = sbits(m_data[(m_addr-m_offset+2)&AMSK]>>IWIDTH, IWIDTH);
		ii1 = sbits(m_data[(m_addr-m_offset+2)&AMSK], IWIDTH);

		sumr = ir3 + ir1;
		sumi = ii3 + ii1;
		difr = ir5 - ir3;
		difi = ii5 - ii3;

		or0 = sbits(m_qstage->o_data >> OWIDTH, OWIDTH);
		oi0 = sbits(m_qstage->o_data, OWIDTH);

		if (0==((locn)&2)) {
			if (or0 != sumr)	{
				printf("FAIL 1: or0 != sumr (%x(exp) != %x(sut))\n", sumr, or0); fail = true;
			}
			if (oi0 != sumi)	{
				printf("FAIL 2: oi0 != sumi (%x(exp) != %x(sut))\n", sumi, oi0); fail = true;}
		} else if (2==((m_addr-m_offset)&3)) {
			if (or0 != difr)	{
				printf("FAIL 3: or0 != difr (%x(exp) != %x(sut))\n", difr, or0); fail = true;}
			if (oi0 != difi)	{
				printf("FAIL 4: oi0 != difi (%x(exp) != %x(sut))\n", difi, oi0); fail = true;}
		} else if (3==((m_addr-m_offset)&3)) {
			if (or0 != difi)	{
				printf("FAIL 3: or0 != difr (%x(exp) != %x(sut))\n", difr, or0); fail = true;}
			if (oi0 != -difr)	{
				printf("FAIL 4: oi0 != difi (%x(exp) != %x(sut))\n", difi, oi0); fail = true;}
		}

//		if (m_qstage->o_sync != ((((m_addr-m_offset)&127) == 0)?1:0)) {
//			printf("BAD O-SYNC, m_addr = %d, m_offset = %d\n", m_addr, m_offset); fail = true;
//		}
#endif


		if (fail)
			exit(-1);
	}

	void	sync(void) {
		m_qstage->i_sync = 1;
		m_addr = 0;
	}

	void	test(unsigned int data) {
		int	isync = m_qstage->i_sync;
		m_qstage->i_ce = 1;
		m_qstage->i_data = data;
		// m_qstage->i_sync = (((m_addr&127)==2)?1:0);
		// printf("DATA[%08x] = %08x ... ", m_addr, data);
		m_data[ (m_addr++)&AMSK] = data;
		tick();

		printf("k=%4d: ISYNC=%d, IN = %08x, OUT =%09lx, SYNC=%d\t%5x,%5x,%5x,%5x\t%x %4x %8x %d\n",
			(m_addr-m_offset), isync, m_qstage->i_data,
			m_qstage->o_data, m_qstage->o_sync,

			m_qstage->sum_r,
			m_qstage->sum_i,
			m_qstage->diff_r,
			m_qstage->diff_i,
			m_qstage->pipeline,
			m_qstage->iaddr,
#ifdef	DBLCLKFFT
			m_qstage->imem,
#else
			m_qstage->imem[1],
#endif
			m_qstage->wait_for_sync);

		check_results();
	}

	void	test(int ir0, int ii0) {
		unsigned int	data;

		data = (((ir0&((1<<IWIDTH)-1)) << IWIDTH) | (ii0 & ((1<<IWIDTH)-1)));
		// printf("%d,%d -> %8x\n", ir0, ii0, data);
		test(data);
	}

	void	random_test(void) {
		int	ir0, ii0;

		// Let's pick some random values
		ir0 = rand(); if (ir0&4) ir0 = -ir0;
		ii0 = rand(); if (ii0&2) ii0 = -ii0;
		test(ir0, ii0);
	}
};

int	main(int argc, char **argv, char **envp) {
	Verilated::commandArgs(argc, argv);
	QTRTEST_TB	*tb = new QTRTEST_TB;
	int16_t		ir0, ii0, ir1, ii1, ir2, ii2;
	int32_t		sumr, sumi, difr, difi;

	// tb->opentrace("qtrstage.vcd");
	tb->reset();

	tb->test( 16, 0);
	tb->test( 16, 0);
	tb->sync();

	tb->test(  8,  0);
	tb->test(  0,  0);
	tb->test(  0,  0);
	tb->test(  0,  0);

	tb->test(  0, 4);
	tb->test(  0,  0);
	tb->test(  0,  0);
	tb->test(  0,  0);

	tb->test(  0,  0);
	tb->test( 32,  0);
	tb->test(  0,  0);
	tb->test(  0,  0);

	tb->test(  0,  0);
	tb->test(  0, 64);
	tb->test(  0,  0);
	tb->test(  0,  0);

	tb->test(  0,  0);
	tb->test(  0,  0);
	tb->test(128,  0);
	tb->test(  0,  0);

	tb->test(  0,  0);
	tb->test(  0,  0);
	tb->test(  0,256);
	tb->test(  0,  0);

	tb->test(  0,  0);
	tb->test(  0,  0);
	tb->test(  0,  0);
	tb->test(  2,  0);

	tb->test(  0,  0);
	tb->test(  0,  0);
	tb->test(  0,  0);
	tb->test(  0,  1);

	tb->test(  0, 16);
	tb->test(  0, 16);
	tb->test( 16,  0);
	tb->test(-16,  0);

	for(int k=0; k<1060; k++) {
		tb->random_test();
	}

	delete	tb;

	printf("SUCCESS!\n");
	exit(0);
}





