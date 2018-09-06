////////////////////////////////////////////////////////////////////////////////
//
// Filename: 	mpy_tb.cpp
//
// Project:	A General Purpose Pipelined FFT Implementation
//
// Purpose:	A test-bench for the shift and add shiftaddmpy.v subfile of
//		the double clocked FFT.  This file may be run autonomously. 
//	If so, the last line output will either read "SUCCESS" on success, or
//	some other failure message otherwise.
//
//	This file depends upon verilator to both compile, run, and therefore
//	test shiftaddmpy.v
//
// Creator:	Dan Gisselquist, Ph.D.
//		Gisselquist Technology, LLC
//
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2015-2018, Gisselquist Technology, LLC
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

#ifdef	USE_OLD_MULTIPLY
#include "Vshiftaddmpy.h"
typedef	Vshiftaddmpy	Vmpy;
#define	AW	TST_SHIFTADDMPY_AW
#define	BW	TST_SHIFTADDMPY_BW
#define	DELAY	(TST_SHIFTADDMPY_AW+2)
#else
#include "Vlongbimpy.h"
typedef	Vlongbimpy	Vmpy;
#define	AW	TST_LONGBIMPY_AW
#define	BW	TST_LONGBIMPY_BW
#define	DELAY	((AW/2)+(AW&1)+2)
#endif

#include "twoc.h"

class	MPYTB {
public:
	Vmpy		*m_mpy;
	VerilatedVcdC	*m_trace;
	long		vals[32];
	int		m_addr;
	uint64_t	m_tickcount;

	MPYTB(void) {
		Verilated::traceEverOn(true);
		m_mpy = new Vmpy;
		m_tickcount = 0;

		for(int i=0; i<32; i++)
			vals[i] = 0;
		m_addr = 0;
	}
	~MPYTB(void) {
		closetrace();
		delete m_mpy;
	}

	void    opentrace(const char *vcdname) {
		if (!m_trace) {
			m_trace = new VerilatedVcdC;
			m_mpy->trace(m_trace, 99);
			m_trace->open(vcdname);
		}
	}

        void    closetrace(void) {
                if (m_trace) {
                        m_trace->close();
                        delete  m_trace;
                        m_trace = NULL;
                }
        }

        void    tick(void) {
                m_tickcount++;

		m_mpy->i_clk = 0;
		m_mpy->eval();
		if (m_trace)	m_trace->dump((uint64_t)(10ul*m_tickcount-2));
		m_mpy->i_clk = 1;
		m_mpy->eval();
		if (m_trace)	m_trace->dump((uint64_t)(10ul*m_tickcount));
		m_mpy->i_clk = 0;
		m_mpy->eval();
		if (m_trace)	{
			m_trace->dump((uint64_t)(10ul*m_tickcount+5));
			m_trace->flush();
		}
	}

	void	cetick(void) {
		int	ce = m_mpy->i_ce, nkce;

		tick();
		nkce = (rand()&1);
#ifdef	FFT_CKPCE
		nkce += FFT_CKPCE;
#endif
		if ((ce)&&(nkce>0)) {
			m_mpy->i_ce = 0;
			for(int kce=1; kce<nkce; kce++)
				tick();
		}

		m_mpy->i_ce = ce;
	}

	void	reset(void) {
		m_mpy->i_clk = 0;
		m_mpy->i_ce = 1;
		m_mpy->i_a_unsorted = 0;
		m_mpy->i_b_unsorted = 0;

		for(int k=0; k<20; k++)
			cetick();
	}

	bool	test(const int ia, const int ib) {
		bool	success;
		long	a, b, out;

		a = sbits(ia, AW);
		b = sbits(ib, BW);
		m_mpy->i_ce = 1;
		m_mpy->i_a_unsorted = ubits(a, AW);
		m_mpy->i_b_unsorted = ubits(b, BW);

		vals[m_addr&31] = a * b;

		cetick();

		printf("k=%3d: A =%0*x, B =%0*x -> O = %*lx (ANS=%*lx)\n",
			m_addr, (AW+3)/4, (int)ubits(a,AW),
			(BW+3)/4, (int)ubits(b,BW),
			(AW+BW+3)/4, (long)m_mpy->o_r,
			(AW+BW+7)/4, ubits(vals[m_addr&31], AW+BW+4));

		out = sbits(m_mpy->o_r, AW+BW);

		m_addr++;

		success = (m_addr < (DELAY+2))||(out == vals[(m_addr-DELAY)&31]);
		if (!success) {
			printf("WRONG ANSWER: %8lx (exp) != %8lx (sut)\n", vals[(m_addr-DELAY)&0x01f], out);
			exit(-1);
		}
		
		return success;
	}
};

int	main(int argc, char **argv, char **envp) {
	Verilated::commandArgs(argc, argv);
	MPYTB		*tb = new MPYTB;

	// tb->opentrace("mpy.vcd");
	tb->reset();

	for(int k=0; k<15; k++) {
		int	a, b;

		a = (1<<k);
		b = 1;
		tb->test(a, b);
	}

	for(int k=0; k<15; k++) {
		int	a, b, out;

		a = (1<<15);
		b = (1<<k);
		tb->test(a, b);
	}

	if (AW+BW <= 20) {
		// Exhaustive test
		for(int a=0; a< (1<<AW); a++)
		for(int b=0; b< (1<<BW); b++)
			tb->test(a, b);
		printf("Exhaust complete\n");
	} else {
		// Pseudorandom test
		for(int k=0; k<2048; k++)
			tb->test(rand(), rand());
	}

	delete	tb;

	printf("SUCCESS!\n");
	exit(0);
}
