////////////////////////////////////////////////////////////////////////////
//
// Filename: 	mpy_tb.cpp
//
// Project:	A Doubletime Pipelined FFT
//
// Purpose:	A test-bench for the shift and add shiftaddmpy.v subfile of
//		the double clocked FFT.  This file may be run autonomously. 
//		If so, the last line output will either read "SUCCESS" on
//		success, or some other failure message otherwise.
//
//		This file depends upon verilator to both compile, run, and
//		therefore test shiftaddmpy.v
//
// Creator:	Dan Gisselquist, Ph.D.
//		Gisselquist Technology, LLC
//
///////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2015, Gisselquist Technology, LLC
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

#include "verilated.h"
#include "twoc.h"

class	MPYTB {
public:
	Vmpy	*mpy;
	long	vals[32];
	int	m_addr;

	MPYTB(void) {
		mpy = new Vmpy;

		for(int i=0; i<32; i++)
			vals[i] = 0;
		m_addr = 0;
	}
	~MPYTB(void) {
		delete mpy;
	}

	void	tick(void) {
		mpy->i_clk = 0;
		mpy->eval();
		mpy->i_clk = 1;
		mpy->eval();
	}

	void	reset(void) {
		mpy->i_clk = 0;
		mpy->i_ce = 1;
		mpy->i_a = 0;
		mpy->i_b = 0;

		for(int k=0; k<20; k++)
			tick();
	}

	bool	test(const int ia, const int ib) {
		bool	success;
		long	a, b, out;

		a = sbits(ia, AW);
		b = sbits(ib, BW);
		mpy->i_ce = 1;
		mpy->i_a = ubits(a, AW);
		mpy->i_b = ubits(b, BW);

		vals[m_addr&31] = a * b;

		tick();
		if (rand()&1) {
			mpy->i_ce = 0;
			tick();
		}

		printf("k=%3d: A =%04x, B =%05x -> O = %9lx (ANS=%10lx)\n",
			m_addr, (int)ubits(a,AW), (int)ubits(b,BW),
			(long)mpy->o_r, ubits(vals[m_addr&31], AW+BW+4));

		out = sbits(mpy->o_r, AW+BW);

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

	for(int k=0; k<2048; k++) {
		int	a, b, out;

		tb->test(rand(), rand());
	}

	delete	tb;

	printf("SUCCESS!\n");
	exit(0);
}
