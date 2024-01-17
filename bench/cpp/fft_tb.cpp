////////////////////////////////////////////////////////////////////////////////
//
// Filename: 	fft_tb.cpp
// {{{
// Project:	A General Purpose Pipelined FFT Implementation
//
// Purpose:	A test-bench for the main program, fftmain.v, of the double
//		clocked FFT.  This file may be run autonomously  (when
//	fully functional).  If so, the last line output will either read
//	"SUCCESS" on success, or some other failure message otherwise.
//
//	This file depends upon verilator to both compile, run, and therefore
//	test fftmain.v
//
// Creator:	Dan Gisselquist, Ph.D.
//		Gisselquist Technology, LLC
//
////////////////////////////////////////////////////////////////////////////////
// }}}
// Copyright (C) 2015-2024, Gisselquist Technology, LLC
// {{{
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
// }}}
// License:	GPL, v3, as defined and found on www.gnu.org,
// {{{
//		http://www.gnu.org/licenses/gpl.html
//
////////////////////////////////////////////////////////////////////////////////
// }}}
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>

#include "verilated.h"
#include "verilated_vcd_c.h"
#include "Vfftmain.h"
#include "twoc.h"

#include "fftsize.h"


#ifdef	ROOT_VERILATOR

#include "Vfftmain___024root.h"

#define	VVAR(A)	rootp->fftmain__DOT_ ## A

#elif	defined(NEW_VERILATOR)
#define	VVAR(A)	fftmain__DOT_ ## A
#else
#define	VVAR(A)	v__DOT_ ## A
#endif

#ifdef	DBLCLKFFT
#define	revstage_iaddr	VVAR(_revstage__DOT__iaddr)
#else
#define	revstage_iaddr	VVAR(_revstage__DOT__wraddr)
#endif
#define	br_sync		VVAR(_br_sync)
#define	br_started	VVAR(_r_br_started)
#define	w_s2048		VVAR(_w_s2048)
#define	w_s1024		VVAR(_w_s1024)
#define	w_s512		VVAR(_w_s512)
#define	w_s256		VVAR(_w_s256)
#define	w_s128		VVAR(_w_s128)
#define	w_s64		VVAR(_w_s64)
#define	w_s32		VVAR(_w_s32)
#define	w_s16		VVAR(_w_s16)
#define	w_s8		VVAR(_w_s8)
#define	w_s4		VVAR(_w_s4)


#define	IWIDTH	FFT_IWIDTH
#define	OWIDTH	FFT_OWIDTH
#define	LGWIDTH	FFT_LGWIDTH

#if (IWIDTH > 16)
typedef	unsigned long	ITYP;
#else
typedef	unsigned int	ITYP;
#endif

#if (OWIDTH > 16)
typedef	unsigned long	OTYP;
#else
typedef	unsigned int	OTYP;
#endif

#define	NFTLOG	16
#define	FFTLEN	(1<<LGWIDTH)

#ifdef	FFT_SKIPS_BIT_REVERSE
#define	APPLY_BITREVERSE_LOCALLY
#endif

unsigned long bitrev(const int nbits, const unsigned long vl) {
	unsigned long	r = 0;
	unsigned long	val = vl;

	for(int k=0; k<nbits; k++) {
		r<<= 1;
		r |= (val & 1);
		val >>= 1;
	}

	return r;
}

class	FFT_TB {
public:
	Vfftmain	*m_fft;
	OTYP		m_data[FFTLEN];
	ITYP		m_log[NFTLOG*FFTLEN];
	int		m_iaddr, m_oaddr, m_ntest, m_logbase;
	FILE		*m_dumpfp;
	fftw_plan	m_plan;
	double		*m_fft_buf;
	bool		m_syncd;
	unsigned long	m_tickcount;
	VerilatedVcdC*	m_trace;

	FFT_TB(void) {
		m_fft = new Vfftmain;
		Verilated::traceEverOn(true);
		m_iaddr = m_oaddr = 0;
		m_dumpfp = NULL;

		m_fft_buf = (double *)fftw_malloc(sizeof(fftw_complex)*(FFTLEN));
		m_plan = fftw_plan_dft_1d(FFTLEN, (fftw_complex *)m_fft_buf,
				(fftw_complex *)m_fft_buf,
				FFTW_FORWARD, FFTW_MEASURE);
		m_syncd = false;
		m_ntest = 0;
	}

	~FFT_TB(void) {
		closetrace();
		delete m_fft;
		m_fft = NULL;
	}

	virtual void opentrace(const char *vcdname) {
		if (!m_trace) {
			m_trace = new VerilatedVcdC;
			m_fft->trace(m_trace, 99);
			m_trace->open(vcdname);
		}
	}

	virtual void closetrace(void) {
		if (m_trace) {
			m_trace->close();
			delete m_trace;
			m_trace = NULL;
		}
	}

	void	tick(void) {
		m_tickcount++;
		if (m_fft->i_reset)
			printf("TICK(%s,%s)\n",
				(m_fft->i_reset)?"RST":"   ",
				(m_fft->i_ce)?"CE":"  ");

		m_fft->i_clk = 0;
		m_fft->eval();
		if (m_trace)
			m_trace->dump((vluint64_t)(10*m_tickcount-2));
		m_fft->i_clk = 1;
		m_fft->eval();
		if (m_trace)
			m_trace->dump((vluint64_t)(10*m_tickcount));
		m_fft->i_clk = 0;
		m_fft->eval();
		if (m_trace) {
			m_trace->dump((vluint64_t)(10*m_tickcount+5));
			m_trace->flush();
		}
	}

	void	cetick(void) {
		int	ce = m_fft->i_ce, nkce;
		tick();

		nkce = (rand()&1);
#ifdef	FFT_CKPCE
		nkce += FFT_CKPCE;
#endif
		if ((ce)&&(nkce>0)) {
			m_fft->i_ce = 0;
			for(int kce=1; kce < nkce; kce++)
				tick();
		}

		m_fft->i_ce = ce;
	}

	void	reset(void) {
		m_fft->i_ce  = 0;
		m_fft->i_reset = 1;
		tick();
		m_fft->i_reset = 0;
		tick();

		m_iaddr = m_oaddr = m_logbase = 0;
		m_syncd = false;
		m_tickcount = 0l;
	}

	long	twos_complement(const long val, const int bits) {
		return sbits(val, bits);
	}

	void	checkresults(void) {
		double	*dp, *sp; // Complex array
		double	vout[FFTLEN*2];
		double	isq=0.0, osq = 0.0;
		ITYP	*lp;

		// Fill up our test array from the log array
		printf("%3d : CHECK: %8d %5x m_log[-%x=%x]\n", m_ntest, m_iaddr, m_iaddr,
			m_logbase, (m_iaddr-m_logbase)&((NFTLOG*FFTLEN-1)&(-FFTLEN)));

		// Convert our logged data into doubles, in an FFT buffer
		dp = m_fft_buf; lp = &m_log[(m_iaddr-m_logbase)&((NFTLOG*FFTLEN-1)&(-FFTLEN))];
		for(int i=0; i<FFTLEN; i++) {
			ITYP	tv = *lp++;

			dp[0] = sbits((long)tv >> IWIDTH, IWIDTH);
			dp[1] = sbits((long)tv, IWIDTH);

			// printf("IN[%4d = %4x] = %9.1f %9.1f\n",
				// i+((m_iaddr-FFTLEN*3)&((4*FFTLEN-1)&(-FFTLEN))),
				// i+((m_iaddr-FFTLEN*3)&((4*FFTLEN-1)&(-FFTLEN))),
				// dp[0], dp[1]);
			dp += 2;
		}

		// Let's measure ... are we the zero vector?  If not, how close?
		dp = m_fft_buf;
		for(int i=0; i<FFTLEN*2; i++) {
			isq += (*dp) * (*dp); dp++;
		}

		fftw_execute(m_plan);

		// Let's load up the output we received into double valued
		// array vout
		dp = vout;
		for(int i=0; i<FFTLEN; i++) {
			*dp = rdata(i);
			osq += (*dp) * (*dp); dp++;
			*dp = idata(i);
			osq += (*dp) * (*dp); dp++;
		}


		// Let's figure out if there's a scale factor difference ...
		double	scale = 0.0, wt = 0.0;
		sp = m_fft_buf;  dp = vout;
		for(int i=0; i<FFTLEN*2; i++) {
			scale += (*sp) * (*dp++);
			wt += (*sp) * (*sp); sp++;
		} scale = scale / wt;

		if (fabs(scale) <= 1./4./FFTLEN)
			scale = 2./(FFTLEN);
		else if (wt == 0.0) scale = 1.0;

		double xisq = 0.0;
		sp = m_fft_buf;  dp = vout;

		if ((true)&&(m_dumpfp)) {
			double	tmp[FFTLEN*2], nscl;

			if (fabs(scale) < 1e-4)
				nscl = 1.0;
			else
				nscl = scale;
			for(int i=0; i<FFTLEN*2; i++)
				tmp[i] = m_fft_buf[i] * nscl;
			fwrite(tmp, sizeof(double), FFTLEN*2, m_dumpfp);
		}

		for(int i=0; i<FFTLEN*2; i++) {
			double vl = (*sp++) * scale - (*dp++);
			xisq += vl * vl;
		}

		printf("%3d : SCALE = %12.6f, WT = %18.1f, ISQ = %15.1f, ",
			m_ntest, scale, wt, isq);
		printf("OSQ = %18.1f, ", osq);
		printf("XISQ = %18.1f, sqrt = %9.2f\n", xisq, sqrt(xisq));
		if (xisq > 1.4 * FFTLEN/2) {
			printf("TEST FAIL!!  Result is out of bounds from ");
			printf("expected result with FFTW3.\n");
			// exit(EXIT_FAILURE);
		}
		m_ntest++;
	}

#ifdef	DBLCLKFFT
	bool	test(ITYP lft, ITYP rht) {
		m_fft->i_ce    = 1;
		m_fft->i_reset = 0;
		m_fft->i_left  = lft;
		m_fft->i_right = rht;

		m_log[(m_iaddr++)&(NFTLOG*FFTLEN-1)] = lft;
		m_log[(m_iaddr++)&(NFTLOG*FFTLEN-1)] = rht;

		cetick();

		if (m_fft->o_sync) {
			if (!m_syncd) {
				m_syncd = true;
				printf("ORIGINAL SYNC AT 0x%lx, m_oaddr set to 0x%x\n", m_tickcount, m_oaddr);
				m_logbase = m_iaddr;
			} else printf("RESYNC AT %lx\n", m_tickcount);
			m_oaddr &= (-1<<LGWIDTH);
		} else m_oaddr += 2;

		printf("%8x,%5d: %08x,%08x -> %011lx,%011lx\t",
			m_iaddr, m_oaddr,
			lft, rht, m_fft->o_left, m_fft->o_right);

#ifndef	APPLY_BITREVERSE_LOCALLY
		printf(" [%3x]%s", m_fft->revstage_iaddr,
			(m_fft->br_sync)?"S"
				:((m_fft->br_started)?".":"x"));
#endif

		printf(" ");
#if (FFT_SIZE>=2048)
		printf("%s", (m_fft->w_s2048)?"S":"-");
#endif
#if (FFT_SIZE>1024)
		printf("%s", (m_fft->w_s1024)?"S":"-");
#endif
#if (FFT_SIZE>512)
		printf("%s", (m_fft->w_s512)?"S":"-");
#endif
#if (FFT_SIZE>256)
		printf("%s", (m_fft->w_s256)?"S":"-");
#endif
#if (FFT_SIZE>128)
		printf("%s", (m_fft->w_s128)?"S":"-");
#endif
#if (FFT_SIZE>64)
		printf("%s", (m_fft->w_s64)?"S":"-");
#endif
#if (FFT_SIZE>32)
		printf("%s", (m_fft->w_s32)?"S":"-");
#endif
#if (FFT_SIZE>16)
		printf("%s", (m_fft->w_s16)?"S":"-");
#endif
#if (FFT_SIZE>8)
		printf("%s", (m_fft->w_s8)?"S":"-");
#endif
#if (FFT_SIZE>4)
		printf("%s", (m_fft->w_s4)?"S":"-");
#endif

		printf(" %s%s\n",
			(m_fft->o_sync)?"\t(SYNC!)":"",
			(m_fft->o_left | m_fft->o_right)?"  (NZ)":"");

		m_data[(m_oaddr  )&(FFTLEN-1)] = m_fft->o_left;
		m_data[(m_oaddr+1)&(FFTLEN-1)] = m_fft->o_right;

		if ((m_syncd)&&((m_oaddr&(FFTLEN-1)) == FFTLEN-2)) {
			dumpwrite();
			checkresults();
		}

		return (m_fft->o_sync);
	}
#else
	bool	test(ITYP data) {
		m_fft->i_ce    = 1;
		m_fft->i_reset = 0;
		m_fft->i_sample  = data;

		m_log[(m_iaddr++)&(NFTLOG*FFTLEN-1)] = data;

		cetick();

		if (m_fft->o_sync) {
			if (!m_syncd) {
				m_syncd = true;
				printf("ORIGINAL SYNC AT 0x%lx, m_oaddr set to 0x%x\n", m_tickcount, m_oaddr);
				m_logbase = m_iaddr;
			} else printf("RESYNC AT %lx\n", m_tickcount);
			m_oaddr &= (-1<<LGWIDTH);
		} else m_oaddr += 1;

		printf("%8x,%5d: %08x -> %011lx\t",
			m_iaddr, m_oaddr, data, m_fft->o_result);

#ifndef	APPLY_BITREVERSE_LOCALLY
		printf(" [%3x]%s", m_fft->revstage_iaddr,
			(m_fft->br_sync)?"S"
				:((m_fft->br_started)?".":"x"));
#endif

		printf(" ");
#if (FFT_SIZE>=2048)
		printf("%s", (m_fft->w_s2048)?"S":"-");
#endif
#if (FFT_SIZE>1024)
		printf("%s", (m_fft->w_s1024)?"S":"-");
#endif
#if (FFT_SIZE>512)
		printf("%s", (m_fft->w_s512)?"S":"-");
#endif
#if (FFT_SIZE>256)
		printf("%s", (m_fft->w_s256)?"S":"-");
#endif
#if (FFT_SIZE>128)
		printf("%s", (m_fft->w_s128)?"S":"-");
#endif
#if (FFT_SIZE>64)
		printf("%s", (m_fft->w_s64)?"S":"-");
#endif
#if (FFT_SIZE>32)
		printf("%s", (m_fft->w_s32)?"S":"-");
#endif
#if (FFT_SIZE>16)
		printf("%s", (m_fft->w_s16)?"S":"-");
#endif
#if (FFT_SIZE>8)
		printf("%s", (m_fft->w_s8)?"S":"-");
#endif
#if (FFT_SIZE>4)
		printf("%s", (m_fft->w_s4)?"S":"-");
#endif

		printf(" %s%s\n",
			(m_fft->o_sync)?"\t(SYNC!)":"",
			(m_fft->o_result)?"  (NZ)":"");

		m_data[(m_oaddr  )&(FFTLEN-1)] = m_fft->o_result;

		if ((m_syncd)&&((m_oaddr&(FFTLEN-1)) == FFTLEN-1)) {
			dumpwrite();
			checkresults();
		}

		return (m_fft->o_sync);
	}
#endif

	bool	test(double lft_r, double lft_i, double rht_r, double rht_i) {
		ITYP	ilft, irht, ilft_r, ilft_i, irht_r, irht_i;

		ilft_r = (ITYP)(lft_r) & ((1<<IWIDTH)-1);
		ilft_i = (ITYP)(lft_i) & ((1<<IWIDTH)-1);
		irht_r = (ITYP)(rht_r) & ((1<<IWIDTH)-1);
		irht_i = (ITYP)(rht_i) & ((1<<IWIDTH)-1);

		ilft = (ilft_r << IWIDTH) | ilft_i;
		irht = (irht_r << IWIDTH) | irht_i;

#ifdef	DBLCLKFFT
		return test(ilft, irht);
#else
		test(ilft);
		return test(irht);
#endif
	}

	double	rdata(int addr) {
		int	index = addr & (FFTLEN-1);

#ifdef	APPLY_BITREVERSE_LOCALLY
		index = bitrev(LGWIDTH, index);
#endif
		return (double)sbits(m_data[index]>>OWIDTH, OWIDTH);
	}

	double	idata(int addr) {
		int	index = addr & (FFTLEN-1);

#ifdef	APPLY_BITREVERSE_LOCALLY
		index = bitrev(LGWIDTH, index);
#endif
		return (double)sbits(m_data[index], OWIDTH);
	}

	void	dump(FILE *fp) {
		m_dumpfp = fp;
	}

	void	dumpwrite(void) {
		if (!m_dumpfp)
			return;

		double	*buf;

		buf = new double[FFTLEN * 2];
		for(int i=0; i<FFTLEN; i++) {
			buf[i*2] = rdata(i);
			buf[i*2+1] = idata(i);
		}

		fwrite(buf, sizeof(double), FFTLEN*2, m_dumpfp);
		delete[] buf;
	}
};


int	main(int argc, char **argv, char **envp) {
	Verilated::commandArgs(argc, argv);
	FFT_TB *fft = new FFT_TB;
	FILE	*fpout;

	fpout = fopen("fft_tb.dbl", "w");
	if (NULL == fpout) {
		fprintf(stderr, "Cannot write output file, fft_tb.dbl\n");
		exit(-1);
	}

	fft->opentrace("fft.vcd");
	fft->reset();

	{
		int	ftlen = FFTLEN;
		fwrite(&ftlen, 1, sizeof(int), fpout);
	}
	fft->dump(fpout);

	// 1.
	double	maxv = ((1l<<(IWIDTH-1))-1l);
	fft->test(0.0, 0.0, maxv, 0.0);
	for(int k=0; k<FFTLEN/2-1; k++)
		fft->test(0.0,0.0,0.0,0.0);

	// 2. Try placing a pulse at the very end location
	for(int k=0; k<FFTLEN/2; k++) {
		double cl, cr, sl, sr, W;
		W = - 2.0 * M_PI / FFTLEN * (1);
		cl = cos(W * (2*k  )) * (double)((1l<<(IWIDTH-2))-1l);
		sl = sin(W * (2*k  )) * (double)((1l<<(IWIDTH-2))-1l);
		cr = cos(W * (2*k+1)) * (double)((1l<<(IWIDTH-2))-1l);
		sr = sin(W * (2*k+1)) * (double)((1l<<(IWIDTH-2))-1l);
		fft->test(cl, sl, cr, sr);
	}

	// 2. 
	fft->test(maxv, 0.0, maxv, 0.0);
	for(int k=0; k<FFTLEN/2-1; k++)
		fft->test(0.0,0.0,0.0,0.0);

	// 3. 
	fft->test(0.0,0.0,0.0,0.0);
	fft->test(maxv, 0.0, 0.0, 0.0);
	for(int k=0; k<FFTLEN/2-1; k++)
		fft->test(0.0,0.0,0.0,0.0);

	// 4.
	for(int k=0; k<8; k++)
		fft->test(maxv, 0.0, maxv, 0.0);
	for(int k=8; k<FFTLEN/2; k++)
		fft->test(0.0,0.0,0.0,0.0);

	// 5.
	if (FFTLEN/2 >= 16) {
		for(int k=0; k<16; k++)
			fft->test(maxv, 0.0, maxv, 0.0);
		for(int k=16; k<FFTLEN/2; k++)
			fft->test(0.0,0.0,0.0,0.0);
	}

	// 6.
	if (FFTLEN/2 >= 32) {
		for(int k=0; k<32; k++)
			fft->test(maxv, 0.0, maxv, 0.0);
		for(int k=32; k<FFTLEN/2; k++)
			fft->test(0.0,0.0,0.0,0.0);
	}

	// 7.
	if (FFTLEN/2 >= 64) {
		for(int k=0; k<64; k++)
			fft->test(maxv, 0.0, maxv, 0.0);
		for(int k=64; k<FFTLEN/2; k++)
			fft->test(0.0,0.0,0.0,0.0);
	}

	if (FFTLEN/2 >= 128) {
		for(int k=0; k<128; k++)
			fft->test(maxv, 0.0, maxv, 0.0);
		for(int k=128; k<FFTLEN/2; k++)
			fft->test(0.0,0.0,0.0,0.0);
	}

	if (FFTLEN/2 >= 256) {
		for(int k=0; k<256; k++)
			fft->test(maxv, 0.0, maxv, 0.0);
		for(int k=256; k<FFTLEN/2; k++)
			fft->test(0.0,0.0,0.0,0.0);
	}

	if (FFTLEN/2 >= 512) {
		for(int k=0; k<256+128; k++)
			fft->test(maxv, 0.0, maxv, 0.0);
		for(int k=256+128; k<FFTLEN/2; k++)
			fft->test(0.0,0.0,0.0,0.0);
	}

	/*
	for(int k=0; k<FFTLEN/2; k++)
		fft->test(0.0,0.0,0.0,0.0);

	for(int k=0; k<FFTLEN/2; k++)
		fft->test(0.0,0.0,0.0,0.0);

	for(int k=0; k<FFTLEN/2; k++)
		fft->test(0.0,0.0,0.0,0.0);

	for(int k=0; k<FFTLEN/2; k++)
		fft->test(0.0,0.0,0.0,0.0);

	for(int k=0; k<FFTLEN/2; k++)
		fft->test(0.0,0.0,0.0,0.0);

	for(int k=0; k<FFTLEN/2; k++)
		fft->test(0.0,0.0,0.0,0.0);
	*/

#ifndef	NO_JUNK
	// 7.

	//     1 -> 0x0001 
	//     2 -> 0x0002 
	//     4 -> 0x0004 
	//     8 -> 0x0008 
	//    16 -> 0x0010 
	//    32 -> 0x0020 
	//    64 -> 0x0040 
	//   128 -> 0x0080 
	//   256 -> 0x0100 
	//   512 -> 0x0200 
	//  1024 -> 0x0400 
	//  2048 -> 0x0800
	//  4096 -> 0x1000
	//  8192 -> 0x2000
	// 16384 -> 0x4000
	for(int v=1; v<32768; v<<=1) for(int k=0; k<FFTLEN/2; k++)
		fft->test((double)v,0.0,(double)v,0.0);
	//     1 -> 0xffff 	
	//     2 -> 0xfffe
	//     4 -> 0xfffc
	//     8 -> 0xfff8
	//    16 -> 0xfff0
	//    32 -> 0xffe0
	//    64 -> 0xffc0
	//   128 -> 0xff80
	//   256 -> 0xff00
	//   512 -> 0xfe00
	//  1024 -> 0xfc00
	//  2048 -> 0xf800
	//  4096 -> 0xf000
	//  8192 -> 0xe000
	// 16384 -> 0xc000
	// 32768 -> 0x8000
	fft->test(0.0,0.0,16384.0,0.0);
	for(int k=0; k<FFTLEN/2-1; k++)
		fft->test(0.0,0.0,0.0,0.0);

	for(int v=1; v<=32768; v<<=1) for(int k=0; k<FFTLEN/2; k++)
		fft->test(-(double)v,0.0,-(double)v,0.0);
	//     1 -> 0x000040 	CORRECT!!
	//     2 -> 0x000080 
	//     4 -> 0x000100 
	//     8 -> 0x000200
	//    16 -> 0x000400
	//    32 -> 0x000800
	//    64 -> 0x001000
	//   128 -> 0x002000
	//   256 -> 0x004000
	//   512 -> 0x008000
	//  1024 -> 0x010000
	//  2048 -> 0x020000
	//  4096 -> 0x040000
	//  8192 -> 0x080000
	// 16384 -> 0x100000
	for(int v=1; v<32768; v<<=1) for(int k=0; k<FFTLEN/2; k++)
		fft->test(0.0,(double)v,0.0,(double)v);
	//     1 -> 0x3fffc0
	//     2 -> 0x3fff80
	//     4 -> 0x3fff00
	//     8 -> 0x3ffe00
	//    16 -> 0x3ffc00
	//    32 -> 0x3ff800
	//    64 -> 0x3ff000
	//   128 -> 0x3fe000
	//   256 -> 0x3fc000
	//   512 -> 0x3f8000
	//  1024 -> 0x3f0000
	//  2048 -> 0x3e0000
	//  4096 -> 0x3c0000
	//  8192 -> 0x380000
	// 16384 -> 0x300000
	for(int v=1; v<32768; v<<=1) for(int k=0; k<FFTLEN/2; k++)
		fft->test(0.0,-(double)v,0.0,-(double)v);

	// 61. Now, how about the smallest alternating real signal
	for(int k=0; k<FFTLEN/2; k++)
		fft->test(2.0,0.0,0.0,0.0); // Don't forget to expect a bias!
	// 62. Now, how about the smallest alternating imaginary signal
	for(int k=0; k<FFTLEN/2; k++)
		fft->test(0.0,2.0,0.0,0.0); // Don't forget to expect a bias!
	// 63. Now, how about the smallest alternating real signal,2nd phase
	for(int k=0; k<FFTLEN/2; k++)
		fft->test(0.0,0.0,2.0,0.0); // Don't forget to expect a bias!
	// 64.Now, how about the smallest alternating imaginary signal,2nd phase
	for(int k=0; k<FFTLEN/2; k++)
		fft->test(0.0,0.0,0.0,2.0); // Don't forget to expect a bias!

	// 65.
	for(int k=0; k<FFTLEN/2; k++)
		fft->test(maxv,0.0,-maxv,0.0);
	// 66.
	for(int k=0; k<FFTLEN/2; k++)
		fft->test(0.0,-maxv,0.0,maxv);
	// 67.
	for(int k=0; k<FFTLEN/2; k++)
		fft->test(-maxv,-maxv,-maxv,-maxv);
	// 68.
	for(int k=0; k<FFTLEN/2; k++)
		fft->test(0.0,-maxv,0.0,maxv);
	// 69.
	for(int k=0; k<FFTLEN/2; k++)
		fft->test(0.0,maxv,0.0,-maxv);
	// 70. 
	for(int k=0; k<FFTLEN/2; k++)
		fft->test(-maxv,-maxv,-maxv,-maxv);

	// 71. Now let's go for an impulse (SUCCESS)
	fft->test(16384.0, 0.0, 0.0, 0.0);
	for(int k=0; k<FFTLEN/2-1; k++)
		fft->test(0.0,0.0,0.0,0.0);

	// 72. And another one on the next clock (FAILS, ugly)
	fft->test(0.0, 0.0, 16384.0, 0.0);
	for(int k=0; k<FFTLEN/2-1; k++)
		fft->test(0.0,0.0,0.0,0.0);

	// 72. And another one on the next clock (FAILS, ugly)
	fft->test(0.0, 0.0,  8192.0, 0.0);
	for(int k=0; k<FFTLEN/2-1; k++)
		fft->test(0.0,0.0,0.0,0.0);

	// 72. And another one on the next clock (FAILS, ugly)
	fft->test(0.0, 0.0,   512.0, 0.0);
	for(int k=0; k<FFTLEN/2-1; k++)
		fft->test(0.0,0.0,0.0,0.0);

	// 73. And an imaginary one on the second clock
	fft->test(0.0, 0.0, 0.0, 16384.0);
	for(int k=0; k<FFTLEN/2-1; k++)
		fft->test(0.0,0.0,0.0,0.0);

	// 74. Likewise the next clock
	fft->test(0.0,0.0,0.0,0.0);
	fft->test(16384.0, 0.0, 0.0, 0.0);
	for(int k=0; k<FFTLEN/2-2; k++)
		fft->test(0.0,0.0,0.0,0.0);

	// 75. And it's imaginary counterpart
	fft->test(0.0,0.0,0.0,0.0);
	fft->test(0.0, 16384.0, 0.0, 0.0);
	for(int k=0; k<FFTLEN/2-2; k++)
		fft->test(0.0,0.0,0.0,0.0);

	// 76. Likewise the next clock
	fft->test(0.0,0.0,0.0,0.0);
	fft->test(0.0, 0.0, 16384.0, 0.0);
	for(int k=0; k<FFTLEN/2-2; k++)
		fft->test(0.0,0.0,0.0,0.0);

	// 77. And it's imaginary counterpart
	fft->test(0.0,0.0,0.0,0.0);
	fft->test(0.0, 0.0, 0.0, 16384.0);
	for(int k=0; k<FFTLEN/2-2; k++)
		fft->test(0.0,0.0,0.0,0.0);


	// 78. Now let's try some exponentials
	for(int k=0; k<FFTLEN/2; k++) {
		double cl, cr, sl, sr, W;
		W = - 2.0 * M_PI / FFTLEN;
		cl = cos(W * (2*k  )) * 16383.0;
		sl = sin(W * (2*k  )) * 16383.0;
		cr = cos(W * (2*k+1)) * 16383.0;
		sr = sin(W * (2*k+1)) * 16383.0;
		fft->test(cl, sl, cr, sr);
	}

	// 72.
	for(int k=0; k<FFTLEN/2; k++) {
		double cl, cr, sl, sr, W;
		W = - 2.0 * M_PI / FFTLEN * 5;
		cl = cos(W * (2*k  )) * 16383.0;
		sl = sin(W * (2*k  )) * 16383.0;
		cr = cos(W * (2*k+1)) * 16383.0;
		sr = sin(W * (2*k+1)) * 16383.0;
		fft->test(cl, sl, cr, sr);
	}

	// 73.
	for(int k=0; k<FFTLEN/2; k++) {
		double cl, cr, sl, sr, W;
		W = - 2.0 * M_PI / FFTLEN * 8;
		cl = cos(W * (2*k  )) * 8190.0;
		sl = sin(W * (2*k  )) * 8190.0;
		cr = cos(W * (2*k+1)) * 8190.0;
		sr = sin(W * (2*k+1)) * 8190.0;
		fft->test(cl, sl, cr, sr);
	}

	// 74.
	for(int k=0; k<FFTLEN/2; k++) {
		double cl, cr, sl, sr, W;
		W = - 2.0 * M_PI / FFTLEN * 25;
		cl = cos(W * (2*k  )) * 4.0;
		sl = sin(W * (2*k  )) * 4.0;
		cr = cos(W * (2*k+1)) * 4.0;
		sr = sin(W * (2*k+1)) * 4.0;
		fft->test(cl, sl, cr, sr);
	}
#endif
	// 19.--24. And finally, let's clear out our results / buffer
	for(int k=0; k<(FFTLEN/2) * 5; k++)
		fft->test(0.0,0.0,0.0,0.0);



	fclose(fpout);

	if (!fft->m_syncd) {
		printf("FAIL -- NO SYNC\n");
		goto test_failure;
	}

	printf("SUCCESS!!\n");
	exit(0);
test_failure:
	printf("TEST FAILED!!\n");
	exit(0);
}


