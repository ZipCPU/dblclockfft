////////////////////////////////////////////////////////////////////////////////
//
// Filename: 	bldstage.cpp
//
// Project:	A General Purpose Pipelined FFT Implementation
//
// Purpose:	
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
// with this program.  (It's in the $(ROOT)/doc directory, run make with no
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
#define _CRT_SECURE_NO_WARNINGS   //  ms vs 2012 doesn't like fopen
#include <stdio.h>
#include <stdlib.h>

#ifdef _MSC_VER //  added for ms vs compatibility

#include <io.h>
#include <direct.h>
#define _USE_MATH_DEFINES

#else
// And for G++/Linux environment

#include <unistd.h>	// Defines the R_OK/W_OK/etc. macros
#endif

#include <string.h>
#include <string>
#include <math.h>
#include <ctype.h>
#include <assert.h>

#include "defaults.h"
#include "legal.h"
#include "fftlib.h"
#include "rounding.h"
#include "bldstage.h"

void	build_dblstage(const char *fname, ROUND_T rounding,
			const bool async_reset, const bool dbg) {
	FILE	*fp = fopen(fname, "w");
	if (NULL == fp) {
		fprintf(stderr, "Could not open \'%s\' for writing\n", fname);
		perror("O/S Err was:");
		return;
	}

	const	char	*rnd_string;
	if (rounding == RND_TRUNCATE)
		rnd_string = "truncate";
	else if (rounding == RND_FROMZERO)
		rnd_string = "roundfromzero";
	else if (rounding == RND_HALFUP)
		rnd_string = "roundhalfup";
	else
		rnd_string = "convround";

	std::string	resetw("i_reset");
	if (async_reset)
		resetw = std::string("i_areset_n");


	fprintf(fp,
SLASHLINE
"//\n"
"// Filename:\tlaststage%s.v\n"
"//\n"
"// Project:\t%s\n"
"//\n"
"// Purpose:\tThis is part of an FPGA implementation that will process\n"
"//		the final stage of a decimate-in-frequency FFT, running\n"
"//	through the data at two samples per clock.  If you notice from the\n"
"//	derivation of an FFT, the only time both even and odd samples are\n"
"//	used at the same time is in this stage.  Therefore, other than this\n"
"//	stage and these twiddles, all of the other stages can run two stages\n"
"//	at a time at one sample per clock.\n"
"//\n"
"// Operation:\n"
"// 	Given a stream of values, operate upon them as though they were\n"
"// 	value pairs, x[2n] and x[2n+1].  The stream begins when n=0, and ends\n"
"// 	when n=1.  When the first x[0] value enters, the synchronization\n"
"//	input, i_sync, must be true as well.\n"
"//\n"
"// 	For this stream, produce outputs\n"
"// 	y[2n  ] = x[2n] + x[2n+1], and\n"
"// 	y[2n+1] = x[2n] - x[2n+1]\n"
"//\n"
"// 	When y[0] is output, a synchronization bit o_sync will be true as\n"
"// 	well, otherwise it will be zero.\n"
"//\n"
"//\n"
"//	In this implementation, the output is valid one clock after the input\n"
"//	is valid.  The output also accumulates one bit above and beyond the\n"
"//	number of bits in the input.\n"
"//\n"
"//		i_clk	A system clock\n", (dbg)?"_dbg":"", prjname);
	if (async_reset)
		fprintf(fp,
"//		i_areset_n	An active low asynchronous reset\n");
	else
		fprintf(fp,
"//		i_reset	A synchronous reset\n");

	fprintf(fp,
"//		i_ce	Circuit enable--nothing happens unless this line is high\n"
"//		i_sync	A synchronization signal, high once per FFT at the start\n"
"//		i_left	The first (even) complex sample input.  The higher order\n"
"//			bits contain the real portion, low order bits the\n"
"//			imaginary portion, all in two\'s complement.\n"
"//		i_right	The next (odd) complex sample input, same format as\n"
"//			i_left.\n"
"//		o_left	The first (even) complex output.\n"
"//		o_right	The next (odd) complex output.\n"
"//		o_sync	Output synchronization signal.\n"
"//\n%s"
"//\n", creator);

	fprintf(fp, "%s", cpyleft);
	fprintf(fp, "//\n//\n`default_nettype\tnone\n//\n");
	fprintf(fp,
"module\tlaststage%s(i_clk, %s, i_ce, i_sync, i_left, i_right, o_left, o_right, o_sync%s);\n"
	"\tparameter\tIWIDTH=%d,OWIDTH=IWIDTH+1, SHIFT=%d;\n"
	"\tinput\t\ti_clk, %s, i_ce, i_sync;\n"
	"\tinput\t\t[(2*IWIDTH-1):0]\ti_left, i_right;\n"
	"\toutput\treg\t[(2*OWIDTH-1):0]\to_left, o_right;\n"
	"\toutput\treg\t\t\to_sync;\n"
	"\n", (dbg)?"_dbg":"", resetw.c_str(), (dbg)?", o_dbg":"",
	TST_DBLSTAGE_IWIDTH, TST_DBLSTAGE_SHIFT,
		resetw.c_str());

	if (dbg) { fprintf(fp, "\toutput\twire\t[33:0]\t\t\to_dbg;\n"
		"\tassign\to_dbg = { ((o_sync)&&(i_ce)), i_ce, o_left[(2*OWIDTH-1):(2*OWIDTH-16)],\n"
			"\t\t\t\t\to_left[(OWIDTH-1):(OWIDTH-16)] };\n"
"\n");
	}
	fprintf(fp,
	"\twire\tsigned\t[(IWIDTH-1):0]\ti_in_0r, i_in_0i, i_in_1r, i_in_1i;\n"
	"\tassign\ti_in_0r = i_left[(2*IWIDTH-1):(IWIDTH)];\n"
	"\tassign\ti_in_0i = i_left[(IWIDTH-1):0];\n"
	"\tassign\ti_in_1r = i_right[(2*IWIDTH-1):(IWIDTH)];\n"
	"\tassign\ti_in_1i = i_right[(IWIDTH-1):0];\n"
	"\twire\t[(OWIDTH-1):0]\t\to_out_0r, o_out_0i,\n"
				"\t\t\t\t\to_out_1r, o_out_1i;\n"
"\n"
"\n"
	"\t// Handle a potential rounding situation, when IWIDTH>=OWIDTH.\n"
"\n"
"\n");
	fprintf(fp,
	"\n"
	"\t// As with any register connected to the sync pulse, these must\n"
	"\t// have initial values and be reset on the %s signal.\n"
	"\t// Other data values need only restrict their updates to i_ce\n"
	"\t// enabled clocks, but sync\'s must obey resets and initial\n"
	"\t// conditions as well.\n"
	"\treg\trnd_sync, r_sync;\n"
"\n"
	"\tinitial\trnd_sync      = 1\'b0; // Sync into rounding\n"
	"\tinitial\tr_sync        = 1\'b0; // Sync coming out\n",
		resetw.c_str());
	if (async_reset)
		fprintf(fp, "\talways @(posedge i_clk, negdge i_areset_n)\n\t\tif (!i_areset_n)\n");
	else
		fprintf(fp, "\talways @(posedge i_clk)\n\t\tif (i_reset)\n");
	fprintf(fp,
		"\t\tbegin\n"
			"\t\t\trnd_sync <= 1\'b0;\n"
			"\t\t\tr_sync <= 1\'b0;\n"
		"\t\tend else if (i_ce)\n"
		"\t\tbegin\n"
			"\t\t\trnd_sync <= i_sync;\n"
			"\t\t\tr_sync <= rnd_sync;\n"
		"\t\tend\n"
"\n"
	"\t// As with other variables, these are really only updated when in\n"
	"\t// the processing pipeline, after the first i_sync.  However, to\n"
	"\t// eliminate as much unnecessary logic as possible, we toggle\n"
	"\t// these any time the i_ce line is enabled, and don\'t reset.\n"
	"\t// them on %s.\n", resetw.c_str());
	fprintf(fp,
	"\t// Don't forget that we accumulate a bit by adding two values\n"
	"\t// together. Therefore our intermediate value must have one more\n"
	"\t// bit than the two originals.\n"
	"\treg\tsigned\t[(IWIDTH):0]\trnd_in_0r, rnd_in_0i;\n"
	"\treg\tsigned\t[(IWIDTH):0]\trnd_in_1r, rnd_in_1i;\n\n"
	"\talways @(posedge i_clk)\n"
		"\t\tif (i_ce)\n"
		"\t\tbegin\n"
			"\t\t\t//\n"
			"\t\t\trnd_in_0r <= i_in_0r + i_in_1r;\n"
			"\t\t\trnd_in_0i <= i_in_0i + i_in_1i;\n"
			"\t\t\t//\n"
			"\t\t\trnd_in_1r <= i_in_0r - i_in_1r;\n"
			"\t\t\trnd_in_1i <= i_in_0i - i_in_1i;\n"
			"\t\t\t//\n"
		"\t\tend\n"
"\n");
	fprintf(fp,
	"\t%s #(IWIDTH+1,OWIDTH,SHIFT) do_rnd_0r(i_clk, i_ce,\n"
	"\t\t\t\t\t\t\trnd_in_0r, o_out_0r);\n\n", rnd_string);
	fprintf(fp,
	"\t%s #(IWIDTH+1,OWIDTH,SHIFT) do_rnd_0i(i_clk, i_ce,\n"
	"\t\t\t\t\t\t\trnd_in_0i, o_out_0i);\n\n", rnd_string);
	fprintf(fp,
	"\t%s #(IWIDTH+1,OWIDTH,SHIFT) do_rnd_1r(i_clk, i_ce,\n"
	"\t\t\t\t\t\t\trnd_in_1r, o_out_1r);\n\n", rnd_string);
	fprintf(fp,
	"\t%s #(IWIDTH+1,OWIDTH,SHIFT) do_rnd_1i(i_clk, i_ce,\n"
	"\t\t\t\t\t\t\trnd_in_1i, o_out_1i);\n\n", rnd_string);

	fprintf(fp, "\n"
	"\t// Prior versions of this routine did not include the extra\n"
	"\t// clock and register/flip-flops that this routine requires.\n"
	"\t// These are placed in here to correct a bug in Verilator, that\n"
	"\t// otherwise struggles.  (Hopefully this will fix the problem ...)\n"
	"\talways @(posedge i_clk)\n"
		"\t\tif (i_ce)\n"
		"\t\tbegin\n"
			"\t\t\to_left  <= { o_out_0r, o_out_0i };\n"
			"\t\t\to_right <= { o_out_1r, o_out_1i };\n"
		"\t\tend\n"
"\n"
	"\tinitial\to_sync = 1\'b0; // Final sync coming out of module\n");
	if (async_reset)
		fprintf(fp, "\talways @(posedge i_clk, negdge i_areset_n)\n\t\tif (!i_areset_n)\n");
	else
		fprintf(fp, "\talways @(posedge i_clk)\n\t\tif (i_reset)\n");
	fprintf(fp,
		"\t\t\to_sync <= 1\'b0;\n"
		"\t\telse if (i_ce)\n"
		"\t\t\to_sync <= r_sync;\n"
"\n"
"endmodule\n");
	fclose(fp);
}

void	build_stage(const char *fname,
		int stage, int nwide, int offset,
		int nbits, bool inv, int xtra, int ckpce,
		const bool async_reset, const bool dbg) {
	FILE	*fstage = fopen(fname, "w");
	int	cbits = nbits + xtra;

	std::string	resetw("i_reset");
	if (async_reset)
		resetw = std::string("i_areset_n");

	if (((unsigned)cbits * 2u) >= sizeof(long long)*8) {
		fprintf(stderr, "ERROR: CMEM Coefficient precision requested overflows long long data type.\n");
		exit(-1);
	}

	if (fstage == NULL) {
		fprintf(stderr, "ERROR: Could not open %s for writing!\n", fname);
		perror("O/S Err was:");
		fprintf(stderr, "Attempting to continue, but this file will be missing.\n");
		return;
	}

	fprintf(fstage,
SLASHLINE
"//\n"
"// Filename:\t%sfftstage%s.v\n"
"//\n"
"// Project:\t%s\n"
"//\n"
"// Purpose:\tThis file is (almost) a Verilog source file.  It is meant to\n"
"//		be used by a FFT core compiler to generate FFTs which may be\n"
"//	used as part of an FFT core.  Specifically, this file encapsulates\n"
"//	the options of an FFT-stage.  For any 2^N length FFT, there shall be\n"
"//	(N-1) of these stages.\n"
"//\n"
"//\n"
"// Operation:\n"
"// 	Given a stream of values, operate upon them as though they were\n"
"// 	value pairs, x[n] and x[n+N/2].  The stream begins when n=0, and ends\n"
"// 	when n=N/2-1 (i.e. there's a full set of N values).  When the value\n"
"// 	x[0] enters, the synchronization input, i_sync, must be true as well.\n"
"//\n"
"// 	For this stream, produce outputs\n"
"// 	y[n    ] = x[n] + x[n+N/2], and\n"
"// 	y[n+N/2] = (x[n] - x[n+N/2]) * c[n],\n"
"// 			where c[n] is a complex coefficient found in the\n"
"// 			external memory file COEFFILE.\n"
"// 	When y[0] is output, a synchronization bit o_sync will be true as\n"
"// 	well, otherwise it will be zero.\n"
"//\n"
"// 	Most of the work to do this is done within the butterfly, whether the\n"
"// 	hardware accelerated butterfly (uses a DSP) or not.\n"
"//\n%s"
"//\n",
		(inv)?"i":"", (dbg)?"_dbg":"", prjname, creator);
	fprintf(fstage, "%s", cpyleft);
	fprintf(fstage, "//\n//\n`default_nettype\tnone\n//\n");
	fprintf(fstage, "module\t%sfftstage%s(i_clk, %s, i_ce, i_sync, i_data, o_data, o_sync%s);\n",
		(inv)?"i":"", (dbg)?"_dbg":"", resetw.c_str(),
		(dbg)?", o_dbg":"");
	// These parameter values are useless at this point--they are to be
	// replaced by the parameter values in the calling program.  Only
	// problem is, the CWIDTH needs to match exactly!
	fprintf(fstage, "\tparameter\tIWIDTH=%d,CWIDTH=%d,OWIDTH=%d;\n",
		nbits, cbits, nbits+1);
	fprintf(fstage,
"\t// Parameters specific to the core that should be changed when this\n"
"\t// core is built ... Note that the minimum LGSPAN (the base two log\n"
"\t// of the span, or the base two log of the current FFT size) is 3.\n"
"\t// Smaller spans (i.e. the span of 2) must use the dbl laststage module.\n"
"\tparameter\tLGWIDTH=%d, LGSPAN=%d, LGBDLY=%d, BFLYSHIFT=0;\n"
"\tparameter\t[0:0]	OPT_HWMPY = 1;\n",
		lgval(stage), (nwide <= 1) ? lgval(stage)-1 : lgval(stage)-2,
		lgdelay(nbits,cbits-nbits));
	fprintf(fstage,
"\tparameter\t\tMPYDELAY = %d;\n", bflydelay(nbits,cbits-nbits));
	fprintf(fstage,
"\t// Clocks per CE.  If your incoming data rate is less than 50%% of your\n"
"\t// clock speed, you can set CKPCE to 2\'b10, make sure there's at least\n"
"\t// one clock between cycles when i_ce is high, and then use two\n"
"\t// multiplies instead of three.  Setting CKPCE to 2\'b11, and insisting\n"
"\t// on at least two clocks with i_ce low between cycles with i_ce high,\n"
"\t// then the hardware optimized butterfly code will used one multiply\n"
"\t// instead of two.\n"
"\tparameter\t[1:0]	CKPCE = 2'h%d;\n", ckpce);

	fprintf(fstage,
"\t// The COEFFILE parameter contains the name of the file containing the\n"
"\t// FFT twiddle factors\n");
	if (nwide == 2) {
		fprintf(fstage, "\tparameter\tCOEFFILE=\"%scmem_%c%d.hex\";\n",
			(inv)?"i":"", (offset)?'o':'e', stage*2);
	} else
		fprintf(fstage, "\tparameter\tCOEFFILE=\"%scmem_%d.hex\";\n",
			(inv)?"i":"", stage);
	fprintf(fstage,
"\tinput					i_clk, %s, i_ce, i_sync;\n"
"\tinput		[(2*IWIDTH-1):0]	i_data;\n"
"\toutput	reg	[(2*OWIDTH-1):0]	o_data;\n"
"\toutput	reg				o_sync;\n"
"\n", resetw.c_str());
	if (dbg) { fprintf(fstage, "\toutput\twire\t[33:0]\t\t\to_dbg;\n"
		"\tassign\to_dbg = { ((o_sync)&&(i_ce)), i_ce, o_data[(2*OWIDTH-1):(2*OWIDTH-16)],\n"
			"\t\t\t\t\to_data[(OWIDTH-1):(OWIDTH-16)] };\n"
"\n");
	}
	fprintf(fstage,
"\treg	wait_for_sync;\n"
"\treg	[(2*IWIDTH-1):0]	ib_a, ib_b;\n"
"\treg	[(2*CWIDTH-1):0]	ib_c;\n"
"\treg	ib_sync;\n"
"\n"
"\treg	b_started;\n"
"\twire	ob_sync;\n"
"\twire	[(2*OWIDTH-1):0]\tob_a, ob_b;\n");
	fprintf(fstage,
"\n"
"\t// cmem is defined as an array of real and complex values,\n"
"\t// where the top CWIDTH bits are the real value and the bottom\n"
"\t// CWIDTH bits are the imaginary value.\n"
"\t//\n"
"\t// cmem[i] = { (2^(CWIDTH-2)) * cos(2*pi*i/(2^LGWIDTH)),\n"
"\t//		(2^(CWIDTH-2)) * sin(2*pi*i/(2^LGWIDTH)) };\n"
"\t//\n"
"\treg	[(2*CWIDTH-1):0]	cmem [0:((1<<LGSPAN)-1)];\n"
"\tinitial\t$readmemh(COEFFILE,cmem);\n\n");

	// gen_coeff_file(coredir, fname, stage, cbits, nwide, offset, inv);

	fprintf(fstage,
"\treg	[(LGSPAN):0]		iaddr;\n"
"\treg	[(2*IWIDTH-1):0]	imem	[0:((1<<LGSPAN)-1)];\n"
"\n"
"\treg	[LGSPAN:0]		oB;\n"
"\treg	[(2*OWIDTH-1):0]	omem	[0:((1<<LGSPAN)-1)];\n"
"\n"
"\tinitial wait_for_sync = 1\'b1;\n"
"\tinitial iaddr = 0;\n");
	if (async_reset)
		fprintf(fstage, "\talways @(posedge i_clk, negedge i_areset_n)\n\t\tif (!i_areset_n)\n");
	else
		fprintf(fstage, "\talways @(posedge i_clk)\n\t\tif (i_reset)\n");

	fprintf(fstage,
	"\tbegin\n"
		"\t\t\twait_for_sync <= 1\'b1;\n"
		"\t\t\tiaddr <= 0;\n"
	"\tend else if ((i_ce)&&((!wait_for_sync)||(i_sync)))\n"
	"\tbegin\n"
		"\t\t//\n"
		"\t\t// First step: Record what we\'re not ready to use yet\n"
		"\t\t//\n"
		"\t\tiaddr <= iaddr + { {(LGSPAN){1\'b0}}, 1\'b1 };\n"
		"\t\twait_for_sync <= 1\'b0;\n"
	"\tend\n"
	"\talways @(posedge i_clk) // Need to make certain here that we don\'t read\n"
	"\tif ((i_ce)&&(!iaddr[LGSPAN])) // and write the same address on\n"
		"\t\timem[iaddr[(LGSPAN-1):0]] <= i_data; // the same clk\n"
	"\n");

	fprintf(fstage,
	"\t//\n"
	"\t// Now, we have all the inputs, so let\'s feed the butterfly\n"
	"\t//\n"
	"\tinitial ib_sync = 1\'b0;\n");
	if (async_reset)
		fprintf(fstage, "\talways @(posedge i_clk, negedge i_areset_n)\n\t\tif (!i_areset_n)\n");
	else
		fprintf(fstage, "\talways @(posedge i_clk)\n\t\tif (i_reset)\n");
	fprintf(fstage,
			"\t\tib_sync <= 1\'b0;\n"
		"\telse if (i_ce)\n"
		"\tbegin\n"
			"\t\t// Set the sync to true on the very first\n"
				"\t\t// valid input in, and hence on the very\n"
				"\t\t// first valid data out per FFT.\n"
				"\t\tib_sync <= (iaddr==(1<<(LGSPAN)));\n"
			"\t\t\tend\n"
	"\talways\t@(posedge i_clk)\n"
	"\tif (i_ce)\n"
	"\tbegin\n"
		"\t\t\t// One input from memory, ...\n"
		"\tib_a <= imem[iaddr[(LGSPAN-1):0]];\n"
		"\t\t// One input clocked in from the top\n"
		"\t\tib_b <= i_data;\n"
		"\t\t// and the coefficient or twiddle factor\n"
		"\t\tib_c <= cmem[iaddr[(LGSPAN-1):0]];\n"
	"\tend\n\n");

	fprintf(fstage,
"\tgenerate if (OPT_HWMPY)\n"
"\tbegin : HWBFLY\n"
"\t\thwbfly #(.IWIDTH(IWIDTH),.CWIDTH(CWIDTH),.OWIDTH(OWIDTH),\n"
			"\t\t\t\t.CKPCE(CKPCE), .SHIFT(BFLYSHIFT))\n"
		"\t\t\tbfly(i_clk, %s, i_ce, ib_c,\n"
			"\t\t\t\tib_a, ib_b, ib_sync, ob_a, ob_b, ob_sync);\n"
"\tend else begin : FWBFLY\n"
"\t\tbutterfly #(.IWIDTH(IWIDTH),.CWIDTH(CWIDTH),.OWIDTH(OWIDTH),\n"
		"\t\t\t\t.LGDELAY(LGBDLY),.MPYDELAY(MPYDELAY),\n"
		"\t\t\t\t.CKPCE(CKPCE),.SHIFT(BFLYSHIFT))\n"
	"\t\t\tbfly(i_clk, %s, i_ce, ib_c,\n"
		"\t\t\t\tib_a, ib_b, ib_sync, ob_a, ob_b, ob_sync);\n"
"\tend endgenerate\n\n",
			resetw.c_str(), resetw.c_str());

	fprintf(fstage,
	"\t//\n"
	"\t// Next step: recover the outputs from the butterfly\n"
	"\t//\n"
	"\tinitial oB        = 0;\n"
	"\tinitial o_sync    = 0;\n"
	"\tinitial b_started = 0;\n");
	if (async_reset)
		fprintf(fstage, "\talways @(posedge i_clk, negedge i_areset_n)\n\t\tif (!i_areset_n)\n");
	else
		fprintf(fstage, "\talways @(posedge i_clk)\n\t\tif (i_reset)\n");
	fprintf(fstage,
	"\t\tbegin\n"
		"\t\t\toB <= 0;\n"
		"\t\t\to_sync <= 0;\n"
		"\t\t\tb_started <= 0;\n"
	"\t\tend else if (i_ce)\n"
	"\t\tbegin\n"
	"\t\t\to_sync <= (!oB[LGSPAN])?ob_sync : 1\'b0;\n"
	"\t\t\tif (ob_sync||b_started)\n"
		"\t\t\t\toB <= oB + { {(LGSPAN){1\'b0}}, 1\'b1 };\n"
	"\t\t\tif ((ob_sync)&&(!oB[LGSPAN]))\n"
		"\t\t\t// A butterfly output is available\n"
			"\t\t\t\tb_started <= 1\'b1;\n"
	"\t\tend\n\n");
	fprintf(fstage,
	"\treg	[(LGSPAN-1):0]\t\tdly_addr;\n"
	"\treg	[(2*OWIDTH-1):0]\tdly_value;\n"
	"\talways @(posedge i_clk)\n"
	"\t\tif (i_ce)\n"
	"\t\tbegin\n"
	"\t\t\tdly_addr <= oB[(LGSPAN-1):0];\n"
	"\t\t\tdly_value <= ob_b;\n"
	"\t\tend\n"
	"\talways @(posedge i_clk)\n"
	"\t\tif (i_ce)\n"
		"\t\t\tomem[dly_addr] <= dly_value;\n"
"\n");
	fprintf(fstage,
	"\talways @(posedge i_clk)\n"
	"\t\tif (i_ce)\n"
	"\t\t\to_data <= (!oB[LGSPAN])?ob_a : omem[oB[(LGSPAN-1):0]];\n"
"\n");
	fprintf(fstage, "endmodule\n");
}
