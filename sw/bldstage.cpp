////////////////////////////////////////////////////////////////////////////////
//
// Filename: 	bldstage.cpp
// {{{
// Project:	A General Purpose Pipelined FFT Implementation
//
// Purpose:	Builds the logic necessary to implement a single stage of an
//		FFT.  This includes referencing the butterfly, but not the
//	actual butterflies themselves.  Further, this file only contains the
//	code for the general case of an FFT stage: the special cases of the
//	two final stages are described in other files.
//
// Creator:	Dan Gisselquist, Ph.D.
//		Gisselquist Technology, LLC
//
////////////////////////////////////////////////////////////////////////////////
// }}}
// Copyright (C) 2015-2021, Gisselquist Technology, LLC
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
//
////////////////////////////////////////////////////////////////////////////////
//
// }}}
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

// build_dblstage
// {{{
// Builds the penultimate FFT stage, using integer operations only.
// This stage is called laststage elsewhere.
//
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
"// {{{\n" // "}}}"
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
"module\tlaststage%s #(\n"
	"\t\t// {{{\n"
	// (i_clk, %s, i_ce, i_sync, i_left, i_right, o_left, o_right, o_sync%s);\n"
	"\t\tparameter\tIWIDTH=%d,OWIDTH=IWIDTH+1, SHIFT=%d\n"
	"\t\t// }}}\n"
	"\t) (\n"
	"\t\t// {{{\n"
	"\t\tinput\twire\ti_clk, %s, i_ce, i_sync,\n"
	"\t\tinput\twire\t[(2*IWIDTH-1):0]\ti_left, i_right,\n"
	"\t\toutput\treg\t[(2*OWIDTH-1):0]\to_left, o_right,\n"
	"\t\toutput\treg\t\t\to_sync%s\n"
	"\n", (dbg)?"_dbg":"", // resetw.c_str(), (dbg)?", o_dbg":"",
	TST_DBLSTAGE_IWIDTH, TST_DBLSTAGE_SHIFT,
		resetw.c_str(), (dbg) ? ",":"");

	if (dbg)
		fprintf(fp, "\toutput\twire\t[33:0]\t\t\to_dbg;\n");

	fprintf(fp, "\t\t// }}}\n\t);\n\n");

	fprintf(fp, "\t// Local declarations\n\t// {{{\n");
	fprintf(fp,
	"\twire\tsigned\t[(IWIDTH-1):0]\ti_in_0r, i_in_0i, i_in_1r, i_in_1i;\n"
	"\twire\t[(OWIDTH-1):0]\t\to_out_0r, o_out_0i,\n"
				"\t\t\t\t\to_out_1r, o_out_1i;\n"
	"\treg\trnd_sync, r_sync;\n"
	"\treg\tsigned\t[(IWIDTH):0]\trnd_in_0r, rnd_in_0i;\n"
	"\treg\tsigned\t[(IWIDTH):0]\trnd_in_1r, rnd_in_1i;\n\n"
	"\t// }}}\n"
"\n");

	fprintf(fp,
	"\tassign\ti_in_0r = i_left[(2*IWIDTH-1):(IWIDTH)];\n"
	"\tassign\ti_in_0i = i_left[(IWIDTH-1):0];\n"
	"\tassign\ti_in_1r = i_right[(2*IWIDTH-1):(IWIDTH)];\n"
	"\tassign\ti_in_1i = i_right[(IWIDTH-1):0];\n");
	if (dbg) fprintf(fp,
		"\tassign\to_dbg = { ((o_sync)&&(i_ce)), i_ce, o_left[(2*OWIDTH-1):(2*OWIDTH-16)],\n"
			"\t\t\t\t\to_left[(OWIDTH-1):(OWIDTH-16)] };\n"
		"\n");

	fprintf(fp,
	"\n"
	"\t// rnd_sync, r_sync\n"
	"\t// {{{\n"
	"\t// As with any register connected to the sync pulse, these must\n"
	"\t// have initial values and be reset on the %s signal.\n"
	"\t// Other data values need only restrict their updates to i_ce\n"
	"\t// enabled clocks, but sync\'s must obey resets and initial\n"
	"\t// conditions as well.\n"
"\n"
	"\tinitial\trnd_sync      = 1\'b0; // Sync into rounding\n"
	"\tinitial\tr_sync        = 1\'b0; // Sync coming out\n",
		resetw.c_str());
	if (async_reset)
		fprintf(fp, "\talways @(posedge i_clk, negdge i_areset_n)\n\tif (!i_areset_n)\n");
	else
		fprintf(fp, "\talways @(posedge i_clk)\n\tif (i_reset)\n");
	fprintf(fp,
		"\t\tbegin\n"
			"\t\t\trnd_sync <= 1\'b0;\n"
			"\t\t\tr_sync <= 1\'b0;\n"
		"\t\tend else if (i_ce)\n"
		"\t\tbegin\n"
			"\t\t\trnd_sync <= i_sync;\n"
			"\t\t\tr_sync <= rnd_sync;\n"
		"\t\tend\n"
	"\t// }}}\n"
"\n"
	"\t// rnd_in_0r, rnd_in_0i, rnd_in_1r, rnd_in_1i\n"
	"\t// {{{\n"
	"\t// As with other variables, these are really only updated when in\n"
	"\t// the processing pipeline, after the first i_sync.  However, to\n"
	"\t// eliminate as much unnecessary logic as possible, we toggle\n"
	"\t// these any time the i_ce line is enabled, and don\'t reset.\n"
	"\t// them on %s.\n", resetw.c_str());
	fprintf(fp,
	"\t// Don't forget that we accumulate a bit by adding two values\n"
	"\t// together. Therefore our intermediate value must have one more\n"
	"\t// bit than the two originals.\n"
	"\talways @(posedge i_clk)\n"
	"\tif (i_ce)\n"
	"\tbegin\n"
		"\t\t//\n"
		"\t\trnd_in_0r <= i_in_0r + i_in_1r;\n"
		"\t\trnd_in_0i <= i_in_0i + i_in_1i;\n"
		"\t\t//\n"
		"\t\trnd_in_1r <= i_in_0r - i_in_1r;\n"
		"\t\trnd_in_1i <= i_in_0i - i_in_1i;\n"
		"\t\t//\n"
	"\tend\n"
	"\t// }}}\n"
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
	"\t// o_left, o_right\n"
	"\t// {{{\n"
	"\t// Prior versions of this routine did not include the extra\n"
	"\t// clock and register/flip-flops that this routine requires.\n"
	"\t// These are placed in here to correct a bug in Verilator, that\n"
	"\t// otherwise struggles.  (Hopefully this will fix the problem ...)\n"
	"\talways @(posedge i_clk)\n"
	"\tif (i_ce)\n"
	"\tbegin\n"
		"\t\to_left  <= { o_out_0r, o_out_0i };\n"
		"\t\to_right <= { o_out_1r, o_out_1i };\n"
	"\tend\n"
	"\t// }}}\n"
"\n"
	"\t// o_sync\n"
	"\t// {{{\n"
	"\tinitial\to_sync = 1\'b0; // Final sync coming out of module\n");
	if (async_reset)
		fprintf(fp, "\talways @(posedge i_clk, negedge i_areset_n)\n\t\tif (!i_areset_n)\n");
	else
		fprintf(fp, "\talways @(posedge i_clk)\n\tif (i_reset)\n");
	fprintf(fp,
		"\t\to_sync <= 1\'b0;\n"
		"\telse if (i_ce)\n"
		"\t\to_sync <= r_sync;\n"
	"\t// }}}\n"
"\n"
"endmodule\n");
	fclose(fp);
}
// }}}

// build_stage
// {{{
void	build_stage(const char *fname,
		int stage, int nwide, int offset,
		int nbits, int xtra, int ckpce,
		const bool async_reset, const bool dbg) {
	FILE	*fstage = fopen(fname, "w");
	// int	cbits = nbits + xtra;

	std::string	resetw("i_reset");
	if (async_reset)
		resetw = std::string("i_areset_n");

	if (fstage == NULL) {
		fprintf(stderr, "ERROR: Could not open %s for writing!\n", fname);
		perror("O/S Err was:");
		fprintf(stderr, "Attempting to continue, but this file will be missing.\n");
		return;
	}

	fprintf(fstage,
SLASHLINE
"//\n"
"// Filename:\tfftstage%s.v\n"
"// {{{\n" // "}}}"
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
		(dbg)?"_dbg":"", prjname, creator);
	fprintf(fstage, "%s", cpyleft);
	fprintf(fstage, "//\n//\n`default_nettype\tnone\n//\n");
	fprintf(fstage, "module\tfftstage%s #(\n", (dbg)?"_dbg":"");
	fprintf(fstage, "\t\t// {{{\n");
	// These parameter values are useless at this point--they are to be
	// replaced by the parameter values in the calling program.  Only
	// problem is, the CWIDTH needs to match exactly!
	fprintf(fstage, "\t\tparameter\tIWIDTH=%d,CWIDTH=%d,OWIDTH=%d,\n",
		nbits, 20, nbits+1); // 20, not cbits, since the tb depends upon it
	fprintf(fstage,
"\t\t// Parameters specific to the core that should be changed when\n"
"\t\t// this core is built ... Note that the minimum LGSPAN (the base\n"
"\t\t// two log of the span, or the base two log of the current FFT\n"
"\t\t// size) is 3.  Smaller spans (i.e. the span of 2) must use the\n"
"\t\t// dbl laststage module.\n"
"\t\t// Verilator lint_off UNUSED\n"
"\t\tparameter\tLGSPAN=%d, BFLYSHIFT=0, // LGWIDTH=%d\n"
"\t\tparameter [0:0]\tOPT_HWMPY = 1,\n",
		(nwide <= 1) ? lgval(stage)-1 : lgval(stage)-2, lgval(stage));
	fprintf(fstage,
"\t\t// Clocks per CE.  If your incoming data rate is less than 50%%\n"
"\t\t// of your clock speed, you can set CKPCE to 2\'b10, make sure\n"
"\t\t// there\'s at least one clock between cycles when i_ce is high,\n"
"\t\t// and then use two multiplies instead of three.  Setting CKPCE\n"
"\t\t// to 2\'b11, and insisting on at least two clocks with i_ce low\n"
"\t\t// between cycles with i_ce high, then the hardware optimized\n"
"\t\t// butterfly code will used one multiply instead of two.\n"
"\t\tparameter\tCKPCE = %d,\n", ckpce);

	fprintf(fstage,
"\t\t// The COEFFILE parameter contains the name of the file\n"
"\t\t// containing the FFT twiddle factors\n");
	if (nwide == 2) {
		fprintf(fstage, "\t\tparameter\tCOEFFILE=\"cmem_%c%d.hex\",\n",
			(offset)?'o':'e', stage*2);
	} else
		fprintf(fstage,
			"\t\tparameter\tCOEFFILE=\"cmem_%d.hex\",\n",
			stage);
	fprintf(fstage, "\t\t// Verilator lint_on  UNUSED\n");

	fprintf(fstage,"\n"
"`ifdef	VERILATOR\n"
	"\t\tparameter  [0:0]\tZERO_ON_IDLE = 1'b0\n"
"`else\n"
	"\t\tlocalparam [0:0]\tZERO_ON_IDLE = 1'b0\n"
"`endif // VERILATOR\n"
	"\t\t// }}}\n\t) (\n\t\t// {{{\n");

	fprintf(fstage,
	"\t\tinput\twire\t			i_clk, %s,\n"
			"\t\t\t\t\t\t\ti_ce, i_sync,\n"
	"\t\tinput\twire\t[(2*IWIDTH-1):0]	i_data,\n"
	"\t\toutput\treg\t[(2*OWIDTH-1):0]	o_data,\n"
	"\t\toutput\treg\t			o_sync%s\n"
"\n", resetw.c_str(), (dbg) ? ",":"");
	if (dbg) { fprintf(fstage, "\t\toutput\twire\t[33:0]\t\t\to_dbg\n");
	}
	fprintf(fstage, "\t\t// }}}\n\t);\n\n");

	fprintf(fstage,
	"\t// Local signal definitions\n"
	"\t// {{{\n"
	"\t// I am using the prefixes\n"
	"\t// 	ib_*	to reference the inputs to the butterfly, and\n"
	"\t// 	ob_*	to reference the outputs from the butterfly\n"
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
"\treg	[(2*CWIDTH-1):0]	cmem [0:((1<<LGSPAN)-1)];\n");

	if (formal_property_flag)
		fprintf(fstage, 
			"`ifdef	FORMAL\n"
			"// Let the formal tool pick the coefficients\n"
			"`else\n");
	fprintf(fstage, "\tinitial\t$readmemh(COEFFILE,cmem);\n\n");
	if (formal_property_flag)
		fprintf(fstage, "`endif\n\n");

	// gen_coeff_file(coredir, fname, stage, cbits, nwide, offset, inv);

	fprintf(fstage,
	"\treg\t[(LGSPAN):0]		iaddr;\n"
	"\treg\t[(2*IWIDTH-1):0]	imem	[0:((1<<LGSPAN)-1)];\n"
"\n"
	"\treg\t[LGSPAN:0]		oaddr;\n"
	"\treg\t[(2*OWIDTH-1):0]	omem	[0:((1<<LGSPAN)-1)];\n"
"\n"
	"\twire\t\t\t\tidle;\n"
	"\treg	[(LGSPAN-1):0]\t\tnxt_oaddr;\n"
	"\treg	[(2*OWIDTH-1):0]\tpre_ovalue;\n"
	"\t// }}}\n"
"\n");
	if (dbg) fprintf(fstage, 
	"\tassign\to_dbg = { ((o_sync)&&(i_ce)), i_ce, o_data[(2*OWIDTH-1):(2*OWIDTH-16)],\n"
			"\t\t\t\t\to_data[(OWIDTH-1):(OWIDTH-16)] };\n"
"\n");

	fprintf(fstage, 
	"\t// wait_for_sync, iaddr\n"
	"\t// {{{\n"
	"\tinitial wait_for_sync = 1\'b1;\n"
	"\tinitial iaddr = 0;\n");
	if (async_reset)
		fprintf(fstage, "\talways @(posedge i_clk, negedge i_areset_n)\n\tif (!i_areset_n)\n");
	else
		fprintf(fstage, "\talways @(posedge i_clk)\n\tif (i_reset)\n");

	fprintf(fstage,
	"\tbegin\n"
		"\t\twait_for_sync <= 1\'b1;\n"
		"\t\tiaddr <= 0;\n"
	"\tend else if ((i_ce)&&((!wait_for_sync)||(i_sync)))\n"
	"\tbegin\n"
		"\t\t//\n"
		"\t\t// First step: Record what we\'re not ready to use yet\n"
		"\t\t//\n"
		"\t\tiaddr <= iaddr + { {(LGSPAN){1\'b0}}, 1\'b1 };\n"
		"\t\twait_for_sync <= 1\'b0;\n"
	"\tend\n"
	"\t// }}}\n"
"\n"
	"\t// Write to imem\n"
	"\t// {{{\n"
	"\talways @(posedge i_clk) // Need to make certain here that we don\'t read\n"
	"\tif ((i_ce)&&(!iaddr[LGSPAN])) // and write the same address on\n"
		"\t\timem[iaddr[(LGSPAN-1):0]] <= i_data; // the same clk\n"
	"\t// }}}\n"
	"\n");

	fprintf(fstage,
	"\t// ib_sync\n"
	"\t// {{{\n"
	"\t// Now, we have all the inputs, so let\'s feed the butterfly\n"
	"\t//\n"
	"\t// ib_sync is the synchronization bit to the butterfly.  It will\n"
	"\t// be tracked within the butterfly, and used to create the o_sync\n"
	"\t// value when the results from this output are produced\n"
	"\tinitial ib_sync = 1\'b0;\n");
	if (async_reset)
		fprintf(fstage, "\talways @(posedge i_clk, negedge i_areset_n)\n\tif (!i_areset_n)\n");
	else
		fprintf(fstage, "\talways @(posedge i_clk)\n\tif (i_reset)\n");
	fprintf(fstage,
			"\t\tib_sync <= 1\'b0;\n"
		"\telse if (i_ce)\n"
		"\tbegin\n"
			"\t\t// Set the sync to true on the very first\n"
			"\t\t// valid input in, and hence on the very\n"
			"\t\t// first valid data out per FFT.\n"
			"\t\tib_sync <= (iaddr==(1<<(LGSPAN)));\n"
		"\tend\n\t// }}}\n\n"
	"\t// ib_a, ib_b, ib_c\n"
	"\t// {{{\n"
	"\t// Read the values from our input memory, and use them to feed\n"
	"\t// first of two butterfly inputs\n"
	"\talways\t@(posedge i_clk)\n"
	"\tif (i_ce)\n"
	"\tbegin\n"
		"\t\t// One input from memory, ...\n"
		"\t\tib_a <= imem[iaddr[(LGSPAN-1):0]];\n"
		"\t\t// One input clocked in from the top\n"
		"\t\tib_b <= i_data;\n"
		"\t\t// and the coefficient or twiddle factor\n"
		"\t\tib_c <= cmem[iaddr[(LGSPAN-1):0]];\n"
	"\tend\n\t// }}}\n\n");

	fprintf(fstage,
	"\t// idle\n"
	"\t// {{{\n"
	"\t// The idle register is designed to keep track of when an input\n"
	"\t// to the butterfly is important and going to be used.  It's used\n"
	"\t// in a flag following, so that when useful values are placed\n"
	"\t// into the butterfly they'll be non-zero (idle=0), otherwise when\n"
	"\t// the inputs to the butterfly are irrelevant and will be ignored,\n"
	"\t// then (idle=1) those inputs will be set to zero.  This\n"
	"\t// functionality is not designed to be used in operation, but only\n"
	"\t// within a Verilator simulation context when chasing a bug.\n"
	"\t// In this limited environment, the non-zero answers will stand\n"
	"\t// in a trace making it easier to highlight a bug.\n"
	"\tgenerate if (ZERO_ON_IDLE)\n"
	"\tbegin : GEN_ZERO_ON_IDLE\n");

		fprintf(fstage,
		"\t\treg	r_idle;\n\n"
		"\t\tinitial	r_idle = 1;\n");
		if (async_reset)
			fprintf(fstage,
		"\t\talways @(posedge i_clk, negedge i_areset_n)\n"
		"\t\tif (!i_areset_n)\n");
		else
			fprintf(fstage,
		"\t\talways @(posedge i_clk)\n"
		"\t\tif (i_reset)\n");

		fprintf(fstage,
			"\t\t\tr_idle <= 1\'b1;\n"
		"\t\telse if (i_ce)\n"
			"\t\t\tr_idle <= (!iaddr[LGSPAN])&&(!wait_for_sync);\n\n"
		"\t\tassign\tidle = r_idle;\n\n");

	fprintf(fstage,
	"\tend else begin : NO_IDLE_GENERATION\n\n"
	"\t\tassign\tidle = 0;\n\n"
	"\tend endgenerate\n"
	"\t// }}}\n\n");

	fprintf(fstage,
	"\t////////////////////////////////////////////////////////////////////////\n"
	"\t//\n"
	"\t// Instantiate the butterfly\n"
	"\t// {{{\n"
	"\t////////////////////////////////////////////////////////////////////////\n"
	"\t//\n"
	"\t//\n");

	if (formal_property_flag)
		fprintf(fstage,
"// For the formal proof, we'll assume the outputs of hwbfly and/or\n"
"// butterfly, rather than actually calculating them.  This will simplify\n"
"// the proof and (if done properly) will be equivalent.  Be careful of\n"
"// defining FORMAL if you want the full logic!\n"
"`ifndef	FORMAL\n"
	"\t//\n");

	fprintf(fstage,
"\tgenerate if (OPT_HWMPY)\n"
"\tbegin : HWBFLY\n"
"\n"
	"\t\thwbfly #(\n"
		"\t\t\t// {{{\n"
		"\t\t\t.IWIDTH(IWIDTH),\n"
		"\t\t\t.CWIDTH(CWIDTH),\n"
		"\t\t\t.OWIDTH(OWIDTH),\n"
		"\t\t\t.CKPCE(CKPCE),\n"
		"\t\t\t.SHIFT(BFLYSHIFT)\n"
		"\t\t\t// }}}\n"
	"\t\t) bfly(\n"
		"\t\t\t// {{{\n");
		if (async_reset)
			fprintf(fstage,
		"\t\t\t.i_clk(i_clk), .i_areset_n(i_areset_n), .i_ce(i_ce),\n");
		else
			fprintf(fstage,
		"\t\t\t.i_clk(i_clk), .i_reset(i_reset), .i_ce(i_ce),\n");

		fprintf(fstage,
		"\t\t\t.i_coef((idle && !i_ce) ? 0:ib_c),\n"
		"\t\t\t.i_left((idle && !i_ce) ? 0:ib_a),\n"
		"\t\t\t.i_right((idle && !i_ce) ? 0:ib_b),\n"
		"\t\t\t.i_aux(ib_sync && i_ce),\n"
		"\t\t\t.o_left(ob_a), .o_right(ob_b), .o_aux(ob_sync)\n"
		"\t\t\t// }}}\n"
	"\t\t);\n"
"\n"
"\tend else begin : FWBFLY\n"
"\n"
	"\t\tbutterfly #(\n"
		"\t\t\t// {{{\n"
		"\t\t\t.IWIDTH(IWIDTH),\n"
		"\t\t\t.CWIDTH(CWIDTH),\n"
		"\t\t\t.OWIDTH(OWIDTH),\n"
		"\t\t\t.CKPCE(CKPCE),\n"
		"\t\t\t.SHIFT(BFLYSHIFT)\n"
		"\t\t\t// }}}\n"
	"\t\t) bfly(\n"
		"\t\t\t// {{{\n");

		if (async_reset)
			fprintf(fstage,
		"\t\t\t.i_clk(i_clk), .i_areset_n(i_areset_n), .i_ce(i_ce),\n");
		else
			fprintf(fstage,
		"\t\t\t.i_clk(i_clk), .i_reset(i_reset), .i_ce(i_ce),\n");

		fprintf(fstage,
		"\t\t\t.i_coef( (idle && !i_ce)?0:ib_c),\n"
		"\t\t\t.i_left( (idle && !i_ce)?0:ib_a),\n"
		"\t\t\t.i_right((idle && !i_ce)?0:ib_b),\n"
		"\t\t\t.i_aux(ib_sync && i_ce),\n"
		"\t\t\t.o_left(ob_a), .o_right(ob_b), .o_aux(ob_sync)\n"
		"\t\t\t// }}}\n"
	"\t\t);\n"
"\n"
"\tend endgenerate\n");

	if (formal_property_flag) {
		fprintf(fstage, "`else\n"
"\n"
	"\t// Verilator lint_off UNDRIVEN\n"
	"\t(* anyseq *)    wire    [(2*OWIDTH-1):0]        f_ob_a, f_ob_b;\n"
	"\t(* anyseq *)    wire    f_ob_sync;\n"
	"\t// Verilator lint_on  UNDRIVEN\n"
"\n"
	"\tassign  ob_sync = f_ob_sync;\n"
	"\tassign  ob_a    = f_ob_a;\n"
	"\tassign  ob_b    = f_ob_b;\n\n");

		fprintf(fstage, "`endif\n\n");
	}
	fprintf(fstage, "\t// }}}\n\n");

	fprintf(fstage,
	"\t// oaddr, o_sync, b_started\n"
	"\t// {{{\n"
	"\t// Next step: recover the outputs from the butterfly\n"
	"\t//\n"
	"\t// The first output can go immediately to the output of this routine\n"
	"\t// The second output must wait until this time in the idle cycle\n"
	"\t// oaddr is the output memory address, keeping track of where we are\n"
	"\t// in this output cycle.\n"
	"\tinitial oaddr     = 0;\n"
	"\tinitial o_sync    = 0;\n"
	"\tinitial b_started = 0;\n");
	if (async_reset)
		fprintf(fstage, "\talways @(posedge i_clk, negedge i_areset_n)\n\tif (!i_areset_n)\n");
	else
		fprintf(fstage, "\talways @(posedge i_clk)\n\tif (i_reset)\n");
	fprintf(fstage,
	"\tbegin\n"
		"\t\toaddr     <= 0;\n"
		"\t\to_sync    <= 0;\n"
		"\t\t// b_started will be true once we've seen the first ob_sync\n"
		"\t\tb_started <= 0;\n"
	"\tend else if (i_ce)\n"
	"\tbegin\n"
	"\t\to_sync <= (!oaddr[LGSPAN])?ob_sync : 1\'b0;\n"
	"\t\tif (ob_sync||b_started)\n"
		"\t\t\toaddr <= oaddr + 1\'b1;\n"
	"\t\tif ((ob_sync)&&(!oaddr[LGSPAN]))\n"
		"\t\t\t// If b_started is true, then a butterfly output\n"
		"\t\t\t// is available\n"
		"\t\t\tb_started <= 1\'b1;\n"
	"\tend\n\t// }}}\n\n");

	fprintf(fstage,
	"\t// nxt_oaddr\n"
	"\t// {{{\n"
	"\talways @(posedge i_clk)\n"
	"\tif (i_ce)\n"
		"\t\tnxt_oaddr[0] <= oaddr[0];\n"
	"\tgenerate if (LGSPAN>1)\n"
	"\tbegin\n"
"\n"
	"\t\talways @(posedge i_clk)\n"
	"\t\tif (i_ce)\n"
		"\t\t\tnxt_oaddr[LGSPAN-1:1] <= oaddr[LGSPAN-1:1] + 1\'b1;\n"
"\n"
	"\tend endgenerate\n"
	"\t// }}}\n"
"\n"
	"\t// omem\n"
	"\t// {{{\n"
	"\t// Only write to the memory on the first half of the outputs\n"
	"\t// We'll use the memory value on the second half of the outputs\n"
	"\talways @(posedge i_clk)\n"
	"\tif ((i_ce)&&(!oaddr[LGSPAN]))\n"
		"\t\tomem[oaddr[(LGSPAN-1):0]] <= ob_b;\n\t// }}}\n\n");

	fprintf(fstage,
	"\t// pre_ovalue\n"
	"\t// {{{\n"
	"\talways @(posedge i_clk)\n"
	"\tif (i_ce)\n"
		"\t\tpre_ovalue <= omem[nxt_oaddr[(LGSPAN-1):0]];\n"
	"\t// }}}\n"
"\n");
	fprintf(fstage,
	"\t// o_data\n"
	"\t// {{{\n"
	"\talways @(posedge i_clk)\n"
	"\tif (i_ce)\n"
	"\t\to_data <= (!oaddr[LGSPAN]) ? ob_a : pre_ovalue;\n"
	"\t// }}}\n"
"\n");

	fprintf(fstage,
"////////////////////////////////////////////////////////////////////////////////\n"
"////////////////////////////////////////////////////////////////////////////////\n"
"////////////////////////////////////////////////////////////////////////////////\n"
"//\n"
"// Formal properties\n"
"// {{{\n"
"////////////////////////////////////////////////////////////////////////////////\n"
"////////////////////////////////////////////////////////////////////////////////\n"
"////////////////////////////////////////////////////////////////////////////////\n"
"`ifdef	FORMAL\n");


	if (formal_property_flag) {

	fprintf(fstage,
	"\t// Local (formal) declarations\n"
	"\t// {{{\n"
	"\t// An arbitrary processing delay from butterfly input to\n"
	"\t// butterfly output(s)\n"
	"\t// Verilator lint_off UNDRIVEN\n"
	"\t(* anyconst *) reg	[LGSPAN:0]	f_mpydelay;\n"
	"\t(* anyconst *)	reg	[LGSPAN:0]	f_addr;\n"
	"\t// Verilator lint_on  UNDRIVEN\n"
	"\treg	[2*IWIDTH-1:0]			f_left, f_right;\n"
	"\treg	[2*OWIDTH-1:0]	f_oleft, f_oright;\n"
	"\treg	[LGSPAN:0]	f_oaddr;\n"
	"\twire	[LGSPAN:0]	f_oaddr_m1 = f_oaddr - 1'b1;\n"
	"\treg	f_output_active;\n"
	"\t// }}}\n"
"\n"
"\n"
	"\talways @(*)\n"
	"\t\tassume(f_mpydelay > 1);\n"
"\n"
	"\treg	f_past_valid;\n"
	"\tinitial	f_past_valid = 1'b0;\n"
	"\talways @(posedge i_clk)\n"
		"\t\tf_past_valid <= 1'b1;\n"
"\n");

	if (async_reset)
		fprintf(fstage, "\talways @(*)\n\tif ((!f_past_valid)||(!i_areset_n))\n");
	else
		fprintf(fstage, "\talways @(posedge i_clk)\n"
				"\tif ((!f_past_valid)||($past(i_reset)))\n");
	fprintf(fstage,
	"\tbegin\n"
		"\t\tassert(iaddr == 0);\n"
		"\t\tassert(wait_for_sync);\n"
		"\t\tassert(o_sync == 0);\n"
		"\t\tassert(oaddr == 0);\n"
		"\t\tassert(!b_started);\n"
		"\t\tassert(!o_sync);\n"
	"\tend\n\n");

	fprintf(fstage,
	"\t////////////////////////////////////////////////////////////////////////\n"
	"\t//\n"
	"\t// Formally verify the input half, from the inputs to this module\n"
	"\t// to the inputs of the butterfly\n"
	"\t//\n"
	"\t////////////////////////////////////////////////////////////////////////\n"
	"\t//\n\t//\n\n"
	"\t// Let's  verify a specific set of inputs\n"
"\n"
	"\talways @(posedge i_clk)\n"
	"\tif (!$past(i_ce) && !$past(i_ce,2) && !$past(i_ce,3) && !$past(i_ce,4))\n"
		"\t\tassume(!i_ce);\n"
"\n"
	"\talways @(*)\n"
		"\t\tassume(f_addr[LGSPAN]==1'b0);\n"
"\n"
	"\talways @(posedge i_clk)\n"
	"\tif ((i_ce)&&(iaddr[LGSPAN:0] == f_addr))\n"
		"\t\tf_left <= i_data;\n"
"\n"
	"\talways @(*)\n"
	"\tif (wait_for_sync)\n"
		"\t\tassert(iaddr == 0);\n"
"\n"
	"\twire	[LGSPAN:0]\tf_last_addr = iaddr - 1'b1;\n"
"\n"
	"\talways @(posedge i_clk)\n"
	"\tif ((!wait_for_sync)&&(f_last_addr >= { 1'b0, f_addr[LGSPAN-1:0]}))\n"
		"\t\tassert(f_left == imem[f_addr[LGSPAN-1:0]]);\n"
"\n"
	"\talways @(posedge i_clk)\n"
	"\tif ((i_ce)&&(iaddr == { 1'b1, f_addr[LGSPAN-1:0]}))\n"
		"\t\tf_right <= i_data;\n"
"\n"
	"\talways @(posedge i_clk)\n"
	"\tif (i_ce && !wait_for_sync\n"
		"\t\t&& (f_last_addr == { 1'b1, f_addr[LGSPAN-1:0]}))\n"
	"\tbegin\n"
		"\t\tassert(ib_a == f_left);\n"
		"\t\tassert(ib_b == f_right);\n"
		"\t\tassert(ib_c == cmem[f_addr[LGSPAN-1:0]]);\n"
	"\tend\n\n");

	fprintf(fstage,
	"\t////////////////////////////////////////////////////////////////////////\n"
	"\t//\n"
	"\t// Formally verify the output half, from the output of the butterfly\n"
	"\t// to the outputs of this module\n"
	"\t//\n"
	"\t////////////////////////////////////////////////////////////////////////\n"
	"\t//\n\t//\n\n");

	fprintf(fstage,
	"\talways @(*)\n"
		"\t\tf_oaddr = iaddr - f_mpydelay + {1'b1,{(LGSPAN-1){1'b0}} };\n"
"\n"
	"\tassign\tf_oaddr_m1 = f_oaddr - 1'b1;\n"
"\n");

	fprintf(fstage,
	"\tinitial\tf_output_active = 1'b0;\n");
	if (async_reset)
		fprintf(fstage, "\talways @(posedge i_clk, negedge i_areset_n)\n\tif (!i_areset_n)\n");
	else
		fprintf(fstage, "\talways @(posedge i_clk)\n\tif (i_reset)\n");
	fprintf(fstage,
		"\t\tf_output_active <= 1'b0;\n"
	"\telse if ((i_ce)&&(ob_sync))\n"
		"\t\tf_output_active <= 1'b1;\n"
"\n"
	"\talways @(*)\n"
		"\t\tassert(f_output_active == b_started);\n"
"\n"
	"\talways @(*)\n"
	"\tif (wait_for_sync)\n"
		"\t\tassert(!f_output_active);\n\n");

	fprintf(fstage,
	"\talways @(*)\n"
	"\tif (f_output_active)\n"
	"\tbegin\n"
		"\t\tassert(oaddr == f_oaddr);\n"
	"\tend else\n"
		"\t\tassert(oaddr == 0);\n"
"\n");


	fprintf(fstage,
	"\talways @(*)\n"
	"\tif (wait_for_sync)\n"
		"\t\tassume(!ob_sync);\n"
"\n"
	"\talways @(*)\n"
		"\t\tassume(ob_sync == (f_oaddr == 0));\n"
"\n"
	"\talways @(posedge i_clk)\n"
	"\tif ((f_past_valid)&&(!$past(i_ce)))\n"
	"\tbegin\n"
		"\t\tassume($stable(ob_a));\n"
		"\t\tassume($stable(ob_b));\n"
	"\tend\n\n");

	fprintf(fstage,
	"\tinitial	f_oleft  = 0;\n"
	"\tinitial	f_oright = 0;\n"
	"\talways @(posedge i_clk)\n"
	"\tif ((i_ce)&&(f_oaddr == f_addr))\n"
	"\tbegin\n"
		"\t\tf_oleft  <= ob_a;\n"
		"\t\tf_oright <= ob_b;\n"
	"\tend\n"
"\n"
	"\talways @(posedge i_clk)\n"
	"\tif ((f_output_active)&&(f_oaddr_m1 >= { 1'b0, f_addr[LGSPAN-1:0]}))\n"
		"\t\tassert(omem[f_addr[LGSPAN-1:0]] == f_oright);\n"
"\n"
	"\talways @(posedge i_clk)\n"
	"\tif ((i_ce)&&(f_oaddr_m1 == 0)&&(f_output_active))\n"
	"\tbegin\n"
		"\t\tassert(o_sync);\n"
	"\tend else if ((i_ce)||(!f_output_active))\n"
		"\t\tassert(!o_sync);\n"
	"\n"
	"\talways @(posedge i_clk)\n"
	"\tif ((i_ce)&&(f_output_active)&&(f_oaddr_m1 == f_addr))\n"
		"\t\tassert(o_data == f_oleft);\n"
"\n"
	"\talways @(posedge i_clk)\n"
	"\tif ((i_ce)&&(f_output_active)&&(f_oaddr[LGSPAN])\n"
			"\t\t\t&&(f_oaddr[LGSPAN-1:0] == f_addr[LGSPAN-1:0]))\n"
		"\t\tassert(pre_ovalue == f_oright);\n"
"\n"
	"\talways @(posedge i_clk)\n"
	"\tif ((i_ce)&&(f_output_active)&&(f_oaddr_m1[LGSPAN])\n"
			"\t\t\t&&(f_oaddr_m1[LGSPAN-1:0] == f_addr[LGSPAN-1:0]))\n"
		"\t\tassert(o_data == f_oright);\n"
"\n");

	fprintf(fstage,
	"\t// Make Verilator happy\n"
	"\t// {{{\n"
	"\t// Verilator lint_off UNUSED\n"
	"\twire	unused_formal;\n"
	"\tassign unused_formal = &{ 1\'b0, idle, ib_sync };\n"
	"\t// Verilator lint_on  UNUSED\n"
	"\t// }}}\n\n");

	} else { // If no formal properties
		fprintf(fstage, "// Formal properties exist, but are not enabled"
				" in this build\n");
	} // End of the formal properties section

	fprintf(fstage,
"`endif // FORMAL\n"
"// }}}\n");

	fprintf(fstage, "endmodule\n");
}
// }}}
