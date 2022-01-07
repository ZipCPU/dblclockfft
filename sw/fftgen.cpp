////////////////////////////////////////////////////////////////////////////////
//
// Filename: 	fftgen.cpp
// {{{
// Project:	A General Purpose Pipelined FFT Implementation
//
// Purpose:	This is the core generator for the project.  Every part
//		and piece of this project begins and ends in this program.
//	Once built, this program will build an FFT (or IFFT) core of arbitrary
//	width, precision, etc., that will run at two samples per clock.
//	(Incidentally, I didn't pick two samples per clock because it was
//	easier, but rather because there weren't any two-sample per clock
//	FFT's posted on opencores.com.  Further, FFT's running at one sample
//	per aren't that hard to find.)
//
//	You can find the documentation for this program in two places.  One is
//	in the usage() function below.  The second is in the 'doc'uments
//	directory that comes with this package, specifically in the spec.pdf
//	file.  If it's not there, type make in the documents directory to
//	build it.
//
//	20160123 - Thanks to Lesha Birukov, adjusted for MS Visual Studio 2012.
//		(Adjustments are at the top of the file ...)
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
// {{{
#include <io.h>
#include <direct.h>
#define _USE_MATH_DEFINES
#define	R_OK    4       /* Test for read permission.  */
#define	W_OK    2       /* Test for write permission.  */
#define	X_OK    0       /* !!!!!! execute permission - unsupported in windows*/
#define	F_OK    0       /* Test for existence.  */

#if _MSC_VER <= 1700

int lstat(const char *filename, struct stat *buf) { return 1; };
#define	S_ISDIR(A)	0

#else

#define	lstat	_stat
#define S_ISDIR	_S_IFDIR

#endif

#define	mkdir(A,B)	_mkdir(A)

#define access _access
// }}}
#elif defined(__MINGW32__)
// MinGW - Minimalist GNU for Windows -- also MinGW-w64/MSYS2
// {{{
#include <direct.h>	// mkdir
#include <unistd.h>
#include <sys/stat.h>

#define	mkdir(A,B)	_mkdir(A)
#define	lstat(p,s)	stat(p,s)
// }}}
#else
// And for G++/Linux environment
// {{{
#include <unistd.h>	// Defines the R_OK/W_OK/etc. macros
#include <sys/stat.h>
// }}}
#endif

#include <string.h>
#include <string>
#include <math.h>
#include <ctype.h>
#include <assert.h>

#include "defaults.h"
#include "legal.h"
#include "rounding.h"
#include "fftlib.h"
#include "bldstage.h"
#include "bitreverse.h"
#include "softmpy.h"
#include "butterfly.h"

// build_dblquarters
// {{{
void	build_dblquarters(const char *fname, ROUND_T rounding, const bool async_reset=false, const bool dbg=false) {
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


	fprintf(fp,
SLASHLINE
"//\n"
"// Filename:\tqtrstage%s.v\n"
"// {{{\n" // "}}}"
"// Project:\t%s\n"
"//\n"
"// Purpose:	This file encapsulates the 4 point stage of a decimation in\n"
"//		frequency FFT.  This particular implementation is optimized\n"
"//	so that all of the multiplies are accomplished by additions and\n"
"//	multiplexers only.\n"
"//\n"
"//\n%s"
"//\n",
		(dbg)?"_dbg":"", prjname, creator);
	fprintf(fp, "%s", cpyleft);
	fprintf(fp, "//\n//\n`default_nettype\tnone\n//\n");

	std::string	resetw("i_reset");
	if (async_reset)
		resetw = std::string("i_areset_n");

	fprintf(fp,
"module\tqtrstage%s(i_clk, %s, i_ce, i_sync, i_data, o_data, o_sync%s);\n"
	"\tparameter	IWIDTH=%d, OWIDTH=IWIDTH+1;\n"
	"\t// Parameters specific to the core that should be changed when this\n"
	"\t// core is built ... Note that the minimum LGSPAN is 2.  Smaller\n"
	"\t// spans must use the fftdoubles stage.\n"
	"\tparameter\tLGWIDTH=%d, ODD=0, INVERSE=0,SHIFT=0;\n"
	"\tinput\twire				i_clk, %s, i_ce, i_sync;\n"
	"\tinput\twire	[(2*IWIDTH-1):0]	i_data;\n"
	"\toutput\treg	[(2*OWIDTH-1):0]	o_data;\n"
	"\toutput\treg				o_sync;\n"
	"\n", (dbg)?"_dbg":"",
	resetw.c_str(),
	(dbg)?", o_dbg":"", TST_QTRSTAGE_IWIDTH,
	TST_QTRSTAGE_LGWIDTH, resetw.c_str());
	if (dbg) { fprintf(fp, "\toutput\twire\t[33:0]\t\t\to_dbg;\n"
		"\tassign\to_dbg = { ((o_sync)&&(i_ce)), i_ce, o_data[(2*OWIDTH-1):(2*OWIDTH-16)],\n"
			"\t\t\t\t\to_data[(OWIDTH-1):(OWIDTH-16)] };\n"
"\n");
	}
	fprintf(fp,
	"\treg\t	wait_for_sync;\n"
	"\treg\t[3:0]	pipeline;\n"
"\n"
	"\treg\t[(IWIDTH):0]	sum_r, sum_i, diff_r, diff_i;\n"
"\n"
	"\treg\t[(2*OWIDTH-1):0]\tob_a;\n"
	"\twire\t[(2*OWIDTH-1):0]\tob_b;\n"
	"\treg\t[(OWIDTH-1):0]\t\tob_b_r, ob_b_i;\n"
	"\tassign\tob_b = { ob_b_r, ob_b_i };\n"
"\n"
	"\treg\t[(LGWIDTH-1):0]\t\tiaddr;\n"
	"\treg\t[(2*IWIDTH-1):0]\timem;\n"
"\n"
	"\twire\tsigned\t[(IWIDTH-1):0]\timem_r, imem_i;\n"
	"\tassign\timem_r = imem[(2*IWIDTH-1):(IWIDTH)];\n"
	"\tassign\timem_i = imem[(IWIDTH-1):0];\n"
"\n"
	"\twire\tsigned\t[(IWIDTH-1):0]\ti_data_r, i_data_i;\n"
	"\tassign\ti_data_r = i_data[(2*IWIDTH-1):(IWIDTH)];\n"
	"\tassign\ti_data_i = i_data[(IWIDTH-1):0];\n"
"\n"
	"\treg	[(2*OWIDTH-1):0]	omem;\n"
"\n");
	fprintf(fp,
	"\twire\tsigned\t[(OWIDTH-1):0]\trnd_sum_r, rnd_sum_i, rnd_diff_r, rnd_diff_i,\n");
	fprintf(fp,
	"\t\t\t\t\tn_rnd_diff_r, n_rnd_diff_i;\n");
	fprintf(fp,
	"\t%s #(IWIDTH+1,OWIDTH,SHIFT)\n\tdo_rnd_sum_r(i_clk, i_ce, "
	"sum_r, rnd_sum_r);\n\n", rnd_string);
	fprintf(fp,
	"\t%s #(IWIDTH+1,OWIDTH,SHIFT)\n\tdo_rnd_sum_i(i_clk, i_ce, "
	"sum_i, rnd_sum_i);\n\n", rnd_string);
	fprintf(fp,
	"\t%s #(IWIDTH+1,OWIDTH,SHIFT)\n\tdo_rnd_diff_r(i_clk, i_ce, "
	"diff_r, rnd_diff_r);\n\n", rnd_string);
	fprintf(fp,
	"\t%s #(IWIDTH+1,OWIDTH,SHIFT)\n\tdo_rnd_diff_i(i_clk, i_ce, "
	"diff_i, rnd_diff_i);\n\n", rnd_string);
	fprintf(fp, "\tassign n_rnd_diff_r = - rnd_diff_r;\n"
		"\tassign n_rnd_diff_i = - rnd_diff_i;\n");

	fprintf(fp,
	"\tinitial wait_for_sync = 1\'b1;\n"
	"\tinitial iaddr = 0;\n");
	if (async_reset)
		fprintf(fp,
			"\talways @(posedge i_clk, negedge i_areset_n)\n"
			"\tif (!i_areset_n)\n");
	else
		fprintf(fp,
	"\talways @(posedge i_clk)\n"
	"\tif (i_reset)\n");
	fprintf(fp,
	"\tbegin\n"
		"\t\twait_for_sync <= 1\'b1;\n"
		"\t\tiaddr <= 0;\n"
	"\tend else if ((i_ce)&&((!wait_for_sync)||(i_sync)))\n"
	"\tbegin\n"
		"\t\tiaddr <= iaddr + { {(LGWIDTH-1){1\'b0}}, 1\'b1 };\n"
		"\t\twait_for_sync <= 1\'b0;\n"
	"\tend\n\n"
	"\talways @(posedge i_clk)\n"
	"\tif (i_ce)\n"
		"\t\timem <= i_data;\n"
		"\n\n");
	fprintf(fp,
	"\t// Note that we don\'t check on wait_for_sync or i_sync here.\n"
	"\t// Why not?  Because iaddr will always be zero until after the\n"
	"\t// first i_ce, so we are safe.\n"
	"\tinitial pipeline = 4\'h0;\n");
	if (async_reset)
		fprintf(fp,
	"\talways\t@(posedge i_clk, negedge i_areset_n)\n"
	"\tif (!i_areset_n)\n");
	else
		fprintf(fp,
	"\talways\t@(posedge i_clk)\n"
	"\tif (i_reset)\n");

	fprintf(fp,
		"\t\tpipeline <= 4\'h0;\n"
	"\telse if (i_ce) // is our pipeline process full?  Which stages?\n"
		"\t\tpipeline <= { pipeline[2:0], iaddr[0] };\n\n");
	fprintf(fp,
	"\t// This is the pipeline[-1] stage, pipeline[0] will be set next.\n"
	"\talways\t@(posedge i_clk)\n"
	"\tif ((i_ce)&&(iaddr[0]))\n"
	"\tbegin\n"
		"\t\tsum_r  <= imem_r + i_data_r;\n"
		"\t\tsum_i  <= imem_i + i_data_i;\n"
		"\t\tdiff_r <= imem_r - i_data_r;\n"
		"\t\tdiff_i <= imem_i - i_data_i;\n"
	"\tend\n\n");
	fprintf(fp,
	"\t// pipeline[1] takes sum_x and diff_x and produces rnd_x\n\n");
	fprintf(fp,
	"\t// Now for pipeline[2].  We can actually do this at all i_ce\n"
	"\t// clock times, since nothing will listen unless pipeline[3]\n"
	"\t// on the next clock.  Thus, we simplify this logic and do\n"
	"\t// it independent of pipeline[2].\n"
	"\talways\t@(posedge i_clk)\n"
	"\tif (i_ce)\n"
	"\tbegin\n"
		"\t\tob_a <= { rnd_sum_r, rnd_sum_i };\n"
		"\t\t// on Even, W = e^{-j2pi 1/4 0} = 1\n"
		"\t\tif (ODD == 0)\n"
		"\t\tbegin\n"
		"\t\t\tob_b_r <= rnd_diff_r;\n"
		"\t\t\tob_b_i <= rnd_diff_i;\n"
		"\t\tend else if (INVERSE==0) begin\n"
		"\t\t\t// on Odd, W = e^{-j2pi 1/4} = -j\n"
		"\t\t\tob_b_r <=   rnd_diff_i;\n"
		"\t\t\tob_b_i <= n_rnd_diff_r;\n"
		"\t\tend else begin\n"
		"\t\t\t// on Odd, W = e^{j2pi 1/4} = j\n"
		"\t\t\tob_b_r <= n_rnd_diff_i;\n"
		"\t\t\tob_b_i <=   rnd_diff_r;\n"
		"\t\tend\n"
	"\tend\n\n");
	fprintf(fp,
	"\talways\t@(posedge i_clk)\n"
	"\tif (i_ce)\n"
	"\tbegin // In sequence, clock = 3\n"
		"\t\tif (pipeline[3])\n"
		"\t\tbegin\n"
			"\t\t\tomem <= ob_b;\n"
			"\t\t\to_data <= ob_a;\n"
		"\t\tend else\n"
			"\t\t\to_data <= omem;\n"
	"\tend\n\n");

	fprintf(fp,
	"\t// This algorithm takes five clocks to complete, therefore we can\n"
	"\t// set o_sync any time the address counter iaddr == 5.\n"
	"\t//\n"
	"\t// Don\'t forget in the sync check that we are running at two\n"
	"\t// clocks per sample.  Thus we need to produce a sync every\n"
	"\t// 2^(LGWIDTH-1) clocks.\n\t//\n"
	"\tinitial\to_sync = 1\'b0;\n");

	fprintf(fp,
	"\tgenerate if (LGWIDTH == 3)\n"
	"\tbegin\n\n");
		fprintf(fp,
		"\t\treg	o_sync_passed;\n\n"
		"\t\tinitial	o_sync_passed = 1\'b0;\n");

		if (async_reset)
			fprintf(fp,
		"\t\talways\t@(posedge i_clk, negedge i_areset_n)\n"
		"\t\tif (!i_areset_n)\n");
		else
			fprintf(fp,
		"\t\talways\t@(posedge i_clk)\n"
		"\t\tif (i_reset)\n");
		fprintf(fp,
			"\t\t\to_sync_passed <= 1\'b0;\n"
		"\t\telse if (i_ce && o_sync)\n"
			"\t\t\to_sync_passed <= 1\'b1;\n\n");

		if (async_reset)
			fprintf(fp,
		"\t\talways\t@(posedge i_clk, negedge i_areset_n)\n"
		"\t\tif (!i_areset_n)\n");
		else
			fprintf(fp,
		"\t\talways\t@(posedge i_clk)\n"
		"\t\tif (i_reset)\n");
		fprintf(fp,
			"\t\t\to_sync <= 1\'b0;\n"
		"\t\telse if (i_ce && (o_sync_passed || iaddr[2]))\n"
			"\t\t\to_sync <= (iaddr[1:0] == 2'b01);\n");

	fprintf(fp,
	"\n\tend else if (LGWIDTH == 4)\n"
	"\tbegin\n\n");

		if (async_reset)
			fprintf(fp,
		"\t\talways\t@(posedge i_clk, negedge i_areset_n)\n"
		"\t\tif (!i_areset_n)\n");
		else
			fprintf(fp,
		"\t\talways\t@(posedge i_clk)\n"
		"\t\tif (i_reset)\n");
		fprintf(fp,
			"\t\t\to_sync <= 1\'b0;\n"
		"\t\telse if (i_ce)\n"
			"\t\t\to_sync <= (iaddr[2:0] == 3'b101);\n");

	fprintf(fp,
	"\n\tend else begin\n\n");

		if (async_reset)
			fprintf(fp,
		"\t\talways\t@(posedge i_clk, negedge i_areset_n)\n"
		"\t\tif (!i_areset_n)\n");
		else
			fprintf(fp,
		"\t\talways\t@(posedge i_clk)\n"
		"\t\tif (i_reset)\n");
		fprintf(fp,
			"\t\t\to_sync <= 1\'b0;\n"
		"\t\telse if (i_ce)\n"
			"\t\t\t// As currently formulated, this line requires a\n"
			"\t\t\t// transform of 32 points or greater.  Notice\n"
			"\t\t\t// that the top bit is ignored, on purpose.\n"
			"\t\t\to_sync <= (iaddr[(LGWIDTH-2):3] == 0) && (iaddr[2:0] == 3'b101);\n");

	fprintf(fp, "\n\tend endgenerate\n\n");

	fprintf(fp, "endmodule\n");
}
// }}}

// build_snglquarters
// {{{
void	build_snglquarters(const char *fname, ROUND_T rounding, const bool async_reset=false, const bool dbg=false) {
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


	fprintf(fp,
SLASHLINE
"//\n"
"// Filename:\tqtrstage%s.v\n"
"// {{{\n" // "}}}"
"// Project:\t%s\n"
"//\n"
"// Purpose:	This file encapsulates the 4 point stage of a decimation in\n"
"//		frequency FFT.  This particular implementation is optimized\n"
"//	so that all of the multiplies are accomplished by additions and\n"
"//	multiplexers only.\n"
"//\n"
"// Operation:\n"
"// 	The operation of this stage is identical to the regular stages of\n"
"// 	the FFT (see them for details), with one additional and critical\n"
"// 	difference: this stage doesn't require any hardware multiplication.\n"
"// 	The multiplies within it may all be accomplished using additions and\n"
"// 	subtractions.\n"
"//\n"
"// 	Let's see how this is done.  Given x[n] and x[n+2], cause thats the\n"
"// 	stage we are working on, with i_sync true for x[0] being input,\n"
"// 	produce the output:\n"
"//\n"
"// 	y[n  ] = x[n] + x[n+2]\n"
"// 	y[n+2] = (x[n] - x[n+2]) * e^{-j2pi n/2}	(forward transform)\n"
"// 	       = (x[n] - x[n+2]) * -j^n\n"
"//\n"
"// 	y[n].r = x[n].r + x[n+2].r	(This is the easy part)\n"
"// 	y[n].i = x[n].i + x[n+2].i\n"
"//\n"
"// 	y[2].r = x[0].r - x[2].r\n"
"// 	y[2].i = x[0].i - x[2].i\n"
"//\n"
"// 	y[3].r =   (x[1].i - x[3].i)		(forward transform)\n"
"// 	y[3].i = - (x[1].r - x[3].r)\n"
"//\n"
"// 	y[3].r = - (x[1].i - x[3].i)		(inverse transform)\n"
"// 	y[3].i =   (x[1].r - x[3].r)		(INVERSE = 1)\n"
// "//\n"
// "//	When the FFT is run in the two samples per clock mode, this quarter\n"
// "//	stage will operate on either x[0] and x[2] (ODD = 0), or x[1] and\n"
// "//	x[3] (ODD = 1).  In all other cases, it will operate on all four\n"
// "//	values.\n"
"//\n%s"
"//\n",
		(dbg)?"_dbg":"", prjname, creator);
	fprintf(fp, "%s", cpyleft);
	fprintf(fp, "//\n//\n`default_nettype\tnone\n//\n");

	std::string	resetw("i_reset");
	if (async_reset)
		resetw = std::string("i_areset_n");

	fprintf(fp,
"module\tqtrstage%s(i_clk, %s, i_ce, i_sync, i_data, o_data, o_sync%s);\n"
	"\tparameter	IWIDTH=%d, OWIDTH=IWIDTH+1;\n"
	"\tparameter\tLGWIDTH=%d, INVERSE=0,SHIFT=0;\n"
	"\tinput\twire				i_clk, %s, i_ce, i_sync;\n"
	"\tinput\twire	[(2*IWIDTH-1):0]	i_data;\n"
	"\toutput\treg	[(2*OWIDTH-1):0]	o_data;\n"
	"\toutput\treg				o_sync;\n"
		"\t\n", (dbg)?"_dbg":"", resetw.c_str(),
		(dbg)?", o_dbg":"", TST_QTRSTAGE_IWIDTH,
		TST_QTRSTAGE_LGWIDTH, resetw.c_str());
	if (dbg) { fprintf(fp, "\toutput\twire\t[33:0]\t\t\to_dbg;\n"
		"\tassign\to_dbg = { ((o_sync)&&(i_ce)), i_ce, o_data[(2*OWIDTH-1):(2*OWIDTH-16)],\n"
			"\t\t\t\t\to_data[(OWIDTH-1):(OWIDTH-16)] };\n"
"\n");
	}

	fprintf(fp,
	"\treg\t	wait_for_sync;\n"
	"\treg\t[2:0]	pipeline;\n"
"\n"
	"\treg\tsigned [(IWIDTH):0]	sum_r, sum_i, diff_r, diff_i;\n"
"\n"
	"\treg\t[(2*OWIDTH-1):0]\tob_a;\n"
	"\twire\t[(2*OWIDTH-1):0]\tob_b;\n"
	"\treg\t[(OWIDTH-1):0]\t\tob_b_r, ob_b_i;\n"
	"\tassign\tob_b = { ob_b_r, ob_b_i };\n"
"\n"
	"\treg\t[(LGWIDTH-1):0]\t\tiaddr;\n"
	"\treg\t[(2*IWIDTH-1):0]\timem\t[0:1];\n"
"\n"
	"\twire\tsigned\t[(IWIDTH-1):0]\timem_r, imem_i;\n"
	"\tassign\timem_r = imem[1][(2*IWIDTH-1):(IWIDTH)];\n"
	"\tassign\timem_i = imem[1][(IWIDTH-1):0];\n"
"\n"
	"\twire\tsigned\t[(IWIDTH-1):0]\ti_data_r, i_data_i;\n"
	"\tassign\ti_data_r = i_data[(2*IWIDTH-1):(IWIDTH)];\n"
	"\tassign\ti_data_i = i_data[(IWIDTH-1):0];\n"
"\n"
	"\treg	[(2*OWIDTH-1):0]	omem [0:1];\n"
"\n");

	fprintf(fp, "\t//\n"
	"\t// Round our output values down to OWIDTH bits\n"
	"\t//\n");

	fprintf(fp,
	"\twire\tsigned\t[(OWIDTH-1):0]\trnd_sum_r, rnd_sum_i,\n"
	"\t\t\trnd_diff_r, rnd_diff_i, n_rnd_diff_r, n_rnd_diff_i;\n"
	"\t%s #(IWIDTH+1,OWIDTH,SHIFT)\tdo_rnd_sum_r(i_clk, i_ce,\n"
	"\t\t\t\tsum_r, rnd_sum_r);\n\n", rnd_string);
	fprintf(fp,
	"\t%s #(IWIDTH+1,OWIDTH,SHIFT)\tdo_rnd_sum_i(i_clk, i_ce,\n"
	"\t\t\t\tsum_i, rnd_sum_i);\n\n", rnd_string);
	fprintf(fp,
	"\t%s #(IWIDTH+1,OWIDTH,SHIFT)\tdo_rnd_diff_r(i_clk, i_ce,\n"
	"\t\t\t\tdiff_r, rnd_diff_r);\n\n", rnd_string);
	fprintf(fp,
	"\t%s #(IWIDTH+1,OWIDTH,SHIFT)\tdo_rnd_diff_i(i_clk, i_ce,\n"
	"\t\t\t\tdiff_i, rnd_diff_i);\n\n", rnd_string);
	fprintf(fp, "\tassign n_rnd_diff_r = - rnd_diff_r;\n"
		"\tassign n_rnd_diff_i = - rnd_diff_i;\n");
	fprintf(fp,
	"\tinitial wait_for_sync = 1\'b1;\n"
	"\tinitial iaddr = 0;\n");
	if (async_reset)
		fprintf(fp,
			"\talways @(posedge i_clk, negedge i_areset_n)\n"
			"\tif (!i_areset_n)\n");
	else
		fprintf(fp,
	"\talways @(posedge i_clk)\n"
	"\tif (i_reset)\n");

	fprintf(fp, "\tbegin\n"
		"\t\twait_for_sync <= 1\'b1;\n"
		"\t\tiaddr <= 0;\n"
	"\tend else if ((i_ce)&&((!wait_for_sync)||(i_sync)))\n"
	"\tbegin\n"
		"\t\tiaddr <= iaddr + 1\'b1;\n"
		"\t\twait_for_sync <= 1\'b0;\n"
	"\tend\n\n"
	"\talways @(posedge i_clk)\n"
	"\tif (i_ce)\n"
	"\tbegin\n"
		"\t\timem[0] <= i_data;\n"
		"\t\timem[1] <= imem[0];\n"
	"\tend\n"
	"\n\n");
	fprintf(fp,
	"\t// Note that we don\'t check on wait_for_sync or i_sync here.\n"
	"\t// Why not?  Because iaddr will always be zero until after the\n"
	"\t// first i_ce, so we are safe.\n"
	"\tinitial pipeline = 3\'h0;\n");

	if (async_reset)
		fprintf(fp,
	"\talways\t@(posedge i_clk, negedge i_areset_n)\n"
	"\tif (!i_areset_n)\n");
	else
		fprintf(fp,
	"\talways\t@(posedge i_clk)\n"
	"\tif (i_reset)\n");

	fprintf(fp,
		"\t\tpipeline <= 3\'h0;\n"
	"\telse if (i_ce) // is our pipeline process full?  Which stages?\n"
		"\t\tpipeline <= { pipeline[1:0], iaddr[1] };\n\n");
	fprintf(fp,
	"\t// This is the pipeline[-1] stage, pipeline[0] will be set next.\n"
	"\talways\t@(posedge i_clk)\n"
	"\tif ((i_ce)&&(iaddr[1]))\n"
	"\tbegin\n"
		"\t\tsum_r  <= imem_r + i_data_r;\n"
		"\t\tsum_i  <= imem_i + i_data_i;\n"
		"\t\tdiff_r <= imem_r - i_data_r;\n"
		"\t\tdiff_i <= imem_i - i_data_i;\n"
	"\tend\n\n");
	fprintf(fp,
	"\t// pipeline[1] takes sum_x and diff_x and produces rnd_x\n\n");

	fprintf(fp,
	"\t// Now for pipeline[2].  We can actually do this at all i_ce\n"
	"\t// clock times, since nothing will listen unless pipeline[3]\n"
	"\t// on the next clock.  Thus, we simplify this logic and do\n"
	"\t// it independent of pipeline[2].\n"
	"\talways\t@(posedge i_clk)\n"
	"\tif (i_ce)\n"
	"\tbegin\n"
		"\t\tob_a <= { rnd_sum_r, rnd_sum_i };\n"
		"\t\t// on Even, W = e^{-j2pi 1/4 0} = 1\n"
		"\t\tif (!iaddr[0])\n"
		"\t\tbegin\n"
		"\t\t\tob_b_r <= rnd_diff_r;\n"
		"\t\t\tob_b_i <= rnd_diff_i;\n"
		"\t\tend else if (INVERSE==0) begin\n"
		"\t\t\t// on Odd, W = e^{-j2pi 1/4} = -j\n"
		"\t\t\tob_b_r <=   rnd_diff_i;\n"
		"\t\t\tob_b_i <= n_rnd_diff_r;\n"
		"\t\tend else begin\n"
		"\t\t\t// on Odd, W = e^{j2pi 1/4} = j\n"
		"\t\t\tob_b_r <= n_rnd_diff_i;\n"
		"\t\t\tob_b_i <=   rnd_diff_r;\n"
		"\t\tend\n"
	"\tend\n\n");
	fprintf(fp,
	"\talways\t@(posedge i_clk)\n"
	"\tif (i_ce)\n"
	"\tbegin // In sequence, clock = 3\n"
		"\t\tomem[0] <= ob_b;\n"
		"\t\tomem[1] <= omem[0];\n"
		"\t\tif (pipeline[2])\n"
			"\t\t\to_data <= ob_a;\n"
		"\t\telse\n"
			"\t\t\to_data <= omem[1];\n"
	"\tend\n\n");

	fprintf(fp,
	"\tinitial\to_sync = 1\'b0;\n");

	if (async_reset)
		fprintf(fp,
	"\talways\t@(posedge i_clk, negedge i_areset_n)\n"
	"\tif (!i_areset_n)\n");
	else
		fprintf(fp,
	"\talways\t@(posedge i_clk)\n"
	"\tif (i_reset)\n");
	fprintf(fp,
		"\t\to_sync <= 1\'b0;\n"
		"\telse if (i_ce)\n"
			"\t\to_sync <= (iaddr[2:0] == 3'b101);\n\n");

	if (formal_property_flag) {
		fprintf(fp,
"`ifdef	FORMAL\n"
	"\t// Formal declarations\n"
	"\t// {{{\n"
	"\treg				f_past_valid;\n"
	"\treg	signed [IWIDTH-1:0]	f_piped_real	[0:7];\n"
	"\treg	signed [IWIDTH-1:0]	f_piped_imag	[0:7];\n"
	"\treg				f_rsyncd;\n"
	"\twire				f_syncd;\n"
	"\treg	[1:0]			f_state;\n"
	"\twire	signed [OWIDTH-1:0]	f_o_real, f_o_imag;\n"
	"\t// }}}\n"
"\n"
	"\tinitial	f_past_valid = 1'b0;\n"
	"\talways @(posedge i_clk)\n"
	"\t	f_past_valid <= 1'b1;\n"
"\n"
"`ifdef	QTRSTAGE\n"
	"\talways @(posedge i_clk)\n"
	"\t	assume((i_ce)||($past(i_ce))||($past(i_ce,2)));\n"
"`endif\n"
"\n"
	"\t// The below logic only works if the rounding stage does nothing\n"
	"\tinitial	assert(IWIDTH+1 == OWIDTH);\n"
"\n"
"\n"
	"\talways @(posedge i_clk)\n"
	"\tif (i_ce)\n"
	"\tbegin\n"
	"\t	f_piped_real[0] <= i_data[2*IWIDTH-1:IWIDTH];\n"
	"\t	f_piped_imag[0] <= i_data[  IWIDTH-1:0];\n"
"\n"
	"\t	f_piped_real[1] <= f_piped_real[0];\n"
	"\t	f_piped_imag[1] <= f_piped_imag[0];\n"
"\n"
	"\t	f_piped_real[2] <= f_piped_real[1];\n"
	"\t	f_piped_imag[2] <= f_piped_imag[1];\n"
"\n"
	"\t	f_piped_real[3] <= f_piped_real[2];\n"
	"\t	f_piped_imag[3] <= f_piped_imag[2];\n"
"\n"
	"\t	f_piped_real[4] <= f_piped_real[3];\n"
	"\t	f_piped_imag[4] <= f_piped_imag[3];\n"
"\n"
	"\t	f_piped_real[5] <= f_piped_real[4];\n"
	"\t	f_piped_imag[5] <= f_piped_imag[4];\n"
"\n"
	"\t	f_piped_real[6] <= f_piped_real[5];\n"
	"\t	f_piped_imag[6] <= f_piped_imag[5];\n"
"\n"
	"\t	f_piped_real[7] <= f_piped_real[6];\n"
	"\t	f_piped_imag[7] <= f_piped_imag[6];\n"
	"\tend\n"
"\n"
"\n"
	"\tinitial	f_rsyncd = 0;\n");

	if (async_reset)
		fprintf(fp, "\talways @(posedge i_clk, negedge i_areset_n)\n\tif (!i_areset_n)\n");
	else
		fprintf(fp, "\talways @(posedge i_clk)\n\tif (i_reset)\n");
	fprintf(fp,
	"\t	f_rsyncd <= 1'b0;\n"
	"\telse if (!f_rsyncd)\n"
	"\t	f_rsyncd <= (o_sync);\n"
	"\tassign	f_syncd = (f_rsyncd)||(o_sync);\n"
"\n"
"\n"
	"\tinitial	f_state = 0;\n");
	if (async_reset)
		fprintf(fp, "\talways @(posedge i_clk, negedge i_areset_n)\n\tif (!i_areset_n)\n");
	else
		fprintf(fp, "\talways @(posedge i_clk)\n\tif (i_reset)\n");

	fprintf(fp,
	"\t	f_state <= 0;\n"
	"\telse if ((i_ce)&&((!wait_for_sync)||(i_sync)))\n"
	"\t	f_state <= f_state + 1;\n"
"\n"
	"\talways @(*)\n"
	"\tif (f_state != 0)\n"
	"\t	assume(!i_sync);\n"
"\n"
	"\talways @(posedge i_clk)\n"
	"\t	assert(f_state[1:0] == iaddr[1:0]);\n"
"\n"
	"\tassign			f_o_real = o_data[2*OWIDTH-1:OWIDTH];\n"
	"\tassign			f_o_imag = o_data[  OWIDTH-1:0];\n"
"\n"
	"\talways @(posedge i_clk)\n"
	"\tif (f_state == 2'b11)\n"
	"\tbegin\n"
	"\t	assume(f_piped_real[0] != 3'sb100);\n"
	"\t	assume(f_piped_real[2] != 3'sb100);\n"
	"\t	assert(sum_r  == f_piped_real[2] + f_piped_real[0]);\n"
	"\t	assert(sum_i  == f_piped_imag[2] + f_piped_imag[0]);\n"
"\n"
	"\t	assert(diff_r == f_piped_real[2] - f_piped_real[0]);\n"
	"\t	assert(diff_i == f_piped_imag[2] - f_piped_imag[0]);\n"
	"\tend\n"
"\n"
	"\talways @(posedge i_clk)\n"
	"\tif ((f_state == 2'b00)&&((f_syncd)||(iaddr >= 4)))\n"
	"\tbegin\n"
	"\t	assert(rnd_sum_r  == f_piped_real[3]+f_piped_real[1]);\n"
	"\t	assert(rnd_sum_i  == f_piped_imag[3]+f_piped_imag[1]);\n"
	"\t	assert(rnd_diff_r == f_piped_real[3]-f_piped_real[1]);\n"
	"\t	assert(rnd_diff_i == f_piped_imag[3]-f_piped_imag[1]);\n"
	"\tend\n"
"\n"
	"\talways @(posedge i_clk)\n"
	"\tif ((f_state == 2'b10)&&(f_syncd))\n"
	"\tbegin\n"
	"\t	// assert(o_sync);\n"
	"\t	assert(f_o_real == f_piped_real[5] + f_piped_real[3]);\n"
	"\t	assert(f_o_imag == f_piped_imag[5] + f_piped_imag[3]);\n"
	"\tend\n"
"\n"
	"\talways @(posedge i_clk)\n"
	"\tif ((f_state == 2'b11)&&(f_syncd))\n"
	"\tbegin\n"
	"\t	assert(!o_sync);\n"
	"\t	assert(f_o_real == f_piped_real[5] + f_piped_real[3]);\n"
	"\t	assert(f_o_imag == f_piped_imag[5] + f_piped_imag[3]);\n"
	"\tend\n"
"\n"
	"\talways @(posedge i_clk)\n"
	"\tif ((f_state == 2'b00)&&(f_syncd))\n"
	"\tbegin\n"
	"\t	assert(!o_sync);\n"
	"\t	assert(f_o_real == f_piped_real[7] - f_piped_real[5]);\n"
	"\t	assert(f_o_imag == f_piped_imag[7] - f_piped_imag[5]);\n"
	"\tend\n"
"\n"
	"\talways @(*)\n"
	"\tif ((iaddr[2:0] == 0)&&(!wait_for_sync))\n"
	"\t	assume(i_sync);\n"
"\n"
	"\talways @(*)\n"
	"\tif (wait_for_sync)\n"
	"\t	assert((iaddr == 0)&&(f_state == 2'b00)&&(!o_sync)&&(!f_rsyncd));\n"
"\n"
	"\talways @(posedge i_clk)\n");

	if (async_reset)
		fprintf(fp, "\tif ((f_past_valid && i_areset_n)&&($past(i_ce))"
			"&&($past(i_sync))&&($past(i_areset_n)))\n");
	else
		fprintf(fp, "\tif ((f_past_valid)&&($past(i_ce))"
			"&&($past(i_sync))&&(!$past(i_reset)))\n");

	fprintf(fp,
	"\t	assert(!wait_for_sync);\n"
"\n"
	"\talways @(posedge i_clk)\n"
	"\tif ((f_state == 2'b01)&&(f_syncd))\n"
	"\tbegin\n"
	"\t	assert(!o_sync);\n"
	"\t	if (INVERSE)\n"
	"\t	begin\n"
	"\t		assert(f_o_real == -f_piped_imag[7]+f_piped_imag[5]);\n"
	"\t		assert(f_o_imag ==  f_piped_real[7]-f_piped_real[5]);\n"
	"\t	end else begin\n"
	"\t		assert(f_o_real ==  f_piped_imag[7]-f_piped_imag[5]);\n"
	"\t		assert(f_o_imag == -f_piped_real[7]+f_piped_real[5]);\n"
	"\t	end\n"
	"\tend\n"
"\n"
"`endif\n");
	}

	fprintf(fp, "endmodule\n");
}
// }}}

// build_sngllast
// {{{
void	build_sngllast(const char *fname, const bool async_reset = false) {
	FILE	*fp = fopen(fname, "w");
	if (NULL == fp) {
		fprintf(stderr, "Could not open \'%s\' for writing\n", fname);
		perror("O/S Err was:");
		return;
	}

	std::string	resetw("i_reset");
	if (async_reset)
		resetw = std::string("i_areset_n");

	fprintf(fp,
SLASHLINE
"//\n"
"// Filename:\tlaststage.v\n"
"// {{{\n" // "}}}"
"// Project:	%s\n"
"//\n"
"// Purpose:	This is part of an FPGA implementation that will process\n"
"//		the final stage of a decimate-in-frequency FFT, running\n"
"//	through the data at one sample per clock.\n"
"//\n"
"//\n%s"
"//\n", prjname, creator);
	fprintf(fp, "%s", cpyleft);
	fprintf(fp, "//\n//\n`default_nettype\tnone\n//\n");

	fprintf(fp,
"module	laststage #(\n"
	"\t\t// {{{\n"
	"\t\tparameter IWIDTH=16,OWIDTH=IWIDTH+1, SHIFT=0\n"
	"\t\t// }}}\n"
	"\t) (\n"
	"\t\t// {{{\n"
	"\t\tinput\twire			i_clk, %s, i_ce, i_sync,\n"
	"\t\tinput\twire  [(2*IWIDTH-1):0]	i_val,\n"
	"\t\toutput\twire [(2*OWIDTH-1):0]	o_val,\n"
	"\t\toutput\treg			o_sync\n"
	"\t\t// }}}\n"
	"\t);\n",
		resetw.c_str());

	fprintf(fp,
	"\t// Local declarations\n"
	"\t// {{{\n"
	"\treg	signed	[(IWIDTH-1):0]	m_r, m_i;\n"
	"\twire	signed	[(IWIDTH-1):0]	i_r, i_i;\n"
"\n"
	"\t// Don't forget that we accumulate a bit by adding two values\n"
	"\t// together. Therefore our intermediate value must have one more\n"
	"\t// bit than the two originals.\n"
	"\treg	signed	[(IWIDTH):0]	rnd_r, rnd_i, sto_r, sto_i;\n"
	"\treg				wait_for_sync, stage;\n"
	"\treg		[1:0]		sync_pipe;\n"
	"\twire	signed	[(OWIDTH-1):0]	o_r, o_i;\n"
	"\t// }}}\n"
"\n"
"	assign	i_r = i_val[(2*IWIDTH-1):(IWIDTH)]; \n"
"	assign	i_i = i_val[(IWIDTH-1):0]; \n"
"\n"
	"\t// wait_for_sync, stage\n"
	"\t// {{{\n"
"	initial	wait_for_sync = 1'b1;\n"
"	initial	stage         = 1'b0;\n");

	if (async_reset)
		fprintf(fp, "\talways @(posedge i_clk, negedge i_areset_n)\n\tif (!i_areset_n)\n");
	else
		fprintf(fp, "\talways @(posedge i_clk)\n\tif (i_reset)\n");
	fprintf(fp,
	"\tbegin\n"
	"\t	wait_for_sync <= 1'b1;\n"
	"\t	stage         <= 1'b0;\n"
	"\tend else if ((i_ce)&&((!wait_for_sync)||(i_sync))&&(!stage))\n"
	"\tbegin\n"
	"\t	wait_for_sync <= 1'b0;\n"
	"\t	//\n"
	"\t	stage <= 1'b1;\n"
	"\t	//\n"
	"\tend else if (i_ce)\n"
	"\t	stage <= 1'b0;\n"
	"\t// }}}\n\n");

	fprintf(fp,
	"\t// sync_pipe\n"
	"\t// {{{\n"
	"\tinitial\tsync_pipe = 0;\n");
	if (async_reset)
		fprintf(fp,
		"\talways @(posedge i_clk, negedge i_areset_n)\n"
		"\tif (!i_areset_n)\n");
	else
		fprintf(fp,
		"\talways @(posedge i_clk)\n"
		"\tif (i_reset)\n");

	fprintf(fp,
		"\t\tsync_pipe <= 0;\n"
		"\telse if (i_ce)\n"
		"\t\tsync_pipe <= { sync_pipe[0], i_sync };\n"
		"\t// }}}\n\n");

	fprintf(fp,
	"\t// o_sync\n"
	"\t// {{{\n"
	"\tinitial\to_sync = 1\'b0;\n");
	if (async_reset)
		fprintf(fp,
		"\talways @(posedge i_clk, negedge i_areset_n)\n"
		"\tif (!i_areset_n)\n");
	else
		fprintf(fp,
		"\talways @(posedge i_clk)\n"
		"\tif (i_reset)\n");

	fprintf(fp,
		"\t\to_sync <= 1\'b0;\n"
		"\telse if (i_ce)\n"
		"\t\to_sync <= sync_pipe[1];\n"
		"\t// }}}\n\n");

	fprintf(fp,
	"\t// m_r, m_i, rnd_r, rnd_i\n"
	"\t// {{{\n"
	"\talways @(posedge i_clk)\n"
	"\tif (i_ce)\n"
	"\tbegin\n"
	"\t\tif (!stage)\n"
	"\t\tbegin\n"
	"\t\t\t// Clock 1\n"
	"\t\t\tm_r <= i_r;\n"
	"\t\t\tm_i <= i_i;\n"
	"\t\t\t// Clock 3\n"
	"\t\t\trnd_r <= sto_r;\n"
	"\t\t\trnd_i <= sto_i;\n"
	"\t\t\t//\n"
	"\t\tend else begin\n"
	"\t\t\t// Clock 2\n"
	"\t\t\trnd_r <= m_r + i_r;\n"
	"\t\t\trnd_i <= m_i + i_i;\n"
	"\t\t\t//\n"
	"\t\t\tsto_r <= m_r - i_r;\n"
	"\t\t\tsto_i <= m_i - i_i;\n"
	"\t\t\t//\n"
	"\t\tend\n"
	"\tend\n"
	"\t// }}}\n"
"\n"
"	// Now that we have our results, let's round them and report them\n"
"\n"
	"\t// Round the results, generating o_r, o_i, and thus o_val\n"
	"\t// {{{\n"
"	convround #(IWIDTH+1,OWIDTH,SHIFT) do_rnd_r(i_clk, i_ce, rnd_r, o_r);\n"
"	convround #(IWIDTH+1,OWIDTH,SHIFT) do_rnd_i(i_clk, i_ce, rnd_i, o_i);\n"
"\n"
"	assign	o_val  = { o_r, o_i };\n"
	"\t// }}}\n"
"\n");

	fprintf(fp,
SLASHLINE
SLASHLINE
SLASHLINE
"//\n"
"// Formal properties\n"
"// {{{\n"
SLASHLINE
SLASHLINE
SLASHLINE );
	if (formal_property_flag) {
		fprintf(fp,
	"`ifdef	FORMAL\n"
		"\t// Local formal declarations\n"
		"\t// {{{\n"
		"\treg	f_past_valid;\n"
		"\twire	f_syncd;\n"
		"\treg	f_rsyncd;\n"
		"\treg	f_state;\n"
		"\t// }}}\n"
	"\n"
		"\tinitial	f_past_valid = 1'b0;\n"
		"\talways @(posedge i_clk)\n"
		"\t	f_past_valid <= 1'b1;\n"
	"\n"
	"`ifdef	LASTSTAGE\n"
		"\talways @(posedge i_clk)\n"
		"\t	assume((i_ce)||($past(i_ce))||($past(i_ce,2)));\n"
	"`endif\n"
	"\n"
		"\tinitial	assert(IWIDTH+1 == OWIDTH);\n"
	"\n"
		"\treg	signed	[IWIDTH-1:0]	f_piped_real	[0:3];\n"
		"\treg	signed	[IWIDTH-1:0]	f_piped_imag	[0:3];\n"
		"\talways @(posedge i_clk)\n"
		"\tif (i_ce)\n"
		"\tbegin\n"
		"\t	f_piped_real[0] <= i_val[2*IWIDTH-1:IWIDTH];\n"
		"\t	f_piped_imag[0] <= i_val[  IWIDTH-1:0];\n"
	"\n"
		"\t	f_piped_real[1] <= f_piped_real[0];\n"
		"\t	f_piped_imag[1] <= f_piped_imag[0];\n"
	"\n"
		"\t	f_piped_real[2] <= f_piped_real[1];\n"
		"\t	f_piped_imag[2] <= f_piped_imag[1];\n"
	"\n"
		"\t	f_piped_real[3] <= f_piped_real[2];\n"
		"\t	f_piped_imag[3] <= f_piped_imag[2];\n"
		"\tend\n"
	"\n"
		"\tinitial	f_rsyncd	= 0;\n"
		"\talways @(posedge i_clk)\n"
		"\tif (i_reset)\n"
		"\t	f_rsyncd <= 1'b0;\n"
		"\telse if (!f_rsyncd)\n"
		"\t	f_rsyncd <= o_sync;\n"
		"\tassign	f_syncd = (f_rsyncd)||(o_sync);\n"
	"\n"
		"\tinitial	f_state = 0;\n"
		"\talways @(posedge i_clk)\n"
		"\tif (i_reset)\n"
		"\t	f_state <= 0;\n"
		"\telse if ((i_ce)&&((!wait_for_sync)||(i_sync)))\n"
		"\t	f_state <= f_state + 1;\n"
	"\n"
		"\talways @(*)\n"
		"\tif (f_state != 0)\n"
		"\t	assume(!i_sync);\n"
	"\n"
		"\talways @(*)\n"
		"\t	assert(stage == f_state[0]);\n"
	"\n"
		"\talways @(posedge i_clk)\n"
		"\tif ((f_state == 1'b1)&&(f_syncd))\n"
		"\tbegin\n"
		"\t	assert(o_r == f_piped_real[2] + f_piped_real[1]);\n"
		"\t	assert(o_i == f_piped_imag[2] + f_piped_imag[1]);\n"
		"\tend\n"
	"\n"
		"\talways @(posedge i_clk)\n"
		"\tif ((f_state == 1'b0)&&(f_syncd))\n"
		"\tbegin\n"
		"\t	assert(!o_sync);\n"
		"\t	assert(o_r == f_piped_real[3] - f_piped_real[2]);\n"
		"\t	assert(o_i == f_piped_imag[3] - f_piped_imag[2]);\n"
		"\tend\n"
	"\n"
		"\talways @(*)\n"
		"\tif (wait_for_sync)\n"
		"\tbegin\n"
		"\t	assert(!f_rsyncd);\n"
		"\t	assert(!o_sync);\n"
		"\t	assert(f_state == 0);\n"
		"\tend\n\n");
	} else {
		fprintf(fp,
		"//\n"
		"// Formal properties have not been included in this build\n"
		"//\n");
	}

	fprintf(fp,
"`endif // FORMAL\n"
"// }}}\n"
"endmodule\n");

	fclose(fp);
}
// }}}

// usage()
// {{{
void	usage(void) {
	fprintf(stderr,
"USAGE:\tfftgen [-f <size>] [-d dir] [-c cbits] [-n nbits] [-m mxbits] [-s]\n"
// "\tfftgen -i\n"
"\t-1\tBuild a normal FFT, running at one clock per complex sample, or\n"
"\t\t(for a real FFT) at one clock per two real input samples.\n"
"\t-A\t(Experimental) Use a negative edged asynchronous reset.\n"
"\t-a <hdrname>  Create a header of information describing the built-in\n"
"\t\tparameters, useful for module-level testing with Verilator\n"
"\t-c <cbits>\tCauses all internal complex coefficients to be\n"
"\t\tlonger than the corresponding data bits, to help avoid\n"
"\t\tcoefficient truncation errors.  The default is %d bits longer\n"
"\t\tthan the data bits.\n"
"\t-d <dir>  Places all of the generated verilog files into <dir>.\n"
"\t\tThe default is a subdirectory of the current directory\n"
"\t\tnamed %s.\n"
"\t-f <size>  Sets the size of the FFT as the number of complex\n"
"\t\tsamples input to the transform.  (No default value, this is\n"
"\t\ta required parameter.)\n"
"\t-i\tAn inverse FFT, meaning that the coefficients are\n"
"\t\tgiven by e^{ j 2 pi k/N n }.  The default is a forward FFT, with\n"
"\t\tcoefficients given by e^{ -j 2 pi k/N n }.\n"
"\t-k #\tSets # clocks per sample, used to minimize multiplies.  Also\n"
"\t\tsets one sample in per i_ce clock (opt -1)\n"
"\t-m <mxbits>\tSets the maximum bit width that the FFT should ever\n"
"\t\tproduce.  Internal values greater than this value will be\n"
"\t\ttruncated to this value.  (The default value grows the input\n"
"\t\tsize by one bit for every two FFT stages.)\n"
"\t-n <nbits>\tSets the bitwidth for values coming into the (i)FFT.\n"
"\t\tThe default is %d bits input for each component of the two\n"
"\t\tcomplex values into the FFT.\n"
"\t-p <nmpy>  Sets the number of hardware multiplies (DSPs) to use, versus\n"
"\t\tshift-add emulation.  The default is not to use any hardware\n"
"\t\tmultipliers.\n"
"\t-r\tBuild a real-FFT at four input points per sample, rather than a\n"
"\t\tcomplex FFT.  (Default is a Complex FFT.)\n"
"\t\tThis option is a place-holder.  The real-FFT has not (yet) been\n"
"\t\timplemented.\n"
"\t-s\tSkip the final bit reversal stage.  This is useful in\n"
"\t\talgorithms that need to apply a filter without needing to do\n"
"\t\tbin shifting, as these algorithms can, with this option, just\n"
"\t\tmultiply by a bit reversed correlation sequence and then\n"
"\t\tinverse FFT the (still bit reversed) result.  (You would need\n"
"\t\ta decimation in time inverse to do this, which this program does\n"
"\t\tnot yet provide.)\n"
"\t-S\tInclude the final bit reversal stage (default).\n"
"\t-x <xtrabits>\tUse this many extra bits internally, before any final\n"
"\t\trounding or truncation of the answer to the final number of\n"
"\t\tbits.  The default is to use %d extra bits internally.\n",
/*
"\t-0\tA forward FFT (default), meaning that the coefficients are\n"
"\t\tgiven by e^{-j 2 pi k/N n }.\n"
"\t-1\tAn inverse FFT, meaning that the coefficients are\n"
"\t\tgiven by e^{ j 2 pi k/N n }.\n",
*/
	DEF_XTRACBITS, DEF_COREDIR, DEF_NBITSIN, DEF_XTRAPBITS);
}
// }}}

// Features still needed:
//	Interactivity.
int main(int argc, char **argv) {
	int	fftsize = -1, lgsize = -1;
	int	nbitsin = DEF_NBITSIN, xtracbits = DEF_XTRACBITS,
			nummpy=DEF_NMPY, nmpypstage=6, mpy_stages;
	int	nbitsout, maxbitsout = -1, xtrapbits=DEF_XTRAPBITS, ckpce = 0;
	const char *EMPTYSTR = "";
	bool	bitreverse = true, inverse=false,
		verbose_flag = false,
		single_clock = true,
		real_fft = false,
		async_reset = false;
	FILE	*vmain;
	std::string	coredir = DEF_COREDIR, cmdline = "", hdrname = "";
	ROUND_T	rounding = RND_CONVERGENT;
	// ROUND_T	rounding = RND_HALFUP;

	bool	dbg = false;
	int	dbgstage = 128;

	// Argument processing
	// {{{
	if (argc <= 1)
		usage();

	// Copy the original command line before we mess with it
	cmdline = argv[0];
	for(int argn=1; argn<argc; argn++) {
		cmdline += " ";
		cmdline += argv[argn];
	}

	{ int c;
	while((c = getopt(argc, argv, "12Aa:c:d:D:f:hik:m:n:p:rsSx:v")) != -1) {
		switch(c) {
		case '1':	single_clock = true;  break;
		case '2':	single_clock = false; break;
		case 'A':	async_reset  = true;  break;
		case 'a':	hdrname = strdup(optarg);	break;
		case 'c':	xtracbits = atoi(optarg);	break;
		case 'd':	coredir = std::string(optarg);	break;
		case 'D':	dbgstage = atoi(optarg);	break;
		case 'f':	fftsize = atoi(optarg);	
				{ int sln = strlen(optarg);
				if (!isdigit(optarg[sln-1])){
					switch(optarg[sln-1]) {
					case 'k': case 'K':
						fftsize <<= 10;
						break;
					case 'm': case 'M':
						fftsize <<= 20;
						break;
					case 'g': case 'G':
						fftsize <<= 30;
						break;
					default:
						printf("ERR: Unknown FFT size, %s!\n", optarg);
						exit(EXIT_FAILURE);
					}
				}} break;
		case 'h':	usage(); exit(EXIT_SUCCESS);	break;
		case 'i':	inverse = true;			break;
		case 'k':	ckpce = atoi(optarg);
				// After CKPCE=3, there's no advantage
				// only bits lost, so keep CKPCE <=3 here
				if (ckpce > 3)
					ckpce = 3;
				single_clock = true;
				break;
		case 'm':	maxbitsout = atoi(optarg);	break;
		case 'n':	nbitsin = atoi(optarg);		break;
		case 'p':	nummpy = atoi(optarg);		break;
		case 'r':	real_fft = true;		break;
		case 'S':	bitreverse = true;		break;
		case 's':	bitreverse = false;		break;
		case 'x':	xtrapbits = atoi(optarg);	break;
		case 'v':	verbose_flag = true;		break;
		// case 'z':	variable_size = true;		break;
		default:
			printf("Unknown argument, -%c\n", c);
			usage();
			exit(EXIT_FAILURE);
		}
	}}
	// }}}

	// verbose: Repeat back our chosen arguments
	// {{{
	if (verbose_flag) {
		if (inverse)
			printf("Building a %d point inverse FFT module, with %s outputs\n",
				fftsize,
				(real_fft)?"real ":"complex");
		else
			printf("Building a %d point %sforward FFT module\n",
				fftsize,
				(real_fft)?"real ":"");
		if (!single_clock)
			printf("  that accepts two inputs per clock\n");
		if (async_reset)
			printf("  using a negative logic ASYNC reset\n");

		printf("The core will be placed into the %s/ directory\n", coredir.c_str());

		if (hdrname[0])
			printf("A C header file, %s, will be written capturing these\n"
				"options for a Verilator testbench\n", 
					hdrname.c_str());
		// nummpy
		// xtrapbits
	}
	// }}}

	// Argument sanity checks
	// {{{
	if (real_fft) {
		printf("The real FFT option is not implemented yet, but still on\nmy to do list.  Please try again later.\n");
		exit(EXIT_FAILURE);
	}

	if (ckpce < 1)
		ckpce = 1;
	if (!bitreverse) {
		printf("WARNING: While I can skip the bit reverse stage, the code to do\n");
		printf("an inverse FFT on a bit--reversed input has not yet been\n");
		printf("built.\n");
	}

	if ((lgsize < 0)&&(fftsize > 1)) {
		for(lgsize=1; (1<<lgsize) < fftsize; lgsize++)
			;
	}

	if (fftsize <= 0) {
		printf("ERROR: Invalid size.  FFT size (%d) may not be <= 0\n",
			fftsize);
		exit(EXIT_FAILURE);
	}

	if (nbitsin < 1) {
		printf("ERROR: Not enough input bits, %d >= 1\n", nbitsin);
		exit(EXIT_FAILURE);
	} else if (nbitsin>48) {
		printf("ERROR: Too many input bits, %d is greater than 48\n",
			nbitsin);
		exit(EXIT_FAILURE);
	}


	if (nextlg(fftsize) != fftsize) {
		fprintf(stderr, "ERR: FFTSize (%d) *must* be a power of two\n",
				fftsize);
		exit(EXIT_FAILURE);
	} else if (fftsize < 2) {
		fprintf(stderr, "ERR: Minimum FFTSize is 2, not %d\n",
				fftsize);
		if (fftsize == 1) {
			fprintf(stderr, "You do realize that a 1 point FFT makes very little sense\n");
			fprintf(stderr, "in an FFT operation that handles two samples per clock?\n");
			fprintf(stderr, "If you really need to do an FFT of this size, the output\n");
			fprintf(stderr, "can be connected straight to the input.\n");
		} else {
			fprintf(stderr, "Indeed, a size of %d doesn\'t make much sense to me at all.\n", fftsize);
			fprintf(stderr, "Is such an operation even defined?\n");
		}
		exit(EXIT_FAILURE);
	}
	// }}}

	// nbitsout, bitreverse, and tmp_size
	// {{{
	// Calculate how many output bits we'll have, and what the log
	// based two size of our FFT is.
	{
		int	tmp_size = fftsize;

		// The first stage always accumulates one bit, regardless
		// of whether you need to or not.
		nbitsout = nbitsin + 1;
		tmp_size >>= 1;

		while(tmp_size > 4) {
			nbitsout += 1;
			tmp_size >>= 2;
		}

		if (tmp_size > 1)
			nbitsout ++;

		if (fftsize <= 2)
			bitreverse = false;
	} if ((maxbitsout > 0)&&(nbitsout > maxbitsout))
		nbitsout = maxbitsout;
	// }}}

	// Reflect our bit-width calcualtion
	// {{{
	if (verbose_flag) {
		printf("Output samples will be %d bits wide\n", nbitsout);
		printf("This %sFFT will take %d-bit samples in, and produce %d samples out\n", (inverse)?"i":"", nbitsin, nbitsout);
		if (maxbitsout > 0)
			printf("  Internally, it will allow items to accumulate to %d bits\n", maxbitsout);
		printf("  Twiddle-factors of %d bits will be used\n",
			nbitsin+xtracbits);
		if (!bitreverse)
		printf("  The output will be left in bit-reversed order\n");
	}
	// }}}

	// Figure out how many multiply stages to use, and how many to skip
	// {{{
	if (!single_clock) {
		nmpypstage = 6;
	} else if (ckpce <= 1) {
		nmpypstage = 3;
	} else if (ckpce == 2) {
		nmpypstage = 2;
	} else
		nmpypstage = 1;

	mpy_stages = nummpy / nmpypstage;
	if (mpy_stages > lgval(fftsize)-2)
		mpy_stages = lgval(fftsize)-2;
	// }}}

	// Create an output directory
	// {{{
	{
		struct stat	sbuf;
		if (lstat(coredir.c_str(), &sbuf)==0) {
			if (!S_ISDIR(sbuf.st_mode)) {
				fprintf(stderr, "\'%s\' already exists, and is not a directory!\n", coredir.c_str());
				fprintf(stderr, "I will stop now, lest I overwrite something you care about.\n");
				fprintf(stderr, "To try again, please remove this file.\n");
				exit(EXIT_FAILURE);
			}
		} else
			mkdir(coredir.c_str(), 0755);
		if (access(coredir.c_str(), X_OK|W_OK) != 0) {
			fprintf(stderr, "I have no access to the directory \'%s\'.\n", coredir.c_str());
			exit(EXIT_FAILURE);
		}
	}
	// }}}

	// Write a header file with our chosen parameters
	// {{{
	if (hdrname.length() > 0) {
		FILE	*hdr = fopen(hdrname.c_str(), "w");
		if (hdr == NULL) {
			fprintf(stderr, "ERROR: Cannot open %s to create header file\n", hdrname.c_str());
			perror("O/S Err:");
			exit(EXIT_FAILURE);
		}

		fprintf(hdr,
SLASHLINE
"//\n"
"// Filename:\t%s\n"
"// {{{\n" // "}}}"
"// Project:\t%s\n"
"//\n"
"// Purpose:	This simple header file captures the internal constants\n"
"//		within the FFT that were used to build it, for the purpose\n"
"//	of making C++ integration (and test bench testing) simpler.  That is,\n"
"//	should the FFT change size, this will note that size change and thus\n"
"//	any test bench or other C++ program dependent upon either the size of\n"
"//	the FFT, the number of bits in or out of it, etc., can pick up the\n"
"//	changes in the defines found within this file.\n"
"//\n",
		hdrname.c_str(), prjname);
		fprintf(hdr, "%s", creator);
		fprintf(hdr, "//\n");
		fprintf(hdr, "%s", cpyleft);
		fprintf(hdr, "//\n"
		"//\n"
		"#ifndef %sFFTHDR_H\n"
		"#define %sFFTHDR_H\n"
		"\n"
		"#define\t%sFFT_IWIDTH\t%d\n"
		"#define\t%sFFT_OWIDTH\t%d\n"
		"#define\t%sFFT_LGWIDTH\t%d\n"
		"#define\t%sFFT_SIZE\t(1<<%sFFT_LGWIDTH)\n\n",
			(inverse)?"I":"", (inverse)?"I":"",
			(inverse)?"I":"", nbitsin,
			(inverse)?"I":"", nbitsout,
			(inverse)?"I":"", lgsize,
			(inverse)?"I":"", (inverse)?"I":"");
		if (async_reset)
			fprintf(hdr, "#define\tASYNC_RESETN\n");
		else
			fprintf(hdr, "// #define ASYNC_RESETN\n");
		if (ckpce > 0)
			fprintf(hdr, "#define\t%sFFT_CKPCE\t%d\t// Clocks per CE\n",
				(inverse)?"I":"", ckpce);
		else
			fprintf(hdr, "// Two samples per i_ce\n");
		if (!bitreverse)
			fprintf(hdr, "#define\t%sFFT_SKIPS_BIT_REVERSE\n",
				(inverse)?"I":"");
		if (real_fft)
			fprintf(hdr, "#define\tRL%sFFT\n\n", (inverse)?"I":"");
		if (!single_clock)
			fprintf(hdr, "#define\tDBLCLK%sFFT\n\n", (inverse)?"I":"");
		else
			fprintf(hdr, "// #define\tDBLCLK%sFFT // this FFT takes one input sample per clock\n\n", (inverse)?"I":"");
		if (USE_OLD_MULTIPLY)
			fprintf(hdr, "#define\tUSE_OLD_MULTIPLY\n\n");

		fprintf(hdr, "// Parameters for testing the longbimpy\n");
		fprintf(hdr, "#define\tTST_LONGBIMPY_AW\t%d\n", TST_LONGBIMPY_AW);
#ifdef	TST_LONGBIMPY_BW
		fprintf(hdr, "#define\tTST_LONGBIMPY_BW\t%d\n\n", TST_LONGBIMPY_BW);
#else
		fprintf(hdr, "#define\tTST_LONGBIMPY_BW\tTST_LONGBIMPY_AW\n\n");
#endif

		fprintf(hdr, "// Parameters for testing the shift add multiply\n");
		fprintf(hdr, "#define\tTST_SHIFTADDMPY_AW\t%d\n", TST_SHIFTADDMPY_AW);
#ifdef	TST_SHIFTADDMPY_BW
		fprintf(hdr, "#define\tTST_SHIFTADDMPY_BW\t%d\n\n", TST_SHIFTADDMPY_BW);
#else
		fprintf(hdr, "#define\tTST_SHIFTADDMPY_BW\tTST_SHIFTADDMPY_AW\n\n");
#endif

#define	TST_SHIFTADDMPY_AW	16
#define	TST_SHIFTADDMPY_BW	20	// Leave undefined to match AW
		fprintf(hdr, "// Parameters for testing the butterfly\n");
		fprintf(hdr, "#define\tTST_BUTTERFLY_IWIDTH\t%d\n", TST_BUTTERFLY_IWIDTH);
		fprintf(hdr, "#define\tTST_BUTTERFLY_CWIDTH\t%d\n", TST_BUTTERFLY_CWIDTH);
		fprintf(hdr, "#define\tTST_BUTTERFLY_OWIDTH\t%d\n", TST_BUTTERFLY_OWIDTH);
		fprintf(hdr, "#define\tTST_BUTTERFLY_MPYDELAY\t%d\n\n",
				bflydelay(TST_BUTTERFLY_IWIDTH,
					TST_BUTTERFLY_CWIDTH-TST_BUTTERFLY_IWIDTH));

		fprintf(hdr, "// Parameters for testing the quarter stage\n");
		fprintf(hdr, "#define\tTST_QTRSTAGE_IWIDTH\t%d\n", TST_QTRSTAGE_IWIDTH);
		fprintf(hdr, "#define\tTST_QTRSTAGE_LGWIDTH\t%d\n\n", TST_QTRSTAGE_LGWIDTH);

		fprintf(hdr, "// Parameters for testing the double stage\n");
		fprintf(hdr, "#define\tTST_DBLSTAGE_IWIDTH\t%d\n", TST_DBLSTAGE_IWIDTH);
		fprintf(hdr, "#define\tTST_DBLSTAGE_SHIFT\t%d\n\n", TST_DBLSTAGE_SHIFT);

		fprintf(hdr, "// Parameters for testing the bit reversal stage\n");
		fprintf(hdr, "#define\tTST_DBLREVERSE_LGSIZE\t%d\n\n", TST_DBLREVERSE_LGSIZE);
		fprintf(hdr, "\n" "#endif\n\n");
		fclose(hdr);
	}
	// }}}

	////////////////////////////////////////////////////////////////////////
	//
	// Build FFTMAIN
	// {{{
	// Open the file
	// {{{
	{
		std::string	fname_string;

		fname_string = coredir;
		fname_string += "/";
		if (inverse) fname_string += "i";
		fname_string += "fftmain.v";

		vmain = fopen(fname_string.c_str(), "w");
		if (NULL == vmain) {
			fprintf(stderr, "Could not open \'%s\' for writing\n", fname_string.c_str());
			perror("Err from O/S:");
			exit(EXIT_FAILURE);
		}

		if (verbose_flag)
			printf("Opened %s\n", fname_string.c_str());
	}
	// }}}

	// Give it a header
	// {{{
	fprintf(vmain,
SLASHLINE
"//\n"
"// Filename:\t%sfftmain.v\n"
"// {{{\n" // "}}}"
"// Project:	%s\n"
"//\n"
"// Purpose:	This is the main module in the General Purpose FPGA FFT\n"
"//		implementation.  As such, all other modules are subordinate\n"
"//	to this one.  This module accomplish a fixed size Complex FFT on\n"
"//	%d data points.\n",
		(inverse)?"i":"",prjname, fftsize);
	if (single_clock) {
	fprintf(vmain,
"//	The FFT is fully pipelined, and accepts as inputs one complex two\'s\n"
"//	complement sample per clock.\n");
	} else {
	fprintf(vmain,
"//	The FFT is fully pipelined, and accepts as inputs two complex two\'s\n"
"//	complement samples per clock.\n");
	}

	fprintf(vmain,
"//\n"
"// Parameters:\n"
"//	i_clk\tThe clock.  All operations are synchronous with this clock.\n"
"//	i_%sreset%s\tSynchronous reset, active high.  Setting this line will\n"
"//	\t\tforce the reset of all of the internals to this routine.\n"
"//	\t\tFurther, following a reset, the o_sync line will go\n"
"//	\t\thigh the same time the first output sample is valid.\n",
		(async_reset)?"a":"", (async_reset)?"_n":"");
	if (single_clock) {
		fprintf(vmain,
"//	i_ce\tA clock enable line.  If this line is set, this module\n"
"//	\t\twill accept one complex input value, and produce\n"
"//	\t\tone (possibly empty) complex output value.\n"
"//	i_sample\tThe complex input sample.  This value is split\n"
"//	\t\tinto two two\'s complement numbers, %d bits each, with\n"
"//	\t\tthe real portion in the high order bits, and the\n"
"//	\t\timaginary portion taking the bottom %d bits.\n"
"//	o_result\tThe output result, of the same format as i_sample,\n"
"//	\t\tonly having %d bits for each of the real and imaginary\n"
"//	\t\tcomponents, leading to %d bits total.\n"
"//	o_sync\tA one bit output indicating the first sample of the FFT frame.\n"
"//	\t\tIt also indicates the first valid sample out of the FFT\n"
"//	\t\ton the first frame.\n", nbitsin, nbitsin, nbitsout, nbitsout*2);
	} else {
		fprintf(vmain,
"//	i_ce\tA clock enable line.  If this line is set, this module\n"
"//	\t\twill accept two complex values as inputs, and produce\n"
"//	\t\ttwo (possibly empty) complex values as outputs.\n"
"//	i_left\tThe first of two complex input samples.  This value is split\n"
"//	\t\tinto two two\'s complement numbers, %d bits each, with\n"
"//	\t\tthe real portion in the high order bits, and the\n"
"//	\t\timaginary portion taking the bottom %d bits.\n"
"//	i_right\tThis is the same thing as i_left, only this is the second of\n"
"//	\t\ttwo such samples.  Hence, i_left would contain input\n"
"//	\t\tsample zero, i_right would contain sample one.  On the\n"
"//	\t\tnext clock i_left would contain input sample two,\n"
"//	\t\ti_right number three and so forth.\n"
"//	o_left\tThe first of two output samples, of the same format as i_left,\n"
"//	\t\tonly having %d bits for each of the real and imaginary\n"
"//	\t\tcomponents, leading to %d bits total.\n"
"//	o_right\tThe second of two output samples produced each clock.  This has\n"
"//	\t\tthe same format as o_left.\n"
"//	o_sync\tA one bit output indicating the first valid sample produced by\n"
"//	\t\tthis FFT following a reset.  Ever after, this will\n"
"//	\t\tindicate the first sample of an FFT frame.\n",
	nbitsin, nbitsin, nbitsout, nbitsout*2);
	}

	fprintf(vmain,
"//\n"
"// Arguments:\tThis file was computer generated using the following command\n"
"//\t\tline:\n"
"//\n");
	fprintf(vmain, "//\t\t%% %s\n", cmdline.c_str());
	fprintf(vmain, "//\n");
	fprintf(vmain, "//\tThis core will use hardware accelerated multiplies (DSPs)\n"
		"//\tfor %d of the %d stages\n", mpy_stages, lgval(fftsize));
	fprintf(vmain, "//\n");
	fprintf(vmain, "%s", creator);
	fprintf(vmain, "//\n");
	fprintf(vmain, "%s", cpyleft);
	fprintf(vmain, "//\n//\n`default_nettype\tnone\n//\n");
	// }}}

	// Module declaration
	// {{{
	std::string	resetw("i_reset");
	if (async_reset)
		resetw = "i_areset_n";

	fprintf(vmain, "//\n");
	fprintf(vmain, "//\n");
	fprintf(vmain, "module %sfftmain(i_clk, %s, i_ce,\n",
		(inverse)?"i":"", resetw.c_str());
	if (single_clock) {
		fprintf(vmain, "\t\ti_sample, o_result, o_sync%s);\n",
			(dbg)?", o_dbg":"");
	} else {
		fprintf(vmain, "\t\ti_left, i_right,\n");
		fprintf(vmain, "\t\to_left, o_right, o_sync%s);\n",
			(dbg)?", o_dbg":"");
	}
	fprintf(vmain,
	"\t// The bit-width of the input, IWIDTH, output, OWIDTH, and the log\n"
	"\t// of the FFT size.  These are localparams, rather than parameters,\n"
	"\t// because once the core has been generated, they can no longer be\n"
	"\t// changed.  (These values can be adjusted by running the core\n"
	"\t// generator again.)  The reason is simply that these values have\n"
	"\t// been hardwired into the core at several places.\n");
	fprintf(vmain, "\tlocalparam\tIWIDTH=%d, OWIDTH=%d; // LGWIDTH=%d;\n\t//\n", nbitsin, nbitsout, lgsize);
	assert(lgsize > 0);
	fprintf(vmain, "\tinput\twire\t\t\t\ti_clk, %s, i_ce;\n\t//\n",
		resetw.c_str());
	if (single_clock) {
	fprintf(vmain, "\tinput\twire\t[(2*IWIDTH-1):0]\ti_sample;\n");
	fprintf(vmain, "\toutput\treg\t[(2*OWIDTH-1):0]\to_result;\n");
	} else {
	fprintf(vmain, "\tinput\twire\t[(2*IWIDTH-1):0]\ti_left, i_right;\n");
	fprintf(vmain, "\toutput\treg\t[(2*OWIDTH-1):0]\to_left, o_right;\n");
	}
	fprintf(vmain, "\toutput\treg\t\t\t\to_sync;\n");
	if (dbg)
		fprintf(vmain, "\toutput\twire\t[33:0]\t\to_dbg;\n");
	fprintf(vmain, "\n\n");
	// }}}

	fprintf(vmain, "\t// Outputs of the FFT, ready for bit reversal.\n");
	fprintf(vmain, "\twire\t\t\t\tbr_sync;\n");
	if (single_clock)
		fprintf(vmain, "\twire\t[(2*OWIDTH-1):0]\tbr_result;\n");
	else
		fprintf(vmain, "\twire\t[(2*OWIDTH-1):0]\tbr_left, br_right;\n");

	int	tmp_size = fftsize, lgtmp = lgsize;
	if (fftsize == 2) { // Special case
		// {{{
		if (bitreverse) {
			fprintf(vmain, "\treg\tbr_start;\n");
			fprintf(vmain, "\tinitial br_start = 1\'b0;\n");
			if (async_reset) {
				fprintf(vmain, "\talways @(posedge i_clk, negedge i_areset_n)\n");
				fprintf(vmain, "\tif (!i_areset_n)\n");
			} else {
				fprintf(vmain, "\talways @(posedge i_clk)\n");
				fprintf(vmain, "\tif (i_reset)\n");
			}
			fprintf(vmain, "\t\tbr_start <= 1\'b0;\n");
			fprintf(vmain, "\telse if (i_ce)\n");
			fprintf(vmain, "\t\tbr_start <= 1\'b1;\n");
		}
		fprintf(vmain, "\n\n");
		fprintf(vmain, "\t// verilator lint_off UNUSED\n\twire\t\tw_s2;\n\t// verilator lint_on  UNUSED\n");
		if (single_clock) {
			fprintf(vmain, "\twire\t[%d:0]\tw_d2;\n",
				2*nbitsout-1);
			fprintf(vmain, "\tlaststage\t#(\n"
				"\t\t// {{{\n"
				"\t\t.IWIDTH(IWIDTH),\n"
				"\t\t.OWIDTH(OWIDTH)\n"
				"\t\t// }}}\n"
				"\t) stage_2(\n"
				"\t\t// {{{\n"
				"\t\t.i_clk(i_clk),\n"
				"\t\t.%s(%s),\n"
				"\t\t.i_ce(i_ce),\n",
				resetw.c_str(), resetw.c_str());
			fprintf(vmain, "\t\t.i_sync(%s%s),\n"
					"\t\t.i_val(i_sample),\n"
					"\t\t.o_val(w_d2),\n"
					"\t\t.o_sync(w_s2)\n"
					"\t\t// }}}\n"
					"\t);\n",
				(async_reset)?"":"!", resetw.c_str());
		} else {
			fprintf(vmain, "\twire\t[%d:0]\tw_e2, w_o2;\n",
				2*nbitsout-1);
			fprintf(vmain, "\tlaststage\t#(\n"
				"\t\t// {{{\n"
				"\t\t.IWIDTH(IWIDTH),\n"
				"\t\t.OWIDTH(OWIDTH)\n"
				"\t\t// }}}\n"
				"\t) stage_2(\n"
				"\t\t.i_clk(i_clk),\n"
				"\t\t.%s(%s),\n"
				"\t\t.i_ce(i_ce),\n",
				resetw.c_str(), resetw.c_str());
			fprintf(vmain, "\t\t.i_sync(%s%s),\n"
					"\t\t.i_left(i_left), .i_right(i_right),\n"
					"\t\t.o_left(w_e2), .o_right(w_o2),\n"
					"\t\t.o_sync(w_s2)\n"
					"\t\t// }}}\n"
					"\t);\n",
				(async_reset)?"":"!", resetw.c_str());
		}
		fprintf(vmain, "\n\n");
		// }}}
	} else if (fftsize == 4) { // Special case
		// {{{
		if (!single_clock) {
			fprintf(stderr, "ERR: The two-clocks per sample FFT does not support 4-pt FFTs\n");
			exit(EXIT_FAILURE);
		}

		if (bitreverse) {
			fprintf(vmain, "\treg\tbr_start;\n");
			fprintf(vmain, "\tinitial br_start = 1\'b0;\n");
			if (async_reset) {
				fprintf(vmain, "\talways @(posedge i_clk, negedge i_areset_n)\n");
				fprintf(vmain, "\tif (!i_areset_n)\n");
			} else {
				fprintf(vmain, "\talways @(posedge i_clk)\n");
				fprintf(vmain, "\tif (i_reset)\n");
			}
			fprintf(vmain, "\t\tbr_start <= 1\'b0;\n");
			fprintf(vmain, "\telse if (i_ce)\n");
			fprintf(vmain, "\t\tbr_start <= 1\'b1;\n");
		}
		fprintf(vmain, "\n\n");
		fprintf(vmain, "\twire\t\tw_s4;\n");
		fprintf(vmain, "\twire\t[%d:0]\tw_d4;\n", 2*nbitsout-1);
		fprintf(vmain, "\tqtrstage\t#(\n"
				"\t\t// {{{\n"
				"\t\t.IWIDTH(IWIDTH),\n"
				"\t\t.OWIDTH(OWIDTH)\n"
				"\t\t// }}}\n"
				"\t) stage_4(\n"
				"\t\t// {{{\n"
				"\t\t.i_clk(i_clk),\n"
				"\t\t.%s(%s),\n"
				"\t\t.i_ce(i_ce),\n",
				resetw.c_str(), resetw.c_str());
		fprintf(vmain, "\t\t.i_sync(%s%s),\n"
				"\t\t.i_data(i_sample),\n"
				"\t\t.o_data(w_d4),\n"
				"\t\t.o_sync(w_s4)\n"
				"\t\t// }}}\n"
				"\t);\n",
			(async_reset)?"":"!", resetw.c_str());
		fprintf(vmain, "\n\n");

		fprintf(vmain, "\t// verilator lint_off UNUSED\n\twire\t\tw_s2;\n\t// verilator lint_on  UNUSED\n");
		fprintf(vmain, "\twire\t[%d:0]\tw_d2;\n", 2*nbitsout-1);
		fprintf(vmain, "\tlaststage\t#(\n"
				"\t\t// {{{\n"
				"\t\t.IWIDTH(OWIDTH),\n"
				"\t\t.OWIDTH(OWIDTH)\n"
				"\t\t// }}}\n"
				"\t) stage_2(\n"
				"\t\t// {{{\n"
				"\t\t.i_clk(i_clk),\n"
				"\t\t.%s(%s),\n"
				"\t\t.i_ce(i_ce),\n",
				resetw.c_str(),
				resetw.c_str());
		fprintf(vmain, "\t\t.i_sync(w_s4),\n"
				"\t\t.i_val(w_d4),\n"
				"\t\t.o_val(w_d2),\n"
				"\t\t.o_sync(w_s2)\n"
				"\t\t// }}}\n"
				"\t);\n");
		// }}}
	} else { // General case -- build the FFT stages
		// {{{
		int	nbits = nbitsin, dropbit=0;
		int	obits = nbits+1+xtrapbits;
		std::string	cmem;
		FILE	*cmemfp;

		if ((maxbitsout > 0)&&(obits > maxbitsout))
			obits = maxbitsout;

		// Always do a first stage
		// {{{
		{
			bool	mpystage;

			// Last two stages are always non-multiply stages
			// since the multiplies can be done by adds
			mpystage = ((lgtmp-2) <= mpy_stages);

			fprintf(vmain, "\n\n");
			if (mpystage)
				fprintf(vmain, "\t// A hardware optimized FFT stage\n");
			fprintf(vmain, "\twire\t\tw_s%d;\n", fftsize);
			if (single_clock) {
				// {{{
				fprintf(vmain, "\twire\t[%d:0]\tw_d%d;\n", 2*(obits+xtrapbits)-1, fftsize);
				cmem = gen_coeff_fname(coredir.c_str(), fftsize, 1, 0, inverse);
				cmemfp = gen_coeff_open(cmem.c_str());
				gen_coeffs(cmemfp, fftsize,  nbitsin+xtracbits, 1, 0, inverse);
				cmem = gen_coeff_fname(EMPTYSTR, fftsize, 1, 0, inverse);
				fprintf(vmain, "\tfftstage%s\t#(\n"
					"\t\t// {{{\n"
					"\t\t.IWIDTH(IWIDTH),\n"
					"\t\t.CWIDTH(IWIDTH+%d),\n"
					"\t\t.OWIDTH(%d),\n"
					"\t\t.LGSPAN(%d),\n"
					"\t\t.BFLYSHIFT(0),\n"
					"\t\t.OPT_HWMPY(%d),\n"
					"\t\t.CKPCE(%d),\n"
					"\t\t.COEFFILE(\"%s\")\n"
					"\t\t// }}}\n"
					"\t) stage_%d(\n"
					"\t\t// {{{\n"
					"\t\t.i_clk(i_clk),\n"
					"\t\t.%s(%s),\n"
					"\t\t.i_ce(i_ce),\n",
					((dbg)&&(dbgstage == fftsize))?"_dbg":"",
					xtracbits, obits+xtrapbits,
					lgtmp-1, (mpystage)?1:0,
					ckpce, cmem.c_str(),
					fftsize, resetw.c_str(), resetw.c_str());
				fprintf(vmain, "\t\t.i_sync(%s%s),\n"
					"\t\t.i_data(i_sample),\n"
					"\t\t.o_data(w_d%d),\n"
					"\t\t.o_sync(w_s%d%s)\n"
					"\t\t// }}}\n"
					"\t);\n",
					(async_reset)?"":"!", resetw.c_str(),
					fftsize, fftsize,
					((dbg)&&(dbgstage == fftsize))
						? ", o_dbg":"");
				// }}}
			} else {
				// {{{
				fprintf(vmain, "\t// verilator lint_off UNUSED\n\twire\t\tw_os%d;\n\t// verilator lint_on  UNUSED\n", fftsize);
				fprintf(vmain, "\twire\t[%d:0]\tw_e%d, w_o%d;\n", 2*(obits+xtrapbits)-1, fftsize, fftsize);
				cmem = gen_coeff_fname(coredir.c_str(), fftsize, 2, 0, inverse);
				cmemfp = gen_coeff_open(cmem.c_str());
				gen_coeffs(cmemfp, fftsize,  nbitsin+xtracbits, 2, 0, inverse);
				cmem = gen_coeff_fname(EMPTYSTR, fftsize, 2, 0, inverse);
				fprintf(vmain, "\tfftstage%s\t#(\n"
					"\t\t// {{{\n"
					"\t\t.IWIDTH(IWIDTH),\n"
					"\t\t.CWIDTH(IWIDTH+%d),\n"
					"\t\t.OWIDTH(%d),\n"
					"\t\t.LGSPAN(%d),\n"
					"\t\t.BFLYSHIFT(0),\n"
					"\t\t.OPT_HWMPY(%d),\n"
					"\t\t.CKPCE(%d),\n"
					"\t\t.COEFFILE(\"%s\")\n"
					"\t\t// }}}\n"
					"\t) stage_e%d(\n"
					"\t\t// {{{\n"
					"\t\t.i_clk(i_clk),\n"
					"\t\t.%s(%s),\n"
					"\t\t.i_ce(i_ce),\n",
					((dbg)&&(dbgstage == fftsize))?"_dbg":"",
					xtracbits, obits+xtrapbits,
					lgtmp-2, (mpystage)?1:0,
					ckpce, cmem.c_str(),
					fftsize, resetw.c_str(),
					resetw.c_str());
				fprintf(vmain, "\t\t.i_sync(%s%s),\n"
					"\t\t.i_data(i_left),\n"
					"\t\t.o_data(w_e%d),\n"
					"\t\t.o_sync(w_s%d%s)\n"
					"\t\t// }}}\n"
					"\t);\n",
					(async_reset)?"":"!", resetw.c_str(),
					fftsize, fftsize,
					((dbg)&&(dbgstage == fftsize))?", o_dbg":"");
				cmem = gen_coeff_fname(coredir.c_str(), fftsize, 2, 1, inverse);
				cmemfp = gen_coeff_open(cmem.c_str());
				gen_coeffs(cmemfp, fftsize,  nbitsin+xtracbits, 2, 1, inverse);
				cmem = gen_coeff_fname(EMPTYSTR, fftsize, 2, 1, inverse);
				fprintf(vmain, "\tfftstage\t#(\n"
					"\t\t// {{{\n"
					"\t\t.IWIDTH(IWIDTH),\n"
					"\t\t.CWIDTH(IWIDTH+%d),\n"
					"\t\t.OWIDTH(%d),\n"
					"\t\t.LGSPAN(%d),\n"
					"\t\t.BFLYSHIFT(0),\n"
					"\t\t.OPT_HWMPY(%d),\n"
					"\t\t.CKPCE(%d),\n"
					"\t\t.COEFFILE(\"%s\")\n"
					"\t\t// }}}\n"
					"\t) stage_o%d(\n"
					"\t\t// {{{\n"
					"\t\t.i_clk(i_clk),\n"
					"\t\t.%s(%s),\n"
					"\t\t.i_ce(i_ce),\n",
					xtracbits, obits+xtrapbits,
					lgtmp-2, (mpystage)?1:0,
					ckpce, cmem.c_str(),
					fftsize, resetw.c_str(),
					resetw.c_str());
				fprintf(vmain, "\t\t.i_sync(%s%s),\n"
					"\t\t.i_data(i_right),\n"
					"\t\t.o_data(w_o%d),\n"
					"\t\t.o_sync(w_os%d)\n"
					"\t\t// }}}\n"
					"\t);\n",
					(async_reset)?"":"!",resetw.c_str(),
					fftsize, fftsize);
				// }}}
			}

			// Build the logic for the FFT stage
			// {{{
			std::string	fname;

			fname = coredir + "/";
			if (inverse)
				fname += "i";
			fname += "fftstage";
			if (dbg) {
				std::string	dbgname(fname);
				dbgname += "_dbg";
				dbgname += ".v";
				if (single_clock)
					build_stage(fname.c_str(), fftsize, 1, 0, nbits, xtracbits, ckpce, async_reset, true);
				else
					build_stage(fname.c_str(), fftsize, 2, 1, nbits, xtracbits, ckpce, async_reset, true);
			}

			fname += ".v";
			if (single_clock) {
				build_stage(fname.c_str(), fftsize, 1, 0,
					nbits, xtracbits, ckpce, async_reset,
					false);
			} else {
				// All stages use the same Verilog, so we only
				// need to build one
				build_stage(fname.c_str(), fftsize, 2, 1,
					nbits, xtracbits, ckpce, async_reset, false);
			}
			// }}}
		}
		// }}}

		// Build all following stages, up to the two last ones
		// {{{
		nbits = obits;	// New number of input bits
		tmp_size >>= 1; lgtmp--;
		dropbit = 0;
		fprintf(vmain, "\n\n");
		while(tmp_size >= 8) {
			obits = nbits+((dropbit)?0:1);

			if ((maxbitsout > 0)&&(obits > maxbitsout))
				obits = maxbitsout;

			{
				bool		mpystage;

				mpystage = ((lgtmp-2) <= mpy_stages);

				if (mpystage)
					fprintf(vmain, "\t// A hardware optimized FFT stage\n");
				fprintf(vmain, "\twire\t\tw_s%d;\n",
					tmp_size);
				if (single_clock) {
					// {{{
					fprintf(vmain,"\twire\t[%d:0]\tw_d%d;\n",
						2*(obits+xtrapbits)-1,
						tmp_size);
					cmem = gen_coeff_fname(coredir.c_str(), tmp_size, 1, 0, inverse);
					cmemfp = gen_coeff_open(cmem.c_str());
					gen_coeffs(cmemfp, tmp_size,
						nbits+xtracbits+xtrapbits, 1, 0, inverse);
					cmem = gen_coeff_fname(EMPTYSTR, tmp_size, 1, 0, inverse);
					fprintf(vmain, "\tfftstage%s\t#(\n"
						"\t\t// {{{\n"
						"\t\t.IWIDTH(%d),\n"
						"\t\t.CWIDTH(%d),\n"
						"\t\t.OWIDTH(%d),\n"
						"\t\t.LGSPAN(%d),\n"
						"\t\t.BFLYSHIFT(%d),\n"
						"\t\t.OPT_HWMPY(%d),\n"
						"\t\t.CKPCE(%d),\n"
						"\t\t.COEFFILE(\"%s\")\n"
						"\t\t// }}}\n"
						"\t) stage_%d(\n"
						"\t\t// {{{\n"
						"\t\t.i_clk(i_clk),\n"
						"\t\t.%s(%s),\n"
						"\t\t.i_ce(i_ce),\n",
						((dbg)&&(dbgstage==tmp_size))?"_dbg":"",
						nbits+xtrapbits,
						nbits+xtracbits+xtrapbits,
						obits+xtrapbits,
						lgtmp-1, (dropbit)?0:0, (mpystage)?1:0,
						ckpce,
						cmem.c_str(), tmp_size,
						resetw.c_str(),
						resetw.c_str());
					fprintf(vmain, "\t\t.i_sync(w_s%d),\n"
						"\t\t.i_data(w_d%d),\n"
						"\t\t.o_data(w_d%d),\n"
						"\t\t.o_sync(w_s%d%s)\n"
						"\t\t// }}}\n"
						"\t);\n",
						tmp_size<<1, tmp_size<<1,
						tmp_size, tmp_size,
						((dbg)&&(dbgstage == tmp_size))
							?", o_dbg":"");
					// }}}
				} else {
					// {{{
					fprintf(vmain, "\t// verilator lint_off UNUSED\n\twire\t\tw_os%d;\n\t// verilator lint_on  UNUSED\n",
						tmp_size);
					fprintf(vmain,"\twire\t[%d:0]\tw_e%d, w_o%d;\n",
						2*(obits+xtrapbits)-1,
						tmp_size, tmp_size);
					cmem = gen_coeff_fname(coredir.c_str(), tmp_size, 2, 0, inverse);
					cmemfp = gen_coeff_open(cmem.c_str());
					gen_coeffs(cmemfp, tmp_size,
						nbits+xtracbits+xtrapbits, 2, 0, inverse);
					cmem = gen_coeff_fname(EMPTYSTR, tmp_size, 2, 0, inverse);
					fprintf(vmain, "\tfftstage%s\t#(\n"
						"\t\t// {{{\n"
						"\t\t.IWIDTH(%d),\n"
						"\t\t.CWIDTH(%d),\n"
						"\t\t.OWIDTH(%d),\n"
						"\t\t.LGSPAN(%d),\n"
						"\t\t.BFLYSHIFT(%d),\n"
						"\t\t.OPT_HWMPY(%d),\n"
						"\t\t.CKPCE(%d),\n"
						"\t\t.COEFFILE(\"%s\")\n"
						"\t\t// }}}\n"
						"\t) stage_e%d(\n"
						"\t\t// {{{\n"
						"\t\t.i_clk(i_clk),\n"
						"\t\t.%s(%s),\n"
						"\t\t.i_ce(i_ce),\n",
						((dbg)&&(dbgstage==tmp_size))?"_dbg":"",
						nbits+xtrapbits,
						nbits+xtracbits+xtrapbits,
						obits+xtrapbits,
						lgtmp-2, (dropbit)?0:0, (mpystage)?1:0,
						ckpce,
						cmem.c_str(), tmp_size,
						resetw.c_str(),
						resetw.c_str());
					fprintf(vmain, "\t\t.i_sync(w_s%d),\n"
						"\t\t.i_data(w_e%d),\n"
						"\t\t.o_data(w_e%d),\n"
						"\t\t.o_sync(w_s%d%s)\n"
						"\t\t// }}}\n"
						"\t);\n",
						tmp_size<<1, tmp_size<<1,
						tmp_size, tmp_size,
						((dbg)&&(dbgstage == tmp_size))
							?", o_dbg":"");
					cmem = gen_coeff_fname(coredir.c_str(),
						tmp_size, 2, 1, inverse);
					cmemfp = gen_coeff_open(cmem.c_str());
					gen_coeffs(cmemfp, tmp_size,
						nbits+xtracbits+xtrapbits,
						2, 1, inverse);
					cmem = gen_coeff_fname(EMPTYSTR,
						tmp_size, 2, 1, inverse);
					fprintf(vmain, "\tfftstage\t#(\n"
						"\t\t// {{{\n"
						"\t\t.IWIDTH(%d),\n"
						"\t\t.CWIDTH(%d),\n"
						"\t\t.OWIDTH(%d),\n"
						"\t\t.LGSPAN(%d),\n"
						"\t\t.BFLYSHIFT(%d),\n"
						"\t\t.OPT_HWMPY(%d),\n"
						"\t\t.CKPCE(%d),\n"
						"\t\t.COEFFILE(\"%s\")\n"
						"\t\t// }}}\n"
						"\n) \tstage_o%d(\n"
						"\t\t// {{{\n"
						"\t\t.i_clk(i_clk),\n"
						"\t\t.%s(%s),\n"
						"\t\t.i_ce(i_ce),\n",
						nbits+xtrapbits,
						nbits+xtracbits+xtrapbits,
						obits+xtrapbits,
						lgtmp-2, (dropbit)?0:0, (mpystage)?1:0,
						ckpce, cmem.c_str(), tmp_size,
						resetw.c_str(),
						resetw.c_str());
					fprintf(vmain, "\t\t.i_sync(w_s%d),\n"
						"\t\t.i_data(w_o%d),\n"
						"\t\t.o_data(w_o%d),\n"
						"\t\t.o_sync(w_os%d)\n"
						"\t\t// }}}\n"
						"\t);\n",
						tmp_size<<1, tmp_size<<1,
						tmp_size, tmp_size);
					// }}}
				}
				fprintf(vmain, "\n");
			}


			dropbit ^= 1;
			nbits = obits;
			tmp_size >>= 1; lgtmp--;
		}
		// }}}

		// The Quarter stage : 90 degrees, adds and subtracts only
		// {{{
		if (tmp_size == 4) {
			obits = nbits+((dropbit)?0:1);

			if ((maxbitsout > 0)&&(obits > maxbitsout))
				obits = maxbitsout;

			fprintf(vmain, "\twire\t\tw_s4;\n");
			if (single_clock) {
				// {{{
				fprintf(vmain, "\twire\t[%d:0]\tw_d4;\n",
					2*(obits+xtrapbits)-1);
				fprintf(vmain, "\tqtrstage%s\t#(\n"
					"\t\t// {{{\n"
					"\t\t.IWIDTH(%d),\n"
					"\t\t.OWIDTH(%d),\n"
					"\t\t.LGWIDTH(%d),\n"
					"\t\t.INVERSE(%d),\n"
					"\t\t.SHIFT(%d)\n"
					"\t\t// }}}\n"
					"\t) stage_4(\n"
						"\t\t// {{{\n"
						"\t\t.i_clk(i_clk),\n"
						"\t\t.%s(%s),\n"
						"\t\t.i_ce(i_ce),\n",
					((dbg)&&(dbgstage==4))?"_dbg":"",
					nbits+xtrapbits, obits+xtrapbits, lgsize,
					(inverse)?1:0, (dropbit)?0:0,
					resetw.c_str(),
					resetw.c_str());
				fprintf(vmain, "\t\t.i_sync(w_s8),\n"
					"\t\t.i_data(w_d8),\n"
					"\t\t.o_data(w_d4),\n"
					"\t\t.o_sync(w_s4%s)\n"
					"\t\t// }}}\n"
					"\t);\n",
					((dbg)&&(dbgstage==4))?", o_dbg":"");
				// }}}
			} else {
				// {{{
				fprintf(vmain, "\t// verilator lint_off UNUSED\n\twire\t\tw_os4;\n\t// verilator lint_on  UNUSED\n");
				fprintf(vmain, "\twire\t[%d:0]\tw_e4, w_o4;\n", 2*(obits+xtrapbits)-1);
				fprintf(vmain, "\tqtrstage%s\t#(\n"
					"\t\t// {{{\n"
					"\t\t.IWIDTH(%d),\n"
					"\t\t.OWIDTH(%d),\n"
					"\t\t.LGWIDTH(%d),\n"
					"\t\t.ODD(0),\n"
					"\t\t.INVERSE(%d),\n"
					"\t\t.SHIFT(%d)\n"
					"\t\t// }}}\n"
					"\t) stage_e4(\n"
					"\t\t// {{{\n"
					"\t\t.i_clk(i_clk),\n"
					"\t\t.%s(%s),\n"
					"\t\t.i_ce(i_ce),\n",
					((dbg)&&(dbgstage==4))?"_dbg":"",
					nbits+xtrapbits, obits+xtrapbits, lgsize,
					(inverse)?1:0, (dropbit)?0:0,
					resetw.c_str(),
					resetw.c_str());
				fprintf(vmain, "\t\t.i_sync(w_s8),\n"
					"\t\t.i_data(w_e8),\n"
					"\t\t.o_data(w_e4),\n"
					"\t\t.o_sync(w_s4%s)\n"
					"\t\t// }}}\n"
					"\t);\n",
					((dbg)&&(dbgstage==4))?", o_dbg":"");
				fprintf(vmain, "\tqtrstage\t#(\n"
					"\t\t// {{{\n"
					"\t\t.IWIDTH(%d),\n"
					"\t\t.OWIDTH(%d),\n"
					"\t\t.LGWIDTH(%d),\n"
					"\t\t.ODD(1),\n"
					"\t\t.INVERSE(%d),\n"
					"\t\t.SHIFT(%d)\n"
					"\t\t// }}}\n"
					"\t) stage_o4(\n"
					"\t\t// {{{\n"
					"\t\t.i_clk(i_clk),\n"
					"\t\t.%s(%s),\n"
					"\t\t.i_ce(i_ce),\n",
					nbits+xtrapbits, obits+xtrapbits, lgsize, (inverse)?1:0, (dropbit)?0:0,
					resetw.c_str(),
					resetw.c_str());
				fprintf(vmain, "\t\t.i_sync(w_s8),\n"
					"\t\t.i_data(w_o8),\n"
					"\t\t.o_data(w_o4),\n"
					"\t\t.o_sync(w_os4)\n"
					"\t\t// }}}\n"
					"\t);\n");
				// }}}
			}
			dropbit ^= 1;
			nbits = obits;
			tmp_size >>= 1; lgtmp--;
		}
		// }}}

		// The last stage : adds and subtracts only
		// {{{
		{
			obits = nbits+((dropbit)?0:1);
			if (obits > nbitsout)
				obits = nbitsout;
			if ((maxbitsout>0)&&(obits > maxbitsout))
				obits = maxbitsout;
			fprintf(vmain, "\t// verilator lint_off UNUSED\n\twire\t\tw_s2;\n\t// verilator lint_on  UNUSED\n");
			if (single_clock) {
				fprintf(vmain, "\twire\t[%d:0]\tw_d2;\n",
					2*obits-1);
			} else {
				fprintf(vmain, "\twire\t[%d:0]\tw_e2, w_o2;\n",
					2*obits-1);
			}
			/*
			if ((nbits+xtrapbits+1 == obits)&&(!dropbit))
				printf("Warning: Less than optimal scaling\n");
			*/

			if (single_clock) {
				// {{{
				fprintf(vmain, "\tlaststage\t#(\n"
					"\t\t// {{{\n"
					"\t\t.IWIDTH(%d),\n"
					"\t\t.OWIDTH(%d),\n"
					"\t\t.SHIFT(%d)\n"
					"\t\t// }}}\n"
					"\t) stage_2(\n"
					"\t\t// {{{\n"
					"\t\t.i_clk(i_clk),\n"
					"\t\t.%s(%s),\n"
					"\t\t.i_ce(i_ce),\n",
					nbits+xtrapbits, obits,(dropbit)?0:1,
					resetw.c_str(), resetw.c_str());
				fprintf(vmain, "\t\t.i_sync(w_s4),\n"
						"\t\t.i_val(w_d4),\n"
						"\t\t.o_val(w_d2),\n"
						"\t\t.o_sync(w_s2)\n"
						"\t\t// }}}\n"
						"\t);\n");
				// }}}
			} else {
				// {{{
				fprintf(vmain, "\tlaststage\t#(\n"
					"\t\t// {{{\n"
					"\t\t.IWIDTH(%d),\n"
					"\t\t.OWIDTH(%d),\n"
					"\t\t.SHIFT(%d)\n"
					"\t\t// }}}\n"
					"\t) stage_2(\n"
					"\t\t// {{{\n"
					"\t\t.i_clk(i_clk),\n"
					"\t\t.%s(%s),\n"
					"\t\t.i_ce(i_ce),\n",
					nbits+xtrapbits, obits,(dropbit)?0:1,
					resetw.c_str(), resetw.c_str());
				fprintf(vmain, "\t\t.i_sync(w_s4),\n"
					"\t\t.i_left(w_e4), .i_right(w_o4),\n"
					"\t\t.o_left(w_e2), .o_right(w_o2),\n"
					"\t\t.o_sync(w_s2)\n"
					"\t\t// }}}\n"
					"\t);\n");
				// }}}
			}

			fprintf(vmain, "\n\n");
			nbits = obits;
		}
		// }}}

		if (bitreverse) {	// Prep for bit reversal
			// {{{
			fprintf(vmain, "\twire\tbr_start;\n");
			fprintf(vmain, "\treg\tr_br_started;\n");
			fprintf(vmain, "\tinitial\tr_br_started = 1\'b0;\n");
			if (async_reset) {
				fprintf(vmain, "\talways @(posedge i_clk, negedge i_areset_n)\n");
				fprintf(vmain, "\tif (!i_areset_n)\n");
			} else {
				fprintf(vmain, "\talways @(posedge i_clk)\n");
				fprintf(vmain, "\tif (i_reset)\n");
			}
			fprintf(vmain, "\t\tr_br_started <= 1\'b0;\n");
			fprintf(vmain, "\telse if (i_ce)\n");
			fprintf(vmain, "\t\tr_br_started <= r_br_started || w_s2;\n");
			fprintf(vmain, "\tassign\tbr_start = r_br_started || w_s2;\n");
			// }}}
		}
		// }}}
	}

	// Bit-reversal stage
	// {{{
	fprintf(vmain, "\n");
	fprintf(vmain, "\t// Now for the bit-reversal stage.\n");
	if (bitreverse) {
		if (single_clock) {
			fprintf(vmain, "\tbitreverse\t#(\n"
				"\t\t// {{{\n"
				"\t\t.LGSIZE(%d), .WIDTH(%d)\n"
				"\t\t// }}}\n"
				"\t) revstage (\n"
				"\t\t// {{{\n"
				"\t\t.i_clk(i_clk),\n"
				"\t\t.%s(%s),\n", lgsize, nbitsout,
				resetw.c_str(),
				resetw.c_str());
			fprintf(vmain,
				"\t\t.i_ce(i_ce & br_start),\n"
				"\t\t.i_in(w_d2),\n"
				"\t\t.o_out(br_result),\n"
				"\t\t.o_sync(br_sync)\n"
				"\t\t// }}}\n"
				"\t);\n");
		} else {
			fprintf(vmain, "\tbitreverse\t#(\n"
				"\t\t// {{{\n"
				"\t\t.LGSIZE(%d),\n"
				"\t\t.WIDTH(%d)\n"
				"\t\t// }}}\n"
				"\t) revstage (\n"
				"\t\t// {{{\n"
				"\t\t.i_clk(i_clk),\n"
				"\t\t.%s(%s),\n", lgsize, nbitsout,
				resetw.c_str(), resetw.c_str());
			fprintf(vmain,
				"\t\t.i_ce(i_ce & br_start),\n"
				"\t\t.i_in_0(w_e2),\n"
				"\t\t.i_in_1(w_o2),\n"
				"\t\t.o_out_0(br_left),\n"
				"\t\t.o_out_1(br_right),\n"
				"\t\t.o_sync(br_sync)\n"
				"\t\t// }}}\n"
				"\t);\n");
			// fprintf(vmain, "\t\t\t(i_ce & br_start), w_e2, w_o2,\n");
			// fprintf(vmain, "\t\t\tbr_left, br_right, br_sync);\n");
		}
	} else {
		fprintf(vmain, "\t//\n"
"\t// Since the bit-reversal stage isn\'t included, according to the current\n"
"\t// settings, this will just be a stub instead of the actual bit-reversal\n"
"\t// logic.\n"
"\t//\n");
		if (single_clock) {
			fprintf(vmain, "\tassign\tbr_result   = w_d2;\n");
		} else {
			fprintf(vmain, "\tassign\tbr_left  = w_e2;\n");
			fprintf(vmain, "\tassign\tbr_right = w_o2;\n");
		}
		fprintf(vmain, "\tassign\tbr_sync    = w_s2;\n");
	}
	// }}}

	// Register the final outputs and we're done
	// {{{
	fprintf(vmain,
"\n\n"
"\t// Last clock: Register our outputs, we\'re done.\n"
"\tinitial\to_sync  = 1\'b0;\n");
	if (async_reset)
		fprintf(vmain,
"\talways @(posedge i_clk, negedge i_areset_n)\n\tif (!i_areset_n)\n");
	else {
		fprintf(vmain,
"\talways @(posedge i_clk)\n\tif (i_reset)\n");
	}

	fprintf(vmain,
"\t\to_sync  <= 1\'b0;\n"
"\telse if (i_ce)\n"
"\t\to_sync  <= br_sync;\n"
"\n"
"\talways @(posedge i_clk)\n"
"\tif (i_ce)\n");
	if (single_clock) {
		fprintf(vmain, "\t\to_result  <= br_result;\n");
	} else {
		fprintf(vmain,
"\tbegin\n"
"\t\to_left  <= br_left;\n"
"\t\to_right <= br_right;\n"
"\tend\n");
	}

	fprintf(vmain,
"\n\n"
"endmodule\n");
	fclose(vmain);
	// }}}
	// }}}
	////////////////////////////////////////////////////////////////////////
	//
	// Build the component modules
	// {{{
	////////////////////////////////////////////////////////////////////////
	//
	//
	{
		std::string	fname;

		// Butterfly
		// {{{
		fname = coredir + "/butterfly.v";
		build_butterfly(fname.c_str(), xtracbits, rounding,
			ckpce, async_reset);
		// }}}

		// The hardware assisted butterfly
		// {{{
		fname = coredir + "/hwbfly.v";
		build_hwbfly(fname.c_str(), xtracbits, rounding,
			ckpce, async_reset);
		// }}}

		// The binary multiply the hardware assisted multiply depends on
		// {{{
		{
			fname = coredir + "/longbimpy.v";
			build_longbimpy(fname.c_str());
			fname = coredir + "/bimpy.v";
			build_bimpy(fname.c_str());
		}
		// }}}

		// Quarter stage
		// {{{
		if ((dbg)&&(dbgstage == 4)) {
			fname = coredir + "/qtrstage_dbg.v";
			if (single_clock)
				build_snglquarters(fname.c_str(), rounding,
					async_reset, true);
			else
				build_dblquarters(fname.c_str(), rounding,
					async_reset, true);
		}
		fname = coredir + "/qtrstage.v";
		if (single_clock)
			build_snglquarters(fname.c_str(), rounding,
					async_reset, false);
		else
			build_dblquarters(fname.c_str(), rounding,
					async_reset, false);
		// }}}

		// Last stage
		// {{{
		if (single_clock) {
			fname = coredir + "/laststage.v";
			build_sngllast(fname.c_str(), async_reset);
		} else {
			if ((dbg)&&(dbgstage == 2))
				fname = coredir + "/laststage_dbg.v";
			else
				fname = coredir + "/laststage.v";
			build_dblstage(fname.c_str(), rounding,
				async_reset, (dbg)&&(dbgstage==2));
		}
		// }}}

		// Bit reversal logic
		// {{{
		if (bitreverse) {
			fname = coredir + "/bitreverse.v";
			if (single_clock)
				build_snglbrev(fname.c_str(), async_reset);
			else
				build_dblreverse(fname.c_str(), async_reset);
		}
		// }}}

		// Rounding
		// {{{
		const	char	*rnd_string = "";
		switch(rounding) {
			case RND_TRUNCATE:	rnd_string = "/truncate.v"; break;
			case RND_FROMZERO:	rnd_string = "/roundfromzero.v"; break;
			case RND_HALFUP:	rnd_string = "/roundhalfup.v"; break;
			default:
				rnd_string = "/convround.v"; break;
		} fname = coredir + rnd_string;
		switch(rounding) {
			case RND_TRUNCATE: build_truncator(fname.c_str()); break;
			case RND_FROMZERO: build_roundfromzero(fname.c_str()); break;
			case RND_HALFUP: build_roundhalfup(fname.c_str()); break;
			default:
				build_convround(fname.c_str()); break;
		}
		// }}}
	}
	// }}}

	if (verbose_flag)
		printf("All done -- success\n");
}
