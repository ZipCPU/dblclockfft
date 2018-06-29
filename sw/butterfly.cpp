////////////////////////////////////////////////////////////////////////////////
//
// Filename: 	butterfly.cpp
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

#else
// And for G++/Linux environment

#include <unistd.h>	// Defines the R_OK/W_OK/etc. macros
#include <sys/stat.h>
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

void	build_butterfly(const char *fname, int xtracbits, ROUND_T rounding,
			int	ckpce, const bool async_reset) {
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

	//if (ckpce >= 3)
		//ckpce = 3;
	if (ckpce <= 1)
		ckpce = 1;
	if (ckpce > 1) {
		fprintf(stderr, "WARNING: Butterfly code does not yet support CKPCE=%d\n", ckpce);
		fprintf(stderr, "WARNING: Using CKPCE=1 instead\n");
		ckpce = 1;
	}

	std::string	resetw("i_reset");
	if (async_reset)
		resetw = std::string("i_areset_n");


	fprintf(fp,
SLASHLINE
"//\n"
"// Filename:\tbutterfly.v\n"
"//\n"
"// Project:\t%s\n"
"//\n"
"// Purpose:\tThis routine caculates a butterfly for a decimation\n"
"//		in frequency version of an FFT.  Specifically, given\n"
"//	complex Left and Right values together with a coefficient, the output\n"
"//	of this routine is given by:\n"
"//\n"
"//		L' = L + R\n"
"//		R' = (L - R)*C\n"
"//\n"
"//	The rest of the junk below handles timing (mostly), to make certain\n"
"//	that L' and R' reach the output at the same clock.  Further, just to\n"
"//	make certain that is the case, an 'aux' input exists.  This aux value\n"
"//	will come out of this routine synchronized to the values it came in\n"
"//	with.  (i.e., both L', R', and aux all have the same delay.)  Hence,\n"
"//	a caller of this routine may set aux on the first input with valid\n"
"//	data, and then wait to see aux set on the output to know when to find\n"
"//	the first output with valid data.\n"
"//\n"
"//	All bits are preserved until the very last clock, where any more bits\n"
"//	than OWIDTH will be quietly discarded.\n"
"//\n"
"//	This design features no overflow checking.\n"
"//\n"
"// Notes:\n"
"//	CORDIC:\n"
"//		Much as we might like, we can't use a cordic here.\n"
"//		The goal is to accomplish an FFT, as defined, and a\n"
"//		CORDIC places a scale factor onto the data.  Removing\n"
"//		the scale factor would cost two multiplies, which\n"
"//		is precisely what we are trying to avoid.\n"
"//\n"
"//\n"
"//	3-MULTIPLIES:\n"
"//		It should also be possible to do this with three multiplies\n"
"//		and an extra two addition cycles.\n"
"//\n"
"//		We want\n"
"//			R+I = (a + jb) * (c + jd)\n"
"//			R+I = (ac-bd) + j(ad+bc)\n"
"//		We multiply\n"
"//			P1 = ac\n"
"//			P2 = bd\n"
"//			P3 = (a+b)(c+d)\n"
"//		Then\n"
"//			R+I=(P1-P2)+j(P3-P2-P1)\n"
"//\n"
"//		WIDTHS:\n"
"//		On multiplying an X width number by an\n"
"//		Y width number, X>Y, the result should be (X+Y)\n"
"//		bits, right?\n"
"//		-2^(X-1) <= a <= 2^(X-1) - 1\n"
"//		-2^(Y-1) <= b <= 2^(Y-1) - 1\n"
"//		(2^(Y-1)-1)*(-2^(X-1)) <= ab <= 2^(X-1)2^(Y-1)\n"
"//		-2^(X+Y-2)+2^(X-1) <= ab <= 2^(X+Y-2) <= 2^(X+Y-1) - 1\n"
"//		-2^(X+Y-1) <= ab <= 2^(X+Y-1)-1\n"
"//		YUP!  But just barely.  Do this and you'll really want\n"
"//		to drop a bit, although you will risk overflow in so\n"
"//		doing.\n"
"//\n"
"//	20150602 -- The sync logic lines have been completely redone.  The\n"
"//		synchronization lines no longer go through the FIFO with the\n"
"//		left hand sum, but are kept out of memory.  This allows the\n"
"//		butterfly to use more optimal memory resources, while also\n"
"//		guaranteeing that the sync lines can be properly reset upon\n"
"//		any reset signal.\n"
"//\n"
"//\n%s"
"//\n", prjname, creator);
	fprintf(fp, "%s", cpyleft);
	fprintf(fp, "//\n//\n`default_nettype\tnone\n//\n");

	fprintf(fp,
"module\tbutterfly(i_clk, %s, i_ce, i_coef, i_left, i_right, i_aux,\n"
		"\t\to_left, o_right, o_aux);\n"
	"\t// Public changeable parameters ...\n"
	"\tparameter IWIDTH=%d,", resetw.c_str(), TST_BUTTERFLY_IWIDTH);
#ifdef	TST_BUTTERFLY_CWIDTH
	fprintf(fp, "CWIDTH=%d,", TST_BUTTERFLY_CWIDTH);
#else
	fprintf(fp, "CWIDTH=IWIDTH+%d,", xtracbits);
#endif
#ifdef	TST_BUTTERFLY_OWIDTH
	fprintf(fp, "OWIDTH=%d;\n", TST_BUTTERFLY_OWIDTH);
#else
	fprintf(fp, "OWIDTH=IWIDTH+1;\n");
#endif
	fprintf(fp,
	"\t// Parameters specific to the core that should not be changed.\n"
	"\tparameter	MPYDELAY=%d'd%d,\n"
			"\t\t\tSHIFT=0, AUXLEN=(MPYDELAY+3);\n"
	"\t// The LGDELAY should be the base two log of the MPYDELAY.  If\n"
	"\t// this value is fractional, then round up to the nearest\n"
	"\t// integer: LGDELAY=ceil(log(MPYDELAY)/log(2));\n"
	"\tparameter\tLGDELAY=%d;\n"
	"\tparameter\tCKPCE=%d;\n"
	"\tinput\t\ti_clk, %s, i_ce;\n"
	"\tinput\t\t[(2*CWIDTH-1):0] i_coef;\n"
	"\tinput\t\t[(2*IWIDTH-1):0] i_left, i_right;\n"
	"\tinput\t\ti_aux;\n"
	"\toutput\twire	[(2*OWIDTH-1):0] o_left, o_right;\n"
	"\toutput\treg\to_aux;\n"
	"\n", lgdelay(16,xtracbits), bflydelay(16, xtracbits),
		lgdelay(16,xtracbits), ckpce, resetw.c_str());
	fprintf(fp,
	"\treg\t[(2*IWIDTH-1):0]\tr_left, r_right;\n"
	"\treg\t[(2*CWIDTH-1):0]\tr_coef, r_coef_2;\n"
	"\twire\tsigned\t[(IWIDTH-1):0]\tr_left_r, r_left_i, r_right_r, r_right_i;\n"
	"\tassign\tr_left_r  = r_left[ (2*IWIDTH-1):(IWIDTH)];\n"
	"\tassign\tr_left_i  = r_left[ (IWIDTH-1):0];\n"
	"\tassign\tr_right_r = r_right[(2*IWIDTH-1):(IWIDTH)];\n"
	"\tassign\tr_right_i = r_right[(IWIDTH-1):0];\n"
"\n"
	"\treg\tsigned\t[(IWIDTH):0]\tr_sum_r, r_sum_i, r_dif_r, r_dif_i;\n"
"\n"
	"\treg	[(LGDELAY-1):0]	fifo_addr;\n"
	"\twire	[(LGDELAY-1):0]	fifo_read_addr;\n"
	"\tassign\tfifo_read_addr = fifo_addr - MPYDELAY[(LGDELAY-1):0];\n"
	"\treg	[(2*IWIDTH+1):0]	fifo_left [ 0:((1<<LGDELAY)-1)];\n"
"\n");
	fprintf(fp,
	"\t// Set up the input to the multiply\n"
	"\talways @(posedge i_clk)\n"
		"\t\tif (i_ce)\n"
		"\t\tbegin\n"
			"\t\t\t// One clock just latches the inputs\n"
			"\t\t\tr_left <= i_left;	// No change in # of bits\n"
			"\t\t\tr_right <= i_right;\n"
			"\t\t\tr_coef  <= i_coef;\n"
			"\t\t\t// Next clock adds/subtracts\n"
			"\t\t\tr_sum_r <= r_left_r + r_right_r; // Now IWIDTH+1 bits\n"
			"\t\t\tr_sum_i <= r_left_i + r_right_i;\n"
			"\t\t\tr_dif_r <= r_left_r - r_right_r;\n"
			"\t\t\tr_dif_i <= r_left_i - r_right_i;\n"
			"\t\t\t// Other inputs are simply delayed on second clock\n"
			"\t\t\tr_coef_2<= r_coef;\n"
	"\t\tend\n"
"\n");
	fprintf(fp,
	"\t// Don\'t forget to record the even side, since it doesn\'t need\n"
	"\t// to be multiplied, but yet we still need the results in sync\n"
	"\t// with the answer when it is ready.\n"
	"\tinitial fifo_addr = 0;\n");
	if (async_reset)
		fprintf(fp, "\talways @(posedge i_clk, negedge i_areset_n)\n\t\tif (!i_areset_n)\n");
	else
		fprintf(fp, "\talways @(posedge i_clk)\n\t\tif (i_reset)\n");
	fprintf(fp,
			"\t\t\tfifo_addr <= 0;\n"
		"\t\telse if (i_ce)\n"
			"\t\t\t// Need to delay the sum side--nothing else happens\n"
			"\t\t\t// to it, but it needs to stay synchronized with the\n"
			"\t\t\t// right side.\n"
			"\t\t\tfifo_addr <= fifo_addr + 1;\n"
"\n"
	"\talways @(posedge i_clk)\n"
		"\t\tif (i_ce)\n"
			"\t\t\tfifo_left[fifo_addr] <= { r_sum_r, r_sum_i };\n"
"\n"
	"\twire\tsigned\t[(CWIDTH-1):0]	ir_coef_r, ir_coef_i;\n"
	"\tassign\tir_coef_r = r_coef_2[(2*CWIDTH-1):CWIDTH];\n"
	"\tassign\tir_coef_i = r_coef_2[(CWIDTH-1):0];\n"
	"\twire\tsigned\t[((IWIDTH+2)+(CWIDTH+1)-1):0]\tp_one, p_two, p_three;\n"
"\n"
"\n");
	fprintf(fp,
	"\t// Multiply output is always a width of the sum of the widths of\n"
	"\t// the two inputs.  ALWAYS.  This is independent of the number of\n"
	"\t// bits in p_one, p_two, or p_three.  These values needed to\n"
	"\t// accumulate a bit (or two) each.  However, this approach to a\n"
	"\t// three multiply complex multiply cannot increase the total\n"
	"\t// number of bits in our final output.  We\'ll take care of\n"
	"\t// dropping back down to the proper width, OWIDTH, in our routine\n"
	"\t// below.\n"
"\n"
"\n");
	fprintf(fp,
	"\t// We accomplish here \"Karatsuba\" multiplication.  That is,\n"
	"\t// by doing three multiplies we accomplish the work of four.\n"
	"\t// Let\'s prove to ourselves that this works ... We wish to\n"
	"\t// multiply: (a+jb) * (c+jd), where a+jb is given by\n"
	"\t//\ta + jb = r_dif_r + j r_dif_i, and\n"
	"\t//\tc + jd = ir_coef_r + j ir_coef_i.\n"
	"\t// We do this by calculating the intermediate products P1, P2,\n"
	"\t// and P3 as\n"
	"\t//\tP1 = ac\n"
	"\t//\tP2 = bd\n"
	"\t//\tP3 = (a + b) * (c + d)\n"
	"\t// and then complete our final answer with\n"
	"\t//\tac - bd = P1 - P2 (this checks)\n"
	"\t//\tad + bc = P3 - P2 - P1\n"
	"\t//\t        = (ac + bc + ad + bd) - bd - ac\n"
	"\t//\t        = bc + ad (this checks)\n"
"\n"
"\n");
	fprintf(fp,
	"\t// This should really be based upon an IF, such as in\n"
	"\t// if (IWIDTH < CWIDTH) then ...\n"
	"\t// However, this is the only (other) way I know to do it.\n"
	"\tgenerate if (CKPCE <= 1)\n"
	"\tbegin\n"
"\n"
		"\t\twire\t[(CWIDTH):0]\tp3c_in;\n"
		"\t\twire\t[(IWIDTH+1):0]\tp3d_in;\n"
		"\t\tassign\tp3c_in = ir_coef_i + ir_coef_r;\n"
		"\t\tassign\tp3d_in = r_dif_r + r_dif_i;\n"
		"\n"
		"\t\t// We need to pad these first two multiplies by an extra\n"
		"\t\t// bit just to keep them aligned with the third,\n"
		"\t\t// simpler, multiply.\n"
		"\t\t%s #(CWIDTH+1,IWIDTH+2) p1(i_clk, i_ce,\n"
				"\t\t\t\t{ir_coef_r[CWIDTH-1],ir_coef_r},\n"
				"\t\t\t\t{r_dif_r[IWIDTH],r_dif_r}, p_one);\n"
		"\t\t%s #(CWIDTH+1,IWIDTH+2) p2(i_clk, i_ce,\n"
				"\t\t\t\t{ir_coef_i[CWIDTH-1],ir_coef_i},\n"
				"\t\t\t\t{r_dif_i[IWIDTH],r_dif_i}, p_two);\n"
		"\t\t%s #(CWIDTH+1,IWIDTH+2) p3(i_clk, i_ce,\n"
			"\t\t\t\tp3c_in, p3d_in, p_three);\n"
"\n",
		(USE_OLD_MULTIPLY)?"shiftaddmpy":"longbimpy",
		(USE_OLD_MULTIPLY)?"shiftaddmpy":"longbimpy",
		(USE_OLD_MULTIPLY)?"shiftaddmpy":"longbimpy");

	///////////////////////////////////////////
	///
	///	Two clocks per CE, so CE, no-ce, CE, no-ce, etc
	///
	fprintf(fp,
	"\tend else if (CKPCE == 2)\n"
	"\tbegin : CKPCE_TWO\n"
		"\t\t// Coefficient multiply inputs\n"
		"\t\treg		[2*(CWIDTH)-1:0]	mpy_pipe_c;\n"
		"\t\t// Data multiply inputs\n"
		"\t\treg		[2*(IWIDTH+1)-1:0]	mpy_pipe_d;\n"
		"\t\twire	signed	[(CWIDTH-1):0]	mpy_pipe_vc;\n"
		"\t\twire	signed	[(IWIDTH):0]	mpy_pipe_vd;\n"
		"\t\t//\n"
		"\t\treg	signed	[(CWIDTH+1)-1:0]	mpy_cof_sum;\n"
		"\t\treg	signed	[(IWIDTH+2)-1:0]	mpy_dif_sum;\n"
"\n"
		"\t\tassign	mpy_pipe_vc =  mpy_pipe_c[2*(CWIDTH)-1:CWIDTH];\n"
		"\t\tassign	mpy_pipe_vd =  mpy_pipe_d[2*(IWIDTH+1)-1:IWIDTH+1];\n"
"\n"
		"\t\treg			mpy_pipe_v;\n"
		"\t\treg			ce_phase;\n"
"\n"
		"\t\treg	signed	[(CWIDTH+IWIDTH+1)-1:0]	mpy_pipe_out;\n"
		"\t\treg	signed [IWIDTH+CWIDTH+3-1:0]	longmpy;\n"
"\n"
"\n"
		"\t\tinitial	ce_phase = 1'b0;\n"
		"\t\talways @(posedge i_clk)\n"
		"\t\tif (i_reset)\n"
			"\t\t\tce_phase <= 1'b0;\n"
		"\t\telse if (i_ce)\n"
			"\t\t\tce_phase <= 1'b1;\n"
		"\t\telse\n"
			"\t\t\tce_phase <= 1'b0;\n"
"\n"
		"\t\talways @(*)\n"
			"\t\t\tmpy_pipe_v = (i_ce)||(ce_phase);\n"
"\n"
		"\t\talways @(posedge i_clk)\n"
		"\t\tif (ce_phase)\n"
		"\t\tbegin\n"
			"\t\t\t// Pre-clock\n"
			"\t\t\tmpy_pipe_c[2*CWIDTH-1:0] <=\n"
				"\t\t\t\t\t{ ir_coef_r, ir_coef_i };\n"
			"\t\t\tmpy_pipe_d[2*(IWIDTH+1)-1:0] <=\n"
				"\t\t\t\t\t{ r_dif_r, r_dif_i };\n"
"\n"
			"\t\t\tmpy_cof_sum  <= ir_coef_i + ir_coef_r;\n"
			"\t\t\tmpy_dif_sum <= r_dif_r + r_dif_i;\n"
"\n"
		"\t\tend else if (i_ce)\n"
		"\t\tbegin\n"
			"\t\t\t// First clock\n"
			"\t\t\tmpy_pipe_c[2*(CWIDTH)-1:0] <= {\n"
				"\t\t\t\tmpy_pipe_c[(CWIDTH)-1:0], {(CWIDTH){1'b0}} };\n"
			"\t\t\tmpy_pipe_d[2*(IWIDTH+1)-1:0] <= {\n"
				"\t\t\t\tmpy_pipe_d[(IWIDTH+1)-1:0], {(IWIDTH+1){1'b0}} };\n"
		"\t\tend\n"
"\n");
	fprintf(fp,
		"\t\t%s #(CWIDTH+1,IWIDTH+2) mpy0(i_clk, mpy_pipe_v,\n"
			"\t\t\t\tmpy_cof_sum, mpy_dif_sum, longmpy);\n"
"\n",
		(USE_OLD_MULTIPLY)?"shiftaddmpy":"longbimpy");

	fprintf(fp,
		"\t\t%s #(CWIDTH,IWIDTH+1) mpy1(i_clk, mpy_pipe_v,\n"
			"\t\t\t\tmpy_pipe_vc, mpy_pipe_vd, mpy_pipe_out);\n"
"\n",
		(USE_OLD_MULTIPLY)?"shiftaddmpy":"longbimpy");

	fprintf(fp,
		"\t\treg\tsigned\t[((IWIDTH+1)+(CWIDTH)-1):0]	rp_one;\n"
		"\t\treg\tsigned\t[((IWIDTH+1)+(CWIDTH)-1):0]	rp2_one,\n"
				"\t\t\t\t\t\t\t\trp_two;\n"
		"\t\treg\tsigned\t[((IWIDTH+2)+(CWIDTH+1)-1):0]	rp_three;\n"
"\n"
		"\t\talways @(posedge i_clk)\n"
		"\t\tif (ce_phase) // 1.5 clock\n"
			"\t\t\trp_one <= mpy_pipe_out;\n"
		"\t\talways @(posedge i_clk)\n"
		"\t\tif (i_ce) // two clocks\n"
			"\t\t\trp_two <= mpy_pipe_out;\n"
		"\t\talways @(posedge i_clk)\n"
		"\t\tif (i_ce) // Second clock\n"
			"\t\t\trp_three<= longmpy;\n"
		"\t\talways @(posedge i_clk)\n"
		"\t\tif (i_ce)\n"
			"\t\t\trp2_one<= rp_one;\n"
"\n"
		"\t\tassign	p_one	= { {(2){rp2_one[(IWIDTH+1)+(CWIDTH)-1]}}, rp2_one };\n"
		"\t\tassign	p_two	= { {(2){rp_two[(IWIDTH+1)+(CWIDTH)-1]}}, rp_two };\n"
		"\t\tassign	p_three	= rp_three;\n"
"\n");

	/////////////////////////
	///
	///	Three clock per CE, so CE, no-ce, no-ce*, CE
	///
	fprintf(fp,
"\tend else if (CKPCE <= 3)\n\tbegin : CKPCE_THREE\n");

	fprintf(fp,
	"\t\t// Coefficient multiply inputs\n"
	"\t\treg\t\t[3*(CWIDTH+1)-1:0]\tmpy_pipe_c;\n"
	"\t\t// Data multiply inputs\n"
	"\t\treg\t\t[3*(IWIDTH+2)-1:0]\tmpy_pipe_d;\n"
	"\t\twire\tsigned	[(CWIDTH):0]	mpy_pipe_vc;\n"
	"\t\twire\tsigned	[(IWIDTH+1):0]	mpy_pipe_vd;\n"
	"\n"
	"\t\tassign\tmpy_pipe_vc =  mpy_pipe_c[3*(CWIDTH+1)-1:2*(CWIDTH+1)];\n"
	"\t\tassign\tmpy_pipe_vd =  mpy_pipe_d[3*(IWIDTH+2)-1:2*(IWIDTH+2)];\n"
	"\n"
	"\t\treg\t\t\tmpy_pipe_v;\n"
	"\t\treg\t\t[2:0]\tce_phase;\n"
	"\n"
	"\t\treg\tsigned	[  (CWIDTH+IWIDTH+3)-1:0]	mpy_pipe_out;\n"
"\n");
	fprintf(fp,
	"\t\tinitial\tce_phase = 3'b011;\n"
	"\t\talways @(posedge i_clk)\n"
	"\t\tif (i_reset)\n"
		"\t\t\tce_phase <= 3'b011;\n"
	"\t\telse if (i_ce)\n"
		"\t\t\tce_phase <= 3'b000;\n"
	"\t\telse if (ce_phase != 3'b011)\n"
		"\t\t\tce_phase <= ce_phase + 1'b1;\n"
"\n"
	"\t\talways @(*)\n"
		"\t\t\tmpy_pipe_v = (i_ce)||(ce_phase < 3'b010);\n"
"\n");

	fprintf(fp,
	"\t\talways @(posedge i_clk)\n"
		"\t\t\tif (ce_phase == 3\'b000)\n"
		"\t\t\tbegin\n"
			"\t\t\t\t// Second clock\n"
			"\t\t\t\tmpy_pipe_c[3*(CWIDTH+1)-1:(CWIDTH+1)] <= {\n"
			"\t\t\t\t\tir_coef_r[CWIDTH-1], ir_coef_r,\n"
			"\t\t\t\t\tir_coef_i[CWIDTH-1], ir_coef_i };\n"
			"\t\t\t\tmpy_pipe_c[CWIDTH:0] <= ir_coef_i + ir_coef_r;\n"
			"\t\t\t\tmpy_pipe_d[3*(IWIDTH+2)-1:(IWIDTH+2)] <= {\n"
			"\t\t\t\t\tr_dif_r[IWIDTH], r_dif_r,\n"
			"\t\t\t\t\tr_dif_i[IWIDTH], r_dif_i };\n"
			"\t\t\t\tmpy_pipe_d[(IWIDTH+2)-1:0] <= r_dif_r + r_dif_i;\n"
"\n"
		"\t\t\tend else if (mpy_pipe_v)\n"
		"\t\t\tbegin\n"
			"\t\t\t\tmpy_pipe_c[3*(CWIDTH+1)-1:0] <= {\n"
			"\t\t\t\t\tmpy_pipe_c[2*(CWIDTH+1)-1:0], {(CWIDTH+1){1\'b0}} };\n"
			"\t\t\t\tmpy_pipe_d[3*(IWIDTH+2)-1:0] <= {\n"
			"\t\t\t\t\tmpy_pipe_d[2*(IWIDTH+2)-1:0], {(IWIDTH+2){1\'b0}} };\n"
		"\t\t\tend\n"
"\n");
	fprintf(fp,
		"\t\t%s #(CWIDTH+1,IWIDTH+2) mpy(i_clk, mpy_pipe_v,\n"
			"\t\t\t\tmpy_pipe_vc, mpy_pipe_vd, mpy_pipe_out);\n"
"\n",
		(USE_OLD_MULTIPLY)?"shiftaddmpy":"longbimpy");

	fprintf(fp,
	"\t\treg\tsigned\t[((IWIDTH+1)+(CWIDTH)-1):0]\trp_one, rp_two;\n"
	"\t\treg\tsigned\t[((IWIDTH+1)+(CWIDTH)-1):0]\trp2_one, rp2_two;\n"
	"\t\treg\tsigned\t[((IWIDTH+2)+(CWIDTH+1)-1):0]\trp_three, rp2_three;\n"

"\n");

	fprintf(fp,
	"\t\talways @(posedge i_clk)\n"
	"\t\tif(i_ce)\n"
		"\t\t\trp_one <= mpy_pipe_out[(CWIDTH+IWIDTH+3)-3:0];\n"
	"\t\talways @(posedge i_clk)\n"
	"\t\tif(ce_phase == 3'b000)\n"
		"\t\t\trp_two <= mpy_pipe_out[(CWIDTH+IWIDTH+3)-3:0];\n"
	"\t\talways @(posedge i_clk)\n"
	"\t\tif(ce_phase == 3'b001)\n"
		"\t\t\trp_three <= mpy_pipe_out;\n"
	"\t\talways @(posedge i_clk)\n"
	"\t\tif (i_ce)\n"
	"\t\tbegin\n"
		"\t\t\trp2_one<= rp_one;\n"
		"\t\t\trp2_two<= rp_two;\n"
		"\t\t\trp2_three<= rp_three;\n"
	"\t\tend\n");
	fprintf(fp,

	"\t\tassign\tp_one\t= rp2_one;\n"
	"\t\tassign\tp_two\t= rp2_two;\n"
	"\t\tassign\tp_three\t= rp2_three;\n"
"\n");

	fprintf(fp,
"\tend endgenerate\n");

	fprintf(fp,
	"\t// These values are held in memory and delayed during the\n"
	"\t// multiply.  Here, we recover them.  During the multiply,\n"
	"\t// values were multiplied by 2^(CWIDTH-2)*exp{-j*2*pi*...},\n"
	"\t// therefore, the left_x values need to be right shifted by\n"
	"\t// CWIDTH-2 as well.  The additional bits come from a sign\n"
	"\t// extension.\n"
	"\twire\tsigned\t[(IWIDTH+CWIDTH):0]	fifo_i, fifo_r;\n"
	"\treg\t\t[(2*IWIDTH+1):0]	fifo_read;\n"
	"\tassign\tfifo_r = { {2{fifo_read[2*(IWIDTH+1)-1]}}, fifo_read[(2*(IWIDTH+1)-1):(IWIDTH+1)], {(CWIDTH-2){1\'b0}} };\n"
	"\tassign\tfifo_i = { {2{fifo_read[(IWIDTH+1)-1]}}, fifo_read[((IWIDTH+1)-1):0], {(CWIDTH-2){1\'b0}} };\n"
"\n"
"\n"
	"\treg\tsigned\t[(CWIDTH+IWIDTH+3-1):0]	mpy_r, mpy_i;\n"
"\n");
	fprintf(fp,
	"\t// Let's do some rounding and remove unnecessary bits.\n"
	"\t// We have (IWIDTH+CWIDTH+3) bits here, we need to drop down to\n"
	"\t// OWIDTH, and SHIFT by SHIFT bits in the process.  The trick is\n"
	"\t// that we don\'t need (IWIDTH+CWIDTH+3) bits.  We\'ve accumulated\n"
	"\t// them, but the actual values will never fill all these bits.\n"
	"\t// In particular, we only need:\n"
	"\t//\t IWIDTH bits for the input\n"
	"\t//\t     +1 bit for the add/subtract\n"
	"\t//\t+CWIDTH bits for the coefficient multiply\n"
	"\t//\t     +1 bit for the add/subtract in the complex multiply\n"
	"\t//\t ------\n"
	"\t//\t (IWIDTH+CWIDTH+2) bits at full precision.\n"
	"\t//\n"
	"\t// However, the coefficient multiply multiplied by a maximum value\n"
	"\t// of 2^(CWIDTH-2).  Thus, we only have\n"
	"\t//\t   IWIDTH bits for the input\n"
	"\t//\t       +1 bit for the add/subtract\n"
	"\t//\t+CWIDTH-2 bits for the coefficient multiply\n"
	"\t//\t       +1 (optional) bit for the add/subtract in the cpx mpy.\n"
	"\t//\t -------- ... multiply.  (This last bit may be shifted out.)\n"
	"\t//\t (IWIDTH+CWIDTH) valid output bits.\n"
	"\t// Now, if the user wants to keep any extras of these (via OWIDTH),\n"
	"\t// or if he wishes to arbitrarily shift some of these off (via\n"
	"\t// SHIFT) we accomplish that here.\n"
"\n");
	fprintf(fp,
	"\twire\tsigned\t[(OWIDTH-1):0]\trnd_left_r, rnd_left_i, rnd_right_r, rnd_right_i;\n\n");

	fprintf(fp,
	"\twire\tsigned\t[(CWIDTH+IWIDTH+3-1):0]\tleft_sr, left_si;\n"
	"\tassign	left_sr = { {(2){fifo_r[(IWIDTH+CWIDTH)]}}, fifo_r };\n"
	"\tassign	left_si = { {(2){fifo_i[(IWIDTH+CWIDTH)]}}, fifo_i };\n\n");

	fprintf(fp,
	"\t%s #(CWIDTH+IWIDTH+3,OWIDTH,SHIFT+4) do_rnd_left_r(i_clk, i_ce,\n"
	"\t\t\t\tleft_sr, rnd_left_r);\n\n",
		rnd_string);
	fprintf(fp,
	"\t%s #(CWIDTH+IWIDTH+3,OWIDTH,SHIFT+4) do_rnd_left_i(i_clk, i_ce,\n"
	"\t\t\t\tleft_si, rnd_left_i);\n\n",
		rnd_string);
	fprintf(fp,
	"\t%s #(CWIDTH+IWIDTH+3,OWIDTH,SHIFT+4) do_rnd_right_r(i_clk, i_ce,\n"
	"\t\t\t\tmpy_r, rnd_right_r);\n\n", rnd_string);
	fprintf(fp,
	"\t%s #(CWIDTH+IWIDTH+3,OWIDTH,SHIFT+4) do_rnd_right_i(i_clk, i_ce,\n"
	"\t\t\t\tmpy_i, rnd_right_i);\n\n", rnd_string);
	fprintf(fp,
	"\talways @(posedge i_clk)\n"
		"\t\tif (i_ce)\n"
		"\t\tbegin\n"
			"\t\t\t// First clock, recover all values\n"
			"\t\t\tfifo_read <= fifo_left[fifo_read_addr];\n"
			"\t\t\t// These values are IWIDTH+CWIDTH+3 bits wide\n"
			"\t\t\t// although they only need to be (IWIDTH+1)\n"
			"\t\t\t// + (CWIDTH) bits wide.  (We\'ve got two\n"
			"\t\t\t// extra bits we need to get rid of.)\n"
			"\t\t\tmpy_r <= p_one - p_two;\n"
			"\t\t\tmpy_i <= p_three - p_one - p_two;\n"
		"\t\tend\n"
"\n");

	fprintf(fp,
	"\treg\t[(AUXLEN-1):0]\taux_pipeline;\n"
	"\tinitial\taux_pipeline = 0;\n");
	if (async_reset)
		fprintf(fp, "\talways @(posedge i_clk, negedge i_areset_n)\n\t\tif (!i_areset_n)\n");
	else
		fprintf(fp, "\talways @(posedge i_clk)\n\t\tif (i_reset)\n");
	fprintf(fp,
	"\t\t\taux_pipeline <= 0;\n"
	"\t\telse if (i_ce)\n"
	"\t\t\taux_pipeline <= { aux_pipeline[(AUXLEN-2):0], i_aux };\n"
"\n");
	fprintf(fp,
	"\tinitial o_aux = 1\'b0;\n");
	if (async_reset)
		fprintf(fp, "\talways @(posedge i_clk, negedge i_areset_n)\n\t\tif (!i_areset_n)\n");
	else
		fprintf(fp, "\talways @(posedge i_clk)\n\t\tif (i_reset)\n");
	fprintf(fp,
		"\t\t\to_aux <= 1\'b0;\n"
		"\t\telse if (i_ce)\n"
		"\t\tbegin\n"
			"\t\t\t// Second clock, latch for final clock\n"
			"\t\t\to_aux <= aux_pipeline[AUXLEN-1];\n"
		"\t\tend\n"
"\n");

	fprintf(fp,
	"\t// As a final step, we pack our outputs into two packed two\'s\n"
	"\t// complement numbers per output word, so that each output word\n"
	"\t// has (2*OWIDTH) bits in it, with the top half being the real\n"
	"\t// portion and the bottom half being the imaginary portion.\n"
	"\tassign	o_left = { rnd_left_r, rnd_left_i };\n"
	"\tassign	o_right= { rnd_right_r,rnd_right_i};\n"
"\n");

	if (formal_property_flag) {
		fprintf(fp,
"`ifdef	FORMAL\n"
	"\tlocalparam	F_LGDEPTH = 5;\n"
	"\tlocalparam	F_DEPTH = AUXLEN;\n"
	"\tlocalparam	F_D = F_DEPTH-1;\n"
"\n"
	"\treg	signed	[IWIDTH-1:0]	f_dlyleft_r  [0:F_DEPTH-1];\n"
	"\treg	signed	[IWIDTH-1:0]	f_dlyleft_i  [0:F_DEPTH-1];\n"
	"\treg	signed	[IWIDTH-1:0]	f_dlyright_r [0:F_DEPTH-1];\n"
	"\treg	signed	[IWIDTH-1:0]	f_dlyright_i [0:F_DEPTH-1];\n"
	"\treg	signed	[CWIDTH-1:0]	f_dlycoeff_r [0:F_DEPTH-1];\n"
	"\treg	signed	[CWIDTH-1:0]	f_dlycoeff_i [0:F_DEPTH-1];\n"
	"\treg	signed	[F_DEPTH-1:0]	f_dlyaux;\n"
"\n"
	"\talways @(posedge i_clk)\n"
	"\tif (i_ce)\n"
	"\tbegin\n"
	"\t	f_dlyleft_r[0]   <= i_left[ (2*IWIDTH-1):IWIDTH];\n"
	"\t	f_dlyleft_i[0]   <= i_left[ (  IWIDTH-1):0];\n"
	"\t	f_dlyright_r[0]  <= i_right[(2*IWIDTH-1):IWIDTH];\n"
	"\t	f_dlyright_i[0]  <= i_right[(  IWIDTH-1):0];\n"
	"\t	f_dlycoeff_r[0]  <= i_coef[ (2*CWIDTH-1):CWIDTH];\n"
	"\t	f_dlycoeff_i[0]  <= i_coef[ (  CWIDTH-1):0];\n"
	"\t	f_dlyaux[0]      <= i_aux;\n"
	"\tend\n"
"\n"
	"\tgenvar	k;\n"
	"\tgenerate for(k=1; k<F_DEPTH; k=k+1)\n"
"\n"
	"\t	always @(posedge i_clk)\n"
	"\t	if (i_ce)\n"
	"\t	begin\n"
	"\t		f_dlyleft_r[k]  <= f_dlyleft_r[ k-1];\n"
	"\t		f_dlyleft_i[k]  <= f_dlyleft_i[ k-1];\n"
	"\t		f_dlyright_r[k] <= f_dlyright_r[k-1];\n"
	"\t		f_dlyright_i[k] <= f_dlyright_i[k-1];\n"
	"\t		f_dlycoeff_r[k] <= f_dlycoeff_r[k-1];\n"
	"\t		f_dlycoeff_i[k] <= f_dlycoeff_i[k-1];\n"
	"\t		f_dlyaux[k]     <= f_dlyaux[    k-1];\n"
	"\t	end\n"
"\n"
	"\tendgenerate\n"
"\n"
	"\talways @(posedge i_clk)\n"
	"\tif ((!$past(i_ce))&&(!$past(i_ce,2))&&(!$past(i_ce,3))\n"
	"\t		&&(!$past(i_ce,4)))\n"
	"\t	assume(i_ce);\n"
"\n"
	"\tgenerate if (CKPCE <= 1)\n"
	"\tbegin\n"
"\n"
	"\t	// i_ce is allowed to be anything in this mode\n"
"\n"
	"\tend else if (CKPCE == 2)\n"
	"\tbegin : F_CKPCE_TWO\n"
"\n"
	"\t	always @(posedge i_clk)\n"
	"\t		if ($past(i_ce))\n"
	"\t			assume(!i_ce);\n"
"\n"
	"\tend else if (CKPCE == 3)\n"
	"\tbegin : F_CKPCE_THREE\n"
"\n"
	"\t	always @(posedge i_clk)\n"
	"\t		if (($past(i_ce))||($past(i_ce,2)))\n"
	"\t			assume(!i_ce);\n"
"\n"
	"\tend endgenerate\n"
"\n"
	"\treg	[F_LGDEPTH-1:0]	f_startup_counter;\n"
	"\tinitial	f_startup_counter = 0;\n"
	"\talways @(posedge i_clk)\n"
	"\tif (i_reset)\n"
	"\t	f_startup_counter = 0;\n"
	"\telse if ((i_ce)&&(!(&f_startup_counter)))\n"
	"\t	f_startup_counter <= f_startup_counter + 1;\n"
"\n"
	"\twire	signed	[IWIDTH:0]	f_sumr, f_sumi;\n"
	"\talways @(*)\n"
	"\tbegin\n"
	"\t	f_sumr = f_dlyleft_r[F_D] + f_dlyright_r[F_D];\n"
	"\t	f_sumi = f_dlyleft_i[F_D] + f_dlyright_i[F_D];\n"
	"\tend\n"
"\n"
	"\twire	signed	[IWIDTH+CWIDTH-1:0]	f_sumrx, f_sumix;\n"
	"\tassign	f_sumrx = { {(2){f_sumr[IWIDTH]}}, f_sumr, {(CWIDTH-2){1'b0}} };\n"
	"\tassign	f_sumix = { {(2){f_sumi[IWIDTH]}}, f_sumi, {(CWIDTH-2){1'b0}} };\n"
"\n"
	"\twire	signed	[IWIDTH:0]	f_difr, f_difi;\n"
	"\talways @(*)\n"
	"\tbegin\n"
	"\t	f_difr = f_dlyleft_r[F_D] - f_dlyright_r[F_D];\n"
	"\t	f_difi = f_dlyleft_i[F_D] - f_dlyright_i[F_D];\n"
	"\tend\n"
"\n"
	"\twire	signed	[IWIDTH+CWIDTH-1:0]	f_difrx, f_difix;\n"
	"\tassign	f_difrx = { {(CWIDTH){f_difr[IWIDTH]}}, f_difr };\n"
	"\tassign	f_difix = { {(CWIDTH){f_difi[IWIDTH]}}, f_difi };\n"
"\n"
"\n"
	"\talways @(posedge i_clk)\n"
	"\tif (f_startup_counter > F_D)\n"
	"\tbegin\n"
	"\t	assert(left_sr == f_sumrx);\n"
	"\t	assert(left_si == f_sumix);\n"
	"\t	assert(aux_pipeline[AUXLEN-1] == f_dlyaux[F_D]);\n"
"\n"
	"\t	if ((f_difr == 0)&&(f_difi == 0))\n"
	"\t	begin\n"
	"\t		assert(mpy_r == 0);\n"
	"\t		assert(mpy_i == 0);\n"
	"\t	end else if ((f_dlycoeff_r[F_D] == 0)\n"
	"\t			&&(f_dlycoeff_i[F_D] == 0))\n"
	"\t	begin\n"
	"\t		assert(mpy_r == 0);\n"
	"\t		assert(mpy_i == 0);\n"
	"\t	end\n"
"\n"
	"\t	if ((f_dlycoeff_r[F_D] == 1)&&(f_dlycoeff_i[F_D] == 0))\n"
	"\t	begin\n"
	"\t		assert(mpy_r == f_difrx);\n"
	"\t		assert(mpy_i == f_difix);\n"
	"\t	end\n"
"\n"
	"\t	if ((f_dlycoeff_r[F_D] == 0)&&(f_dlycoeff_i[F_D] == 1))\n"
	"\t	begin\n"
	"\t		assert(mpy_r == -f_difix);\n"
	"\t		assert(mpy_i ==  f_difrx);\n"
	"\t	end\n"
"\n"
	"\t	if ((f_difr == 1)&&(f_difi == 0))\n"
	"\t	begin\n"
	"\t		assert(mpy_r == f_dlycoeff_r[F_D]);\n"
	"\t		assert(mpy_i == f_dlycoeff_i[F_D]);\n"
	"\t	end\n"
"\n"
	"\t	if ((f_difr == 0)&&(f_difi == 1))\n"
	"\t	begin\n"
	"\t		assert(mpy_r == -f_dlycoeff_i[F_D]);\n"
	"\t		assert(mpy_i ==  f_dlycoeff_r[F_D]);\n"
	"\t	end\n"
	"\tend\n"
"\n"
"`endif\n");
	}

	fprintf(fp,
"endmodule\n");
	fclose(fp);
}

void	build_hwbfly(const char *fname, int xtracbits, ROUND_T rounding,
		int ckpce, const bool async_reset) {
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
"// Filename:\thwbfly.v\n"
"//\n"
"// Project:\t%s\n"
"//\n"
"// Purpose:\tThis routine is identical to the butterfly.v routine found\n"
"//		in 'butterfly.v', save only that it uses the verilog\n"
"//	operator '*' in hopes that the synthesizer would be able to optimize\n"
"//	it with hardware resources.\n"
"//\n"
"//	It is understood that a hardware multiply can complete its operation in\n"
"//	a single clock.\n"
"//\n"
"//\n%s"
"//\n", prjname, creator);
	fprintf(fp, "%s", cpyleft);
	fprintf(fp, "//\n//\n`default_nettype\tnone\n//\n");
	fprintf(fp,
"module	hwbfly(i_clk, %s, i_ce, i_coef, i_left, i_right, i_aux,\n"
		"\t\to_left, o_right, o_aux);\n"
	"\t// Public changeable parameters ...\n"
	"\tparameter IWIDTH=16,CWIDTH=IWIDTH+%d,OWIDTH=IWIDTH+1;\n"
	"\t// Parameters specific to the core that should not be changed.\n"
	"\tparameter\tSHIFT=0;\n"
	"\tparameter\t[1:0]\tCKPCE=%d;\n"
	"\tinput\t\ti_clk, %s, i_ce;\n"
	"\tinput\t\t[(2*CWIDTH-1):0]\ti_coef;\n"
	"\tinput\t\t[(2*IWIDTH-1):0]\ti_left, i_right;\n"
	"\tinput\t\ti_aux;\n"
	"\toutput\twire\t[(2*OWIDTH-1):0]\to_left, o_right;\n"
	"\toutput\treg\to_aux;\n"
"\n", resetw.c_str(), xtracbits, ckpce, resetw.c_str());
	fprintf(fp,
	"\treg\t[(2*IWIDTH-1):0]	r_left, r_right;\n"
	"\treg\t			r_aux, r_aux_2;\n"
	"\treg\t[(2*CWIDTH-1):0]	r_coef;\n"
	"\twire	signed	[(IWIDTH-1):0]	r_left_r, r_left_i, r_right_r, r_right_i;\n"
	"\tassign\tr_left_r  = r_left[ (2*IWIDTH-1):(IWIDTH)];\n"
	"\tassign\tr_left_i  = r_left[ (IWIDTH-1):0];\n"
	"\tassign\tr_right_r = r_right[(2*IWIDTH-1):(IWIDTH)];\n"
	"\tassign\tr_right_i = r_right[(IWIDTH-1):0];\n"
	"\treg	signed	[(CWIDTH-1):0]	ir_coef_r, ir_coef_i;\n"
"\n"
	"\treg	signed	[(IWIDTH):0]	r_sum_r, r_sum_i, r_dif_r, r_dif_i;\n"
"\n"
	"\treg	[(2*IWIDTH+2):0]	leftv, leftvv;\n"
"\n"
	"\t// Set up the input to the multiply\n"
	"\tinitial r_aux   = 1\'b0;\n"
	"\tinitial r_aux_2 = 1\'b0;\n");
	if (async_reset)
		fprintf(fp, "\talways @(posedge i_clk, negedge i_areset_n)\n\t\tif (!i_areset_n)\n");
	else
		fprintf(fp, "\talways @(posedge i_clk)\n\t\tif (i_reset)\n");
	fprintf(fp,
		"\t\tbegin\n"
			"\t\t\tr_aux <= 1\'b0;\n"
			"\t\t\tr_aux_2 <= 1\'b0;\n"
		"\t\tend else if (i_ce)\n"
		"\t\tbegin\n"
			"\t\t\t// One clock just latches the inputs\n"
			"\t\t\tr_aux <= i_aux;\n"
			"\t\t\t// Next clock adds/subtracts\n"
			"\t\t\t// Other inputs are simply delayed on second clock\n"
			"\t\t\tr_aux_2 <= r_aux;\n"
		"\t\tend\n"
	"\talways @(posedge i_clk)\n"
		"\t\tif (i_ce)\n"
		"\t\tbegin\n"
			"\t\t\t// One clock just latches the inputs\n"
			"\t\t\tr_left <= i_left;	// No change in # of bits\n"
			"\t\t\tr_right <= i_right;\n"
			"\t\t\tr_coef  <= i_coef;\n"
			"\t\t\t// Next clock adds/subtracts\n"
			"\t\t\tr_sum_r <= r_left_r + r_right_r; // Now IWIDTH+1 bits\n"
			"\t\t\tr_sum_i <= r_left_i + r_right_i;\n"
			"\t\t\tr_dif_r <= r_left_r - r_right_r;\n"
			"\t\t\tr_dif_i <= r_left_i - r_right_i;\n"
			"\t\t\t// Other inputs are simply delayed on second clock\n"
			"\t\t\tir_coef_r <= r_coef[(2*CWIDTH-1):CWIDTH];\n"
			"\t\t\tir_coef_i <= r_coef[(CWIDTH-1):0];\n"
		"\t\tend\n"
	"\n\n");
	fprintf(fp,
"\t// See comments in the butterfly.v source file for a discussion of\n"
"\t// these operations and the appropriate bit widths.\n\n");
	fprintf(fp,
	"\twire\tsigned	[((IWIDTH+1)+(CWIDTH)-1):0]	p_one, p_two;\n"
	"\twire\tsigned	[((IWIDTH+2)+(CWIDTH+1)-1):0]	p_three;\n"
"\n"
	"\tinitial leftv    = 0;\n"
	"\tinitial leftvv   = 0;\n");
	if (async_reset)
		fprintf(fp, "\talways @(posedge i_clk, negedge i_areset_n)\n\t\tif (!i_areset_n)\n");
	else
		fprintf(fp, "\talways @(posedge i_clk)\n\t\tif (i_reset)\n");
	fprintf(fp,
		"\t\tbegin\n"
			"\t\t\tleftv <= 0;\n"
			"\t\t\tleftvv <= 0;\n"
		"\t\tend else if (i_ce)\n"
		"\t\tbegin\n"
			"\t\t\t// Second clock, pipeline = 1\n"
			"\t\t\tleftv <= { r_aux_2, r_sum_r, r_sum_i };\n"
"\n"
			"\t\t\t// Third clock, pipeline = 3\n"
			"\t\t\t//   As desired, each of these lines infers a DSP48\n"
			"\t\t\tleftvv <= leftv;\n"
		"\t\tend\n"
"\n");

	// Nominally, we should handle code for 1, 2, or 3 clocks per CE, with
	// one clock per CE meaning CE could be constant.  The code below
	// instead handles 1 or 3 clocks per CE, leaving the two clocks per
	// CE optimization(s) unfulfilled.

//	fprintf(fp,
//"\tend else if (CKPCI == 2'b01)\n\tbegin\n");

	///////////////////////////////////////////
	///
	///	One clock per CE, so CE, CE, CE, CE, CE is possible
	///
	fprintf(fp,
"\tgenerate if (CKPCE <= 2'b01)\n\tbegin : CKPCE_ONE\n");

	fprintf(fp,
	"\t\t// Coefficient multiply inputs\n"
	"\t\treg\tsigned	[(CWIDTH-1):0]	p1c_in, p2c_in;\n"
	"\t\t// Data multiply inputs\n"
	"\t\treg\tsigned	[(IWIDTH):0]	p1d_in, p2d_in;\n"
	"\t\t// Product 3, coefficient input\n"
	"\t\treg\tsigned	[(CWIDTH):0]	p3c_in;\n"
	"\t\t// Product 3, data input\n"
	"\t\treg\tsigned	[(IWIDTH+1):0]	p3d_in;\n"
"\n");
	fprintf(fp,
	"\t\treg\tsigned	[((IWIDTH+1)+(CWIDTH)-1):0]	rp_one, rp_two;\n"
	"\t\treg\tsigned	[((IWIDTH+2)+(CWIDTH+1)-1):0]	rp_three;\n"
"\n");

	fprintf(fp,
	"\t\talways @(posedge i_clk)\n"
	"\t\tif (i_ce)\n"
	"\t\tbegin\n"
		"\t\t\t// Second clock, pipeline = 1\n"
		"\t\t\tp1c_in <= ir_coef_r;\n"
		"\t\t\tp2c_in <= ir_coef_i;\n"
		"\t\t\tp1d_in <= r_dif_r;\n"
		"\t\t\tp2d_in <= r_dif_i;\n"
		"\t\t\tp3c_in <= ir_coef_i + ir_coef_r;\n"
		"\t\t\tp3d_in <= r_dif_r + r_dif_i;\n"
	"\t\tend\n\n");

	if (formal_property_flag)
		fprintf(fp,
"`ifndef	FORMAL\n");

	fprintf(fp,
	"\t\talways @(posedge i_clk)\n"
	"\t\tif (i_ce)\n"
	"\t\tbegin\n"
		"\t\t\t// Third clock, pipeline = 3\n"
		"\t\t\t//   As desired, each of these lines infers a DSP48\n"
		"\t\t\trp_one   <= p1c_in * p1d_in;\n"
		"\t\t\trp_two   <= p2c_in * p2d_in;\n"
		"\t\t\trp_three <= p3c_in * p3d_in;\n"
	"\t\tend\n");

	if (formal_property_flag)
		fprintf(fp,
"`else\n"
		"\t\twire	signed	[((IWIDTH+1)+(CWIDTH)-1):0]	pre_rp_one, pre_rp_two;\n"
		"\t\twire	signed	[((IWIDTH+2)+(CWIDTH+1)-1):0]	pre_rp_three;\n"
"\n"
		"\t\tabs_mpy #(CWIDTH,IWIDTH+1,1'b1)\n"
		"\t\t	onei(p1c_in, p1d_in, pre_rp_one);\n"
		"\t\tabs_mpy #(CWIDTH,IWIDTH+1,1'b1)\n"
		"\t\t	twoi(p2c_in, p2d_in, pre_rp_two);\n"
		"\t\tabs_mpy #(CWIDTH+1,IWIDTH+2,1'b1)\n"
		"\t\t	threei(p3c_in, p3d_in, pre_rp_three);\n"
"\n"
		"\t\talways @(posedge i_clk)\n"
		"\t\tif (i_ce)\n"
		"\t\tbegin\n"
		"\t\t	rp_one   = pre_rp_one;\n"
		"\t\t	rp_two   = pre_rp_two;\n"
		"\t\t	rp_three = pre_rp_three;\n"
		"\t\tend\n"
"`endif // FORMAL\n");

	fprintf(fp,"\n"
	"\t\tassign\tp_one   = rp_one;\n"
	"\t\tassign\tp_two   = rp_two;\n"
	"\t\tassign\tp_three = rp_three;\n"
"\n");

	///////////////////////////////////////////
	///
	///	Two clocks per CE, so CE, no-ce, CE, no-ce, etc
	///
	fprintf(fp,
	"\tend else if (CKPCE <= 2'b10)\n"
	"\tbegin : CKPCE_TWO\n"
		"\t\t// Coefficient multiply inputs\n"
		"\t\treg		[2*(CWIDTH)-1:0]	mpy_pipe_c;\n"
		"\t\t// Data multiply inputs\n"
		"\t\treg		[2*(IWIDTH+1)-1:0]	mpy_pipe_d;\n"
		"\t\twire	signed	[(CWIDTH-1):0]	mpy_pipe_vc;\n"
		"\t\twire	signed	[(IWIDTH):0]	mpy_pipe_vd;\n"
		"\t\t//\n"
		"\t\treg	signed	[(CWIDTH+1)-1:0]	mpy_cof_sum;\n"
		"\t\treg	signed	[(IWIDTH+2)-1:0]	mpy_dif_sum;\n"
"\n"
		"\t\tassign	mpy_pipe_vc =  mpy_pipe_c[2*(CWIDTH)-1:CWIDTH];\n"
		"\t\tassign	mpy_pipe_vd =  mpy_pipe_d[2*(IWIDTH+1)-1:IWIDTH+1];\n"
"\n"
		"\t\treg			mpy_pipe_v;\n"
		"\t\treg			ce_phase;\n"
"\n"
		"\t\treg	signed	[(CWIDTH+IWIDTH+1)-1:0]	mpy_pipe_out;\n"
		"\t\treg	signed [IWIDTH+CWIDTH+3-1:0]	longmpy;\n"
"\n"
"\n"
		"\t\tinitial	ce_phase = 1'b1;\n"
		"\t\talways @(posedge i_clk)\n"
		"\t\tif (i_reset)\n"
			"\t\t\tce_phase <= 1'b1;\n"
		"\t\telse if (i_ce)\n"
			"\t\t\tce_phase <= 1'b0;\n"
		"\t\telse\n"
			"\t\t\tce_phase <= 1'b1;\n"
"\n"
		"\t\talways @(*)\n"
			"\t\t\tmpy_pipe_v = (i_ce)||(!ce_phase);\n"
"\n"
		"\t\talways @(posedge i_clk)\n"
		"\t\tif (!ce_phase)\n"
		"\t\tbegin\n"
			"\t\t\t// Pre-clock\n"
			"\t\t\tmpy_pipe_c[2*CWIDTH-1:0] <=\n"
				"\t\t\t\t\t{ ir_coef_r, ir_coef_i };\n"
			"\t\t\tmpy_pipe_d[2*(IWIDTH+1)-1:0] <=\n"
				"\t\t\t\t\t{ r_dif_r, r_dif_i };\n"
"\n"
			"\t\t\tmpy_cof_sum  <= ir_coef_i + ir_coef_r;\n"
			"\t\t\tmpy_dif_sum <= r_dif_r + r_dif_i;\n"
"\n"
		"\t\tend else if (i_ce)\n"
		"\t\tbegin\n"
			"\t\t\t// First clock\n"
			"\t\t\tmpy_pipe_c[2*(CWIDTH)-1:0] <= {\n"
				"\t\t\t\tmpy_pipe_c[(CWIDTH)-1:0], {(CWIDTH){1'b0}} };\n"
			"\t\t\tmpy_pipe_d[2*(IWIDTH+1)-1:0] <= {\n"
				"\t\t\t\tmpy_pipe_d[(IWIDTH+1)-1:0], {(IWIDTH+1){1'b0}} };\n"
		"\t\tend\n\n");

	if (formal_property_flag)
		fprintf(fp, "`ifndef	FORMAL\n");

	fprintf(fp,
		"\t\talways @(posedge i_clk)\n"
		"\t\tif (i_ce) // First clock\n"
			"\t\t\tlongmpy <= mpy_cof_sum * mpy_dif_sum;\n"
"\n"
		"\t\talways @(posedge i_clk)\n"
		"\t\tif (mpy_pipe_v)\n"
			"\t\t\tmpy_pipe_out <= mpy_pipe_vc * mpy_pipe_vd;\n");

	if (formal_property_flag)
		fprintf(fp, "`else\n"
		"\t\twire	signed [IWIDTH+CWIDTH+3-1:0]	pre_longmpy;\n"
		"\t\twire	signed	[(CWIDTH+IWIDTH+1)-1:0]	pre_mpy_pipe_out;\n"
"\n"
		"\t\tabs_mpy	#(CWIDTH+1,IWIDTH+2,1)\n"
		"\t\t	longmpyi(mpy_cof_sum, mpy_dif_sum, pre_longmpy);\n"
"\n"
		"\t\talways @(posedge i_clk)\n"
		"\t\tif (i_ce)\n"
		"\t\t	longmpy <= pre_longmpy;\n"
"\n"
"\n"
		"\t\tabs_mpy #(CWIDTH,IWIDTH+1,1)\n"
		"\t\t	mpy_pipe_outi(mpy_pipe_vc, mpy_pipe_vd, pre_mpy_pipe_out);\n"
"\n"
		"\t\talways @(posedge i_clk)\n"
		"\t\tif (mpy_pipe_v)\n"
		"\t\t	mpy_pipe_out <= pre_mpy_pipe_out;\n"
"`endif\n");

	fprintf(fp,"\n"
		"\t\treg\tsigned\t[((IWIDTH+1)+(CWIDTH)-1):0]	rp_one;\n"
		"\t\treg\tsigned\t[((IWIDTH+1)+(CWIDTH)-1):0]	rp2_one,\n"
				"\t\t\t\t\t\t\t\trp_two;\n"
		"\t\treg\tsigned\t[((IWIDTH+2)+(CWIDTH+1)-1):0]	rp_three;\n"
"\n"
		"\t\talways @(posedge i_clk)\n"
		"\t\tif (!ce_phase) // 1.5 clock\n"
			"\t\t\trp_one <= mpy_pipe_out;\n"
		"\t\talways @(posedge i_clk)\n"
		"\t\tif (i_ce) // two clocks\n"
			"\t\t\trp_two <= mpy_pipe_out;\n"
		"\t\talways @(posedge i_clk)\n"
		"\t\tif (i_ce) // Second clock\n"
			"\t\t\trp_three<= longmpy;\n"
		"\t\talways @(posedge i_clk)\n"
		"\t\tif (i_ce)\n"
			"\t\t\trp2_one<= rp_one;\n"
"\n"
		"\t\tassign	p_one	= rp2_one;\n"
		"\t\tassign	p_two	= rp_two;\n"
		"\t\tassign	p_three	= rp_three;\n"
"\n");

	/////////////////////////
	///
	///	Three clock per CE, so CE, no-ce, no-ce*, CE
	///
	fprintf(fp,
"\tend else if (CKPCE <= 2'b11)\n\tbegin : CKPCE_THREE\n");

	fprintf(fp,
	"\t\t// Coefficient multiply inputs\n"
	"\t\treg\t\t[3*(CWIDTH+1)-1:0]\tmpy_pipe_c;\n"
	"\t\t// Data multiply inputs\n"
	"\t\treg\t\t[3*(IWIDTH+2)-1:0]\tmpy_pipe_d;\n"
	"\t\twire\tsigned	[(CWIDTH):0]	mpy_pipe_vc;\n"
	"\t\twire\tsigned	[(IWIDTH+1):0]	mpy_pipe_vd;\n"
	"\n"
	"\t\tassign\tmpy_pipe_vc =  mpy_pipe_c[3*(CWIDTH+1)-1:2*(CWIDTH+1)];\n"
	"\t\tassign\tmpy_pipe_vd =  mpy_pipe_d[3*(IWIDTH+2)-1:2*(IWIDTH+2)];\n"
	"\n"
	"\t\treg\t\t\tmpy_pipe_v;\n"
	"\t\treg\t\t[2:0]\tce_phase;\n"
	"\n"
	"\t\treg\tsigned	[  (CWIDTH+IWIDTH+3)-1:0]	mpy_pipe_out;\n"
"\n");
	fprintf(fp,
	"\t\tinitial\tce_phase = 3'b011;\n"
	"\t\talways @(posedge i_clk)\n"
	"\t\tif (i_reset)\n"
		"\t\t\tce_phase <= 3'b011;\n"
	"\t\telse if (i_ce)\n"
		"\t\t\tce_phase <= 3'b000;\n"
	"\t\telse if (ce_phase != 3'b011)\n"
		"\t\t\tce_phase <= ce_phase + 1'b1;\n"
"\n"
	"\t\talways @(*)\n"
		"\t\t\tmpy_pipe_v = (i_ce)||(ce_phase < 3'b010);\n"
"\n");

	fprintf(fp,
	"\t\talways @(posedge i_clk)\n"
		"\t\t\tif (ce_phase == 3\'b000)\n"
		"\t\t\tbegin\n"
			"\t\t\t\t// Second clock\n"
			"\t\t\t\tmpy_pipe_c[3*(CWIDTH+1)-1:(CWIDTH+1)] <= {\n"
			"\t\t\t\t\tir_coef_r[CWIDTH-1], ir_coef_r,\n"
			"\t\t\t\t\tir_coef_i[CWIDTH-1], ir_coef_i };\n"
			"\t\t\t\tmpy_pipe_c[CWIDTH:0] <= ir_coef_i + ir_coef_r;\n"
			"\t\t\t\tmpy_pipe_d[3*(IWIDTH+2)-1:(IWIDTH+2)] <= {\n"
			"\t\t\t\t\tr_dif_r[IWIDTH], r_dif_r,\n"
			"\t\t\t\t\tr_dif_i[IWIDTH], r_dif_i };\n"
			"\t\t\t\tmpy_pipe_d[(IWIDTH+2)-1:0] <= r_dif_r + r_dif_i;\n"
"\n"
		"\t\t\tend else if (mpy_pipe_v)\n"
		"\t\t\tbegin\n"
			"\t\t\t\tmpy_pipe_c[3*(CWIDTH+1)-1:0] <= {\n"
			"\t\t\t\t\tmpy_pipe_c[2*(CWIDTH+1)-1:0], {(CWIDTH+1){1\'b0}} };\n"
			"\t\t\t\tmpy_pipe_d[3*(IWIDTH+2)-1:0] <= {\n"
			"\t\t\t\t\tmpy_pipe_d[2*(IWIDTH+2)-1:0], {(IWIDTH+2){1\'b0}} };\n"
		"\t\t\tend\n\n");

	if (formal_property_flag)
		fprintf(fp, "`ifndef\tFORMAL\n");

	fprintf(fp,
	"\t\talways @(posedge i_clk)\n"
	"\t\t\tif (mpy_pipe_v)\n"
			"\t\t\t\tmpy_pipe_out <= mpy_pipe_vc * mpy_pipe_vd;\n"
"\n");

	if (formal_property_flag)
		fprintf(fp,
"`else\t// FORMAL\n"
		"\t\twire	signed	[  (CWIDTH+IWIDTH+3)-1:0] pre_mpy_pipe_out;\n"
"\n"
		"\t\tabs_mpy #(CWIDTH+1,IWIDTH+2,1)\n"
		"\t\t	mpy_pipe_outi(mpy_pipe_vc, mpy_pipe_vd, pre_mpy_pipe_out);\n"
		"\t\talways @(posedge i_clk)\n"
		"\t\t	if (mpy_pipe_v)\n"
		"\t\t		mpy_pipe_out <= pre_mpy_pipe_out;\n"
"`endif\t// FORMAL\n\n");


	fprintf(fp,
	"\t\treg\tsigned\t[((IWIDTH+1)+(CWIDTH)-1):0]\trp_one, rp_two;\n"
	"\t\treg\tsigned\t[((IWIDTH+1)+(CWIDTH)-1):0]\trp2_one, rp2_two;\n"
	"\t\treg\tsigned\t[((IWIDTH+2)+(CWIDTH+1)-1):0]\trp_three, rp2_three;\n"

"\n");

	fprintf(fp,
	"\t\talways @(posedge i_clk)\n"
	"\t\tif(i_ce)\n"
		"\t\t\trp_one <= mpy_pipe_out[(CWIDTH+IWIDTH+3)-3:0];\n"
	"\t\talways @(posedge i_clk)\n"
	"\t\tif(ce_phase == 3'b000)\n"
		"\t\t\trp_two <= mpy_pipe_out[(CWIDTH+IWIDTH+3)-3:0];\n"
	"\t\talways @(posedge i_clk)\n"
	"\t\tif(ce_phase == 3'b001)\n"
		"\t\t\trp_three <= mpy_pipe_out;\n"
	"\t\talways @(posedge i_clk)\n"
	"\t\tif (i_ce)\n"
	"\t\tbegin\n"
		"\t\t\trp2_one<= rp_one;\n"
		"\t\t\trp2_two<= rp_two;\n"
		"\t\t\trp2_three<= rp_three;\n"
	"\t\tend\n");
	fprintf(fp,

	"\t\tassign\tp_one\t= rp2_one;\n"
	"\t\tassign\tp_two\t= rp2_two;\n"
	"\t\tassign\tp_three\t= rp2_three;\n"
"\n");

	fprintf(fp,
"\tend endgenerate\n");

	fprintf(fp,
	"\twire\tsigned	[((IWIDTH+2)+(CWIDTH+1)-1):0]	w_one, w_two;\n"
	"\tassign\tw_one = { {(2){p_one[((IWIDTH+1)+(CWIDTH)-1)]}}, p_one };\n"
	"\tassign\tw_two = { {(2){p_two[((IWIDTH+1)+(CWIDTH)-1)]}}, p_two };\n"
"\n");

	fprintf(fp,
	"\t// These values are held in memory and delayed during the\n"
	"\t// multiply.  Here, we recover them.  During the multiply,\n"
	"\t// values were multiplied by 2^(CWIDTH-2)*exp{-j*2*pi*...},\n"
	"\t// therefore, the left_x values need to be right shifted by\n"
	"\t// CWIDTH-2 as well.  The additional bits come from a sign\n"
	"\t// extension.\n"
	"\twire\taux_s;\n"
	"\twire\tsigned\t[(IWIDTH+CWIDTH):0]	left_si, left_sr;\n"
	"\treg\t\t[(2*IWIDTH+2):0]	left_saved;\n"
	"\tassign\tleft_sr = { {2{left_saved[2*(IWIDTH+1)-1]}}, left_saved[(2*(IWIDTH+1)-1):(IWIDTH+1)], {(CWIDTH-2){1\'b0}} };\n"
	"\tassign\tleft_si = { {2{left_saved[(IWIDTH+1)-1]}}, left_saved[((IWIDTH+1)-1):0], {(CWIDTH-2){1\'b0}} };\n"
	"\tassign\taux_s = left_saved[2*IWIDTH+2];\n"
"\n"
"\n");

	fprintf(fp,
	"\tinitial left_saved = 0;\n"
	"\tinitial o_aux      = 1\'b0;\n");
	if (async_reset)
		fprintf(fp, "\talways @(posedge i_clk, negedge i_areset_n)\n\t\tif (!i_areset_n)\n");
	else
		fprintf(fp, "\talways @(posedge i_clk)\n\t\tif (i_reset)\n");
	fprintf(fp,
	"\t\tbegin\n"
		"\t\t\tleft_saved <= 0;\n"
		"\t\t\to_aux <= 1\'b0;\n"
	"\t\tend else if (i_ce)\n"
	"\t\tbegin\n"
		"\t\t\t// First clock, recover all values\n"
		"\t\t\tleft_saved <= leftvv;\n"
"\n"
		"\t\t\t// Second clock, round and latch for final clock\n"
		"\t\t\to_aux <= aux_s;\n"
	"\t\tend\n"
	"\talways @(posedge i_clk)\n"
	"\t\tif (i_ce)\n"
	"\t\tbegin\n"
		"\t\t\t// These values are IWIDTH+CWIDTH+3 bits wide\n"
		"\t\t\t// although they only need to be (IWIDTH+1)\n"
		"\t\t\t// + (CWIDTH) bits wide.  (We've got two\n"
		"\t\t\t// extra bits we need to get rid of.)\n"
		"\n"
		"\t\t\t// These two lines also infer DSP48\'s.\n"
		"\t\t\t// To keep from using extra DSP48 resources,\n"
		"\t\t\t// they are prevented from using DSP48\'s\n"
		"\t\t\t// by the (* use_dsp48 ... *) comment above.\n"
		"\t\t\tmpy_r <= w_one - w_two;\n"
		"\t\t\tmpy_i <= p_three - w_one - w_two;\n"
	"\t\tend\n"
	"\n");

	fprintf(fp,
	"\t// Round the results\n"
	"\t(* use_dsp48=\"no\" *)\n"
	"\treg	signed	[(CWIDTH+IWIDTH+3-1):0]	mpy_r, mpy_i;\n");
	fprintf(fp,
	"\twire\tsigned\t[(OWIDTH-1):0]\trnd_left_r, rnd_left_i, rnd_right_r, rnd_right_i;\n\n");
	fprintf(fp,
	"\t%s #(CWIDTH+IWIDTH+1,OWIDTH,SHIFT+2) do_rnd_left_r(i_clk, i_ce,\n"
	"\t\t\t\tleft_sr, rnd_left_r);\n\n",
		rnd_string);
	fprintf(fp,
	"\t%s #(CWIDTH+IWIDTH+1,OWIDTH,SHIFT+2) do_rnd_left_i(i_clk, i_ce,\n"
	"\t\t\t\tleft_si, rnd_left_i);\n\n",
		rnd_string);
	fprintf(fp,
	"\t%s #(CWIDTH+IWIDTH+3,OWIDTH,SHIFT+4) do_rnd_right_r(i_clk, i_ce,\n"
	"\t\t\t\tmpy_r, rnd_right_r);\n\n", rnd_string);
	fprintf(fp,
	"\t%s #(CWIDTH+IWIDTH+3,OWIDTH,SHIFT+4) do_rnd_right_i(i_clk, i_ce,\n"
	"\t\t\t\tmpy_i, rnd_right_i);\n\n", rnd_string);


	fprintf(fp,
	"\t// As a final step, we pack our outputs into two packed two's\n"
	"\t// complement numbers per output word, so that each output word\n"
	"\t// has (2*OWIDTH) bits in it, with the top half being the real\n"
	"\t// portion and the bottom half being the imaginary portion.\n"
	"\tassign\to_left = { rnd_left_r, rnd_left_i };\n"
	"\tassign\to_right= { rnd_right_r,rnd_right_i};\n"
"\n");

	if (formal_property_flag) {
		fprintf(fp,
"`ifdef	FORMAL\n"
	"\tparameter	F_DEPTH = 6;\n"
	"\tlocalparam	F_D = F_DEPTH-2;\n"
"\n"
	"\treg	signed	[IWIDTH-1:0]	f_dlyleft_r  [0:F_DEPTH-1];\n"
	"\treg	signed	[IWIDTH-1:0]	f_dlyleft_i  [0:F_DEPTH-1];\n"
	"\treg	signed	[IWIDTH-1:0]	f_dlyright_r [0:F_DEPTH-1];\n"
	"\treg	signed	[IWIDTH-1:0]	f_dlyright_i [0:F_DEPTH-1];\n"
	"\treg	signed	[CWIDTH-1:0]	f_dlycoeff_r [0:F_DEPTH-1];\n"
	"\treg	signed	[CWIDTH-1:0]	f_dlycoeff_i [0:F_DEPTH-1];\n"
	"\treg	signed	[F_DEPTH-1:0]	f_dlyaux;\n"
"\n"
	"\talways @(posedge i_clk)\n"
	"\tif (i_ce)\n"
	"\tbegin\n"
		"\t\tf_dlyleft_r[0]   <= i_left[ (2*IWIDTH-1):IWIDTH];\n"
		"\t\tf_dlyleft_i[0]   <= i_left[ (  IWIDTH-1):0];\n"
		"\t\tf_dlyright_r[0]  <= i_right[(2*IWIDTH-1):IWIDTH];\n"
		"\t\tf_dlyright_i[0]  <= i_right[(  IWIDTH-1):0];\n"
		"\t\tf_dlycoeff_r[0]  <= i_coef[ (2*CWIDTH-1):CWIDTH];\n"
		"\t\tf_dlycoeff_i[0]  <= i_coef[ (  CWIDTH-1):0];\n"
		"\t\tf_dlyaux[0]      <= i_aux;\n"
	"\tend\n"
"\n"
	"\tgenvar	k;\n"
	"\tgenerate for(k=1; k<F_DEPTH; k=k+1)\n"
"\n"
		"\t\talways @(posedge i_clk)\n"
		"\t\tif (i_ce)\n"
		"\t\tbegin\n"
			"\t\t\tf_dlyleft_r[k]  <= f_dlyleft_r[ k-1];\n"
			"\t\t\tf_dlyleft_i[k]  <= f_dlyleft_i[ k-1];\n"
			"\t\t\tf_dlyright_r[k] <= f_dlyright_r[k-1];\n"
			"\t\t\tf_dlyright_i[k] <= f_dlyright_i[k-1];\n"
			"\t\t\tf_dlycoeff_r[k] <= f_dlycoeff_r[k-1];\n"
			"\t\t\tf_dlycoeff_i[k] <= f_dlycoeff_i[k-1];\n"
			"\t\t\tf_dlyaux[k]     <= f_dlyaux[    k-1];\n"
		"\t\tend\n"
"\n"
	"\tendgenerate\n"
"\n"
	"\talways @(posedge i_clk)\n"
	"\tif ((!$past(i_ce))&&(!$past(i_ce,2))&&(!$past(i_ce,3))\n"
			"\t\t\t&&(!$past(i_ce,4)))\n"
		"\t\tassume(i_ce);\n"
"\n"
	"\tgenerate if (CKPCE <= 1)\n"
	"\tbegin\n"
"\n"
	"\t\t// i_ce is allowed to be anything in this mode\n"
"\n"
	"\tend else if (CKPCE == 2)\n"
	"\tbegin : F_CKPCE_TWO\n"
"\n"
		"\t\talways @(posedge i_clk)\n"
		"\t\t	if ($past(i_ce))\n"
		"\t\t		assume(!i_ce);\n"
	"\n"
	"\tend else if (CKPCE == 3)\n"
	"\tbegin : F_CKPCE_THREE\n"
"\n"
		"\t\talways @(posedge i_clk)\n"
		"\t\t	if (($past(i_ce))||($past(i_ce,2)))\n"
		"\t\t		assume(!i_ce);\n"
"\n"
	"\tend endgenerate\n"
"\n"
	"\treg	[4:0]	f_startup_counter;\n"
	"\tinitial	f_startup_counter = 0;\n"
	"\talways @(posedge i_clk)\n"
	"\tif (i_reset)\n"
		"\t\tf_startup_counter = 0;\n"
	"\telse if ((i_ce)&&(!(&f_startup_counter)))\n"
		"\t\tf_startup_counter <= f_startup_counter + 1;\n"
"\n"
	"\twire	signed	[IWIDTH:0]	f_sumr, f_sumi;\n"
	"\talways @(*)\n"
	"\tbegin\n"
		"\t\tf_sumr = f_dlyleft_r[F_D] + f_dlyright_r[F_D];\n"
		"\t\tf_sumi = f_dlyleft_i[F_D] + f_dlyright_i[F_D];\n"
	"\tend\n"
"\n"
	"\twire	signed	[IWIDTH+CWIDTH-1:0]	f_sumrx, f_sumix;\n"
	"\tassign	f_sumrx = { {(2){f_sumr[IWIDTH]}}, f_sumr, {(CWIDTH-2){1'b0}} };\n"
	"\tassign	f_sumix = { {(2){f_sumi[IWIDTH]}}, f_sumi, {(CWIDTH-2){1'b0}} };\n"
	"\n"
	"\twire	signed	[IWIDTH:0]	f_difr, f_difi;\n"
	"\talways @(*)\n"
	"\tbegin\n"
		"\t\tf_difr = f_dlyleft_r[F_D] - f_dlyright_r[F_D];\n"
		"\t\tf_difi = f_dlyleft_i[F_D] - f_dlyright_i[F_D];\n"
	"\tend\n"
"\n"
	"\twire	signed	[IWIDTH+CWIDTH-1:0]	f_difrx, f_difix;\n"
	"\tassign	f_difrx = { {(CWIDTH){f_difr[IWIDTH]}}, f_difr };\n"
	"\tassign	f_difix = { {(CWIDTH){f_difi[IWIDTH]}}, f_difi };\n"
"\n"
"\n"
	"\talways @(posedge i_clk)\n"
	"\tif (f_startup_counter > F_D)\n"
	"\tbegin\n"
		"\t\tassert(left_sr == f_sumrx);\n"
		"\t\tassert(left_si == f_sumix);\n"
		"\t\tassert(aux_s == f_dlyaux[F_D]);\n"
"\n"
		"\t\tif ((f_difr == 0)&&(f_difi == 0))\n"
		"\t\tbegin\n"
		"\t\t	assert(mpy_r == 0);\n"
		"\t\t	assert(mpy_i == 0);\n"
		"\t\tend else if ((f_dlycoeff_r[F_D] == 0)\n"
		"\t\t		&&(f_dlycoeff_i[F_D] == 0))\n"
		"\t\tbegin\n"
		"\t		assert(mpy_r == 0);\n"
		"\t\t	assert(mpy_i == 0);\n"
		"\t\tend\n"
"\n"
		"\t\tif ((f_dlycoeff_r[F_D] == 1)&&(f_dlycoeff_i[F_D] == 0))\n"
		"\t\tbegin\n"
		"\t\t	assert(mpy_r == f_difrx);\n"
		"\t\t	assert(mpy_i == f_difix);\n"
		"\t\tend\n"
"\n"
		"\t\tif ((f_dlycoeff_r[F_D] == 0)&&(f_dlycoeff_i[F_D] == 1))\n"
		"\t\tbegin\n"
		"\t\t	assert(mpy_r == -f_difix);\n"
		"\t\t	assert(mpy_i ==  f_difrx);\n"
		"\t\tend\n"
"\n"
		"\t\tif ((f_difr == 1)&&(f_difi == 0))\n"
		"\t\tbegin\n"
		"\t\t	assert(mpy_r == f_dlycoeff_r[F_D]);\n"
		"\t\t	assert(mpy_i == f_dlycoeff_i[F_D]);\n"
		"\t\tend\n"
"\n"
		"\t\tif ((f_difr == 0)&&(f_difi == 1))\n"
		"\t\tbegin\n"
		"\t\t	assert(mpy_r == -f_dlycoeff_i[F_D]);\n"
		"\t\t	assert(mpy_i ==  f_dlycoeff_r[F_D]);\n"
		"\t\tend\n"
	"\tend\n"
"\n"
"`endif // FORMAL\n");
	}

	fprintf(fp,
"endmodule\n");

	fclose(fp);
}
