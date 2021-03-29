////////////////////////////////////////////////////////////////////////////////
//
// Filename: 	rounding.cpp
// {{{
// Project:	A General Purpose Pipelined FFT Implementation
//
// Purpose:	To create one of a series of modules to handle dropping bits
// 		within the FFT implementation.
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
#define	_CRT_SECURE_NO_WARNINGS	// ms vs 2012 doesn't like fopen

#include <stdio.h>
#include <stdlib.h>

#include <string.h>
#include <string>
#include <math.h>
#include <ctype.h>
#include <assert.h>

#include "legal.h"
#include "rounding.h"

#define	SLASHLINE "////////////////////////////////////////////////////////////////////////////////\n"

// build_truncator
// {{{
void	build_truncator(const char *fname) {
	printf("TRUNCATING!\n");
	FILE	*fp = fopen(fname, "w");
	if (NULL == fp) {
		fprintf(stderr, "Could not open \'%s\' for writing\n", fname);
		perror("O/S Err was:");
		return;
	}

	fprintf(fp,
SLASHLINE
"//\n"
"// Filename:\ttruncate.v\n"
"// {{{\n" // "}}}"
"// Project:\t%s\n"
"//\n"
"// Purpose:	Truncation is one of several options that can be used\n"
"//		internal to the various FFT stages to drop bits from one\n"
"//	stage to the next.  In general, it is the simplest method of dropping\n"
"//	bits, since it requires only a bit selection.\n"
"//\n"
"//	This form of rounding isn\'t really that great for FFT\'s, since it\n"
"//	tends to produce a DC bias in the result.  (Other less pronounced\n"
"//	biases may also exist.)\n"
"//\n"
"//	This particular version also registers the output with the clock, so\n"
"//	there will be a delay of one going through this module.  This will\n"
"//	keep it in line with the other forms of rounding that can be used.\n"
"//\n"
"//\n%s"
"//\n",
		prjname, creator);

	fprintf(fp, "%s", cpyleft);
	fprintf(fp, "//\n//\n`default_nettype\tnone\n//\n");
	fprintf(fp,
"module	truncate(i_clk, i_ce, i_val, o_val);\n"
	"\tparameter\tIWID=16, OWID=8, SHIFT=0;\n"
	"\tinput\twire\t\t\t\ti_clk, i_ce;\n"
	"\tinput\twire\tsigned\t[(IWID-1):0]\ti_val;\n"
	"\toutput\treg\tsigned\t[(OWID-1):0]\to_val;\n"
"\n"
	"\talways @(posedge i_clk)\n"
		"\t\tif (i_ce)\n"
		"\t\t\to_val <= i_val[(IWID-1-SHIFT):(IWID-SHIFT-OWID)];\n"
"\n"
"endmodule\n");
}
// }}}

// build_roundhalfup
// {{{
void	build_roundhalfup(const char *fname) {
	FILE	*fp = fopen(fname, "w");
	if (NULL == fp) {
		fprintf(stderr, "Could not open \'%s\' for writing\n", fname);
		perror("O/S Err was:");
		return;
	}

	fprintf(fp,
SLASHLINE
"//\n"
"// Filename:\troundhalfup.v\n"
"// {{{\n" // "}}}"
"// Project:\t%s\n"
"//\n"
"// Purpose:\tRounding half up is the way I was always taught to round in\n"
"//		school.  A one half value is added to the result, and then\n"
"//	the result is truncated.  When used in an FFT, this produces less\n"
"//	bias than the truncation method, although a bias still tends to\n"
"//	remain.\n"
"//\n"
"//\n%s"
"//\n",
		prjname, creator);

	fprintf(fp, "%s", cpyleft);
	fprintf(fp, "//\n//\n`default_nettype\tnone\n//\n");
	fprintf(fp,
"module	roundhalfup(i_clk, i_ce, i_val, o_val);\n"
	"\tparameter\tIWID=16, OWID=8, SHIFT=0;\n"
	"\tinput\twire\t\t\t\ti_clk, i_ce;\n"
	"\tinput\twire\tsigned\t[(IWID-1):0]\ti_val;\n"
	"\toutput\treg\tsigned\t[(OWID-1):0]\to_val;\n"
"\n"
	"\t// Let's deal with two cases to be as general as we can be here\n"
	"\t//\n"
	"\t//	1. The desired output would lose no bits at all\n"
	"\t//	2. One or more bits would be dropped, so the rounding is simply\n"
	"\t//\t\ta matter of adding one to the bit about to be dropped,\n"
	"\t//\t\tmoving all halfway and above numbers up to the next\n"
	"\t//\t\tvalue.\n"
	"\tgenerate\n"
	"\tif (IWID-SHIFT == OWID)\n"
	"\tbegin : NO_ROUNDING // No truncation or rounding, output drops no bits\n"
"\n"
		"\t\talways @(posedge i_clk)\n"
			"\t\t\tif (i_ce)\to_val <= i_val[(IWID-SHIFT-1):0];\n"
"\n"
	"\tend else // if (IWID-SHIFT-1 >= OWID)\n"
	"\tbegin : DROP_ONE_BIT // Output drops one bit, can only add one or ... not.\n"
		"\t\twire\t[(OWID-1):0]	truncated_value, rounded_up;\n"
		"\t\twire\t\t\tlast_valid_bit, first_lost_bit;\n"
		"\t\tassign\ttruncated_value=i_val[(IWID-1-SHIFT):(IWID-SHIFT-OWID)];\n"
		"\t\tassign\trounded_up=truncated_value + {{(OWID-1){1\'b0}}, 1\'b1 };\n"
		"\t\tassign\tfirst_lost_bit = i_val[(IWID-SHIFT-OWID-1)];\n"
"\n"
		"\t\talways @(posedge i_clk)\n"
		"\t\t\tif (i_ce)\n"
		"\t\t\tbegin\n"
			"\t\t\t\tif (!first_lost_bit) // Round down / truncate\n"
			"\t\t\t\t\to_val <= truncated_value;\n"
			"\t\t\t\telse\n"
			"\t\t\t\t\to_val <= rounded_up; // even value\n"
		"\t\t\tend\n"
"\n"
	"\tend\n"
	"\tendgenerate\n"
"\n"
"endmodule\n");
}
// }}}

// build_roundfromzero
// {{{
void	build_roundfromzero(const char *fname) {
	FILE	*fp = fopen(fname, "w");
	if (NULL == fp) {
		fprintf(stderr, "Could not open \'%s\' for writing\n", fname);
		perror("O/S Err was:");
		return;
	}

	fprintf(fp,
SLASHLINE
"//\n"
"// Filename:\troundfromzero.v\n"
"// {{{\n" // "}}}"
"// Project:	%s\n"
"//\n"
"// Purpose:	Truncation is one of several options that can be used\n"
"//		internal to the various FFT stages to drop bits from one\n"
"//	stage to the next.  In general, it is the simplest method of dropping\n"
"//	bits, since it requires only a bit selection.\n"
"//\n"
"//	This form of rounding isn\'t really that great for FFT\'s, since it\n"
"//	tends to produce a DC bias in the result.  (Other less pronounced\n"
"//	biases may also exist.)\n"
"//\n"
"//	This particular version also registers the output with the clock, so\n"
"//	clock, so there will be a delay of one going through this module.\n"
"//	This will keep it in line with the other forms of rounding that can\n"
"//	be used.\n"
"//\n"
"//\n%s"
"//\n",
		prjname, creator);

	fprintf(fp, "%s", cpyleft);
	fprintf(fp, "//\n//\n`default_nettype\tnone\n//\n");
	fprintf(fp,
"module	roundfromzero(i_clk, i_ce, i_val, o_val);\n"
	"\tparameter\tIWID=16, OWID=8, SHIFT=0;\n"
	"\tinput\twire\t\t\t\ti_clk, i_ce;\n"
	"\tinput\twire\tsigned\t[(IWID-1):0]\ti_val;\n"
	"\toutput\treg\tsigned\t[(OWID-1):0]\to_val;\n"
"\n"
	"\t// Let's deal with three cases to be as general as we can be here\n"
	"\t//\n"
	"\t//\t1. The desired output would lose no bits at all\n"
	"\t//\t2. One bit would be dropped, so the rounding is simply\n"
	"\t//\t\tadjusting the value to be the closer to zero in\n"
	"\t//\t\tcases of being halfway between two.  If identically\n"
	"\t//\t\tequal to a number, we just leave it as is.\n"
	"\t//\t3. Two or more bits would be dropped.  In this case, we round\n"
	"\t//\t\tnormally unless we are rounding a value of exactly\n"
	"\t//\t\thalfway between the two.  In the halfway case, we\n"
	"\t//\t\tround away from zero.\n"
	"\tgenerate\n"
	"\tif (IWID == OWID)\n"
	"\tbegin : NO_ROUNDING\n"
		"\t\t// In this case, the shift is irrelevant and cannot be\n"
		"\t\t// applied. No truncation or rounding takes place here.\n"
"\n"
		"\t\talways @(posedge i_clk)\n"
		"\t\tif (i_ce)\to_val <= i_val[(IWID-1):0];\n"
"\n"
	"\tend else if (IWID-SHIFT == OWID)\n"
	"\tbegin : SHIFT_ONE_BIT\n"
	"\t\t// No truncation or rounding, output drops no bits\n"
"\n"
		"\t\talways @(posedge i_clk)\n"
		"\t\tif (i_ce)\to_val <= i_val[(IWID-SHIFT-1):0];\n"
"\n"
	"\tend else if (IWID-SHIFT-1 == OWID)\n"
	"\tbegin : DROP_ONE_BIT\n"
	"\t\t// Output drops one bit, can only add one or ... not.\n"
	"\n"
	"\t\twire\t[(OWID-1):0]\ttruncated_value, rounded_up;\n"
	"\t\twire\t\t\tsign_bit, first_lost_bit;\n"
	"\t\tassign\ttruncated_value=i_val[(IWID-1-SHIFT):(IWID-SHIFT-OWID)];\n"
	"\t\tassign\trounded_up=truncated_value + {{(OWID-1){1\'b0}}, 1\'b1 };\n"
	"\t\tassign\tfirst_lost_bit = i_val[0];\n"
	"\t\tassign\tsign_bit = i_val[(IWID-1)];\n"
"\n"
	"\t\talways @(posedge i_clk)\n"
	"\t\tif (i_ce)\n"
	"\t\tbegin\n"
		"\t\t\tif (!first_lost_bit) // Round down / truncate\n"
			"\t\t\t\to_val <= truncated_value;\n"
		"\t\t\telse if (sign_bit)\n"
			"\t\t\t\to_val <= truncated_value;\n"
		"\t\t\telse\n"
			"\t\t\t\to_val <= rounded_up;\n"
	"\t\tend\n"
"\n"
	"\tend else begin : ROUND_RESULT\n"
		"\t\t// If there's more than one bit we are dropping\n"
		"\t\twire\t[(OWID-1):0]\ttruncated_value, rounded_up;\n"
		"\t\twire\t\t\tsign_bit, first_lost_bit;\n"
		"\t\tassign\ttruncated_value=i_val[(IWID-1-SHIFT):(IWID-SHIFT-OWID)];\n"
		"\t\tassign\trounded_up=truncated_value + {{(OWID-1){1\'b0}}, 1\'b1 };\n"
		"\t\tassign\tfirst_lost_bit = i_val[(IWID-SHIFT-OWID-1)];\n"
		"\t\tassign\tsign_bit = i_val[(IWID-1)];\n"
"\n"
		"\t\twire\t[(IWID-SHIFT-OWID-2):0]\tother_lost_bits;\n"
		"\t\tassign\tother_lost_bits = i_val[(IWID-SHIFT-OWID-2):0];\n"
"\n"
		"\t\talways @(posedge i_clk)\n"
		"\t\tif (i_ce)\n"
		"\t\tbegin\n"
			"\t\t\tif (!first_lost_bit) // Round down / truncate\n"
				"\t\t\t\to_val <= truncated_value;\n"
			"\t\t\telse if (|other_lost_bits) // Round up to\n"
				"\t\t\t\to_val <= rounded_up; // closest value\n"
			"\t\t\telse if (sign_bit)\n"
				"\t\t\t\to_val <= truncated_value;\n"
			"\t\t\telse\n"
				"\t\t\t\to_val <= rounded_up;\n"
		"\t\tend\n"
	"\tend\n"
	"\tendgenerate\n"
"\n"
"endmodule\n");
}
// }}}

// build_convround
// {{{
void	build_convround(const char *fname) {
	FILE	*fp = fopen(fname, "w");
	if (NULL == fp) {
		fprintf(stderr, "Could not open \'%s\' for writing\n", fname);
		perror("O/S Err was:");
		return;
	}

	fprintf(fp,
SLASHLINE
"//\n"
"// Filename: 	convround.v\n"
"// {{{\n" // "}}}"
"// Project:	%s\n"
"//\n"
"// Purpose:	A convergent rounding routine, also known as banker\'s\n"
"//		rounding, Dutch rounding, Gaussian rounding, unbiased\n"
"//	rounding, or ... more, at least according to Wikipedia.\n"
"//\n"
"//	This form of rounding works by rounding, when the direction is in\n"
"//	question, towards the nearest even value.\n"
"//\n"
"//\n%s"
"//\n",
		prjname, creator);

	fprintf(fp, "%s", cpyleft);
	fprintf(fp, "//\n//\n`default_nettype\tnone\n//\n");
	fprintf(fp,
"module	convround(i_clk, i_ce, i_val, o_val);\n"
"\tparameter\tIWID=16, OWID=8, SHIFT=0;\n"
"\tinput\twire\t\t\t\ti_clk, i_ce;\n"
"\tinput\twire\tsigned\t[(IWID-1):0]\ti_val;\n"
"\toutput\treg\tsigned\t[(OWID-1):0]\to_val;\n"
"\n"
"\t// Let's deal with three cases to be as general as we can be here\n"
"\t//\n"
"\t//\t1. The desired output would lose no bits at all\n"
"\t//\t2. One bit would be dropped, so the rounding is simply\n"
"\t//\t\tadjusting the value to be the nearest even number in\n"
"\t//\t\tcases of being halfway between two.  If identically\n"
"\t//\t\tequal to a number, we just leave it as is.\n"
"\t//\t3. Two or more bits would be dropped.  In this case, we round\n"
"\t//\t\tnormally unless we are rounding a value of exactly\n"
"\t//\t\thalfway between the two.  In the halfway case we round\n"
"\t//\t\tto the nearest even number.\n"
"\tgenerate\n"
// What if IWID < OWID?  We should expand here ... somehow
	"\tif (IWID == OWID) // In this case, the shift is irrelevant and\n"
	"\tbegin : NO_ROUNDING // cannot be applied.  No truncation or rounding takes\n"
	"\t// effect here.\n"
"\n"
		"\t\talways @(posedge i_clk)\n"
		"\t\tif (i_ce)\to_val <= i_val[(IWID-1):0];\n"
"\n"
// What if IWID-SHIFT < OWID?  Shouldn't we also shift here as well?
"\tend else if (IWID-SHIFT < OWID)\n"
"\tbegin : ADD_BITS_TO_OUTPUT // No truncation or rounding, output drops no bits\n"
"\t// Instead, we need to stuff the bits in the output\n"
"\n"
"\t\talways @(posedge i_clk)\n"
"\t\tif (i_ce)\to_val <= { {(OWID-IWID+SHIFT){i_val[IWID-SHIFT-1]}}, i_val[(IWID-SHIFT-1):0] };\n"
"\n"
"\tend else if (IWID-SHIFT == OWID)\n"
"\tbegin : SHIFT_ONE_BIT\n"
"\t// No truncation or rounding, output drops no bits\n"
"\n"
"\t\talways @(posedge i_clk)\n"
"\t\tif (i_ce)\to_val <= i_val[(IWID-SHIFT-1):0];\n"
"\n"
"\tend else if (IWID-SHIFT-1 == OWID)\n"
// Is there any way to limit the number of bits that are examined here, for the
// purpose of simplifying/reducing logic?  I mean, if we go from 32 to 16 bits,
// must we check all 15 bits for equality to zero?
"\tbegin : DROP_ONE_BIT // Output drops one bit, can only add one or ... not.\n"
"\t\twire\t[(OWID-1):0]	truncated_value, rounded_up;\n"
"\t\twire\t\t\tlast_valid_bit, first_lost_bit;\n"
"\t\tassign\ttruncated_value=i_val[(IWID-1-SHIFT):(IWID-SHIFT-OWID)];\n"
"\t\tassign\trounded_up=truncated_value + {{(OWID-1){1\'b0}}, 1\'b1 };\n"
"\t\tassign\tlast_valid_bit = truncated_value[0];\n"
"\t\tassign\tfirst_lost_bit = i_val[0];\n"
"\n"
"\t\talways @(posedge i_clk)\n"
"\t\tif (i_ce)\n"
"\t\tbegin\n"
"\t\t\tif (!first_lost_bit) // Round down / truncate\n"
"\t\t\t\to_val <= truncated_value;\n"
"\t\t\telse if (last_valid_bit)// Round up to nearest\n"
"\t\t\t\to_val <= rounded_up; // even value\n"
"\t\t\telse // else round down to the nearest\n"
"\t\t\t\to_val <= truncated_value; // even value\n"
"\t\tend\n"
"\n"
"\tend else // If there's more than one bit we are dropping\n"
"\tbegin : ROUND_RESULT\n"
"\t\twire\t[(OWID-1):0]	truncated_value, rounded_up;\n"
"\t\twire\t\t\tlast_valid_bit, first_lost_bit;\n\n"
"\t\tassign\ttruncated_value=i_val[(IWID-1-SHIFT):(IWID-SHIFT-OWID)];\n"
"\t\tassign\trounded_up=truncated_value + {{(OWID-1){1\'b0}}, 1\'b1 };\n"
"\t\tassign\tlast_valid_bit = truncated_value[0];\n"
"\t\tassign\tfirst_lost_bit = i_val[(IWID-SHIFT-OWID-1)];\n"
"\n"
"\t\twire\t[(IWID-SHIFT-OWID-2):0]\tother_lost_bits;\n"
"\t\tassign\tother_lost_bits = i_val[(IWID-SHIFT-OWID-2):0];\n"
"\n"
"\t\talways @(posedge i_clk)\n"
"\t\t\tif (i_ce)\n"
"\t\t\tbegin\n"
"\t\t\t\tif (!first_lost_bit) // Round down / truncate\n"
"\t\t\t\t\to_val <= truncated_value;\n"
"\t\t\t\telse if (|other_lost_bits) // Round up to\n"
"\t\t\t\t\to_val <= rounded_up; // closest value\n"
"\t\t\t\telse if (last_valid_bit) // Round up to\n"
"\t\t\t\t\to_val <= rounded_up; // nearest even\n"
"\t\t\t\telse	// else round down to nearest even\n"
"\t\t\t\t\to_val <= truncated_value;\n"
"\t\t\tend\n"
"\tend\n"
"\tendgenerate\n"
"\n"
"endmodule\n");
}
// }}}

