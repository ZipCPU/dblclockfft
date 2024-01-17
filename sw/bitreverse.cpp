////////////////////////////////////////////////////////////////////////////////
//
// Filename: 	bitreverse.cpp
// {{{
// Project:	A General Purpose Pipelined FFT Implementation
//
// Purpose:	
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
#include <sys/stat.h>
#endif

#include <string.h>
#include <string>
#include <math.h>
#include <ctype.h>
#include <assert.h>

#include "defaults.h"
#include "legal.h"
#include "bitreverse.h"

// build_snglbrev(fname, async_reset)
// {{{
void	build_snglbrev(const char *fname, const bool async_reset) {
	FILE	*fp = fopen(fname, "w");
	if (NULL == fp) {
		fprintf(stderr, "Could not open \'%s\' for writing\n", fname);
		perror("O/S Err was:");
		return;
	}

	std::string	resetw("i_reset");
	if (async_reset)
		resetw = std::string("i_areset_n");

	char	*modulename = strdup(fname), *pslash;
	modulename[strlen(modulename)-2] = '\0';
	pslash = strrchr(modulename, '/');
	if (pslash != NULL)
		strcpy(modulename, pslash+1);

	fprintf(fp,
SLASHLINE
"//\n"
"// Filename:\t%s.v\n"
"// {{{\n" // "}}}"
"// Project:\t%s\n"
"//\n"
"// Purpose:\tThis module bitreverses a pipelined FFT input.  It differes\n"
"//		from the dblreverse module in that this is just a simple and\n"
"//	straightforward bitreverse, rather than one written to handle two\n"
"//	words at once.\n"
"//\n"
"//\n%s"
"//\n", modulename, prjname, creator);
	fprintf(fp, "%s", cpyleft);
	fprintf(fp, "//\n//\n`default_nettype\tnone\n//\n");
	fprintf(fp,
"module	%s #(\n"
	"\t\t// {{{\n"
	"\t\tparameter\t\t\tLGSIZE=%d, WIDTH=24\n"
	"\t\t// }}}\n"
	"\t) (\n"
	"\t\t// {{{\n"
	"\t\tinput\twire\t\t\ti_clk, %s, i_ce,\n"
	"\t\tinput\twire\t[(2*WIDTH-1):0]\ti_in,\n"
	"\t\toutput\treg\t[(2*WIDTH-1):0]\to_out,\n"
	"\t\toutput\treg\t\t\to_sync\n"
	"\t\t// }}}\n"
	"\t);\n\n", modulename, TST_DBLREVERSE_LGSIZE,
		resetw.c_str());

	fprintf(fp,
	"\t// Local declarations\n"
	"\t// {{{\n"
	"\treg	[(LGSIZE):0]	wraddr;\n"
	"\twire	[(LGSIZE):0]	rdaddr;\n"
"\n"
	"\treg	[(2*WIDTH-1):0]	brmem	[0:((1<<(LGSIZE+1))-1)];\n"
"\n"
	"\treg	in_reset;\n"
	"\t// }}}\n"
"\n"
	"\t// bitreverse rdaddr\n"
	"\t// {{{\n"
"	genvar	k;\n"
"	generate for(k=0; k<LGSIZE; k=k+1)\n"
"	begin : DBL\n"
"		assign rdaddr[k] = wraddr[LGSIZE-1-k];\n"
"	end endgenerate\n"
"	assign	rdaddr[LGSIZE] = !wraddr[LGSIZE];\n"
	"\t// }}}\n"
"\n"
	"\t// in_reset\n"
	"\t// {{{\n"
	"\tinitial	in_reset = 1'b1;\n");

	if (async_reset)
		fprintf(fp, "\talways @(posedge i_clk, negedge i_areset_n)\n\tif (!i_areset_n)\n");
	else
		fprintf(fp, "\talways @(posedge i_clk)\n\tif (i_reset)\n");
	fprintf(fp,
	"\t	in_reset <= 1'b1;\n"
	"\telse if ((i_ce)&&(&wraddr[(LGSIZE-1):0]))\n"
	"\t	in_reset <= 1'b0;\n"
	"\t// }}}\n"
"\n"
	"\t// wraddr\n"
	"\t// {{{\n"
	"\tinitial	wraddr = 0;\n");

	if (async_reset)
		fprintf(fp, "\talways @(posedge i_clk, negedge i_areset_n)\n\tif (!i_areset_n)\n");
	else
		fprintf(fp, "\talways @(posedge i_clk)\n\tif (i_reset)\n");
	fprintf(fp,
	"\t	wraddr <= 0;\n"
	"\telse if (i_ce)\n"
	"\tbegin\n"
	"\t	brmem[wraddr] <= i_in;\n"
	"\t	wraddr <= wraddr + 1;\n"
	"\tend\n"
	"\t// }}}\n"
"\n"
	"\t// o_out\n"
	"\t// {{{\n"
	"\talways @(posedge i_clk)\n"
	"\tif (i_ce) // If (i_reset) we just output junk ... not a problem\n"
	"\t	o_out <= brmem[rdaddr]; // w/o a sync pulse\n"
	"\t// }}}\n"
"\n"
	"\t// o_sync\n"
	"\t// {{{\n"
	"\tinitial o_sync = 1'b0;\n");

	if (async_reset)
		fprintf(fp, "\talways @(posedge i_clk, negedge i_areset_n)\n\tif (!i_areset_n)\n");
	else
		fprintf(fp, "\talways @(posedge i_clk)\n\tif (i_reset)\n");
	fprintf(fp,
	"\t	o_sync <= 1'b0;\n"
	"\telse if ((i_ce)&&(!in_reset))\n"
	"\t	o_sync <= (wraddr[(LGSIZE-1):0] == 0);\n"
	"\t// }}}\n"
SLASHLINE
SLASHLINE
SLASHLINE
"//\n"
"// Formal property section\n"
"// {{{\n"
SLASHLINE
SLASHLINE
SLASHLINE
"\n");


	if (formal_property_flag) {
		fprintf(fp,
"`ifdef\tFORMAL\n"
"`define\tASSERT	assert\n"
"`ifdef	BITREVERSE\n"
"`define\tASSUME	assume\n");
		if (async_reset)
			fprintf(fp,
		"\n\talways @($global_clock)\n"
			"\t\tassume(i_clk != $past(i_clk));\n\n");

		fprintf(fp,
"`else\n"
"`define\tASSUME	assert\n"
"`endif\n"
"\n"
	"\treg	f_past_valid;\n"
	"\tinitial	f_past_valid = 1'b0;\n"
	"\talways @(posedge i_clk)\n"
		"\t\tf_past_valid <= 1'b1;\n\n");
	
		if (async_reset)
			fprintf(fp,
	"\tinitial	`ASSUME(!i_areset_n);\n"
	"\talways @($global_clock)\n"
	"\tif (!$rose(i_clk)))\n"
		"\t\t`ASSERT(!$rose(i_areset_n));\n\n"
	"\talways @($global_clock)\n"
	"\tif (!$rose(i_clk))\n"
	"\tbegin\n"
		"\t\t`ASSUME($stable(i_ce));\n"
		"\t\t`ASSUME($stable(i_in));\n"
		"\t\t//\n"
		"\t\tif (i_areset_n)\n"
		"\t\tbegin\n"
		"\t\t\t`ASSERT($stable(o_out));\n"
		"\t\t\t`ASSERT($stable(o_sync));\n"
		"\t\tend\n"
	"\tend\n"
"\n"
	"\talways @(posedge i_clk)\n"
		"\tif ((!f_past_valid)||(!i_areset_n))\n"
		"\tbegin\n");
		else
			fprintf(fp,
	"\tinitial	`ASSUME(i_reset);\n"
	"\talways @(posedge i_clk)\n"
		"\tif ((!f_past_valid)||($past(i_reset)))\n"
		"\tbegin\n");

		fprintf(fp,
		"\t\t`ASSERT(wraddr == 0);\n"
		"\t\t`ASSERT(in_reset);\n"
		"\t\t`ASSERT(!o_sync);\n");
		fprintf(fp, "\tend\n");


		fprintf(fp, "`ifdef	BITREVERSE\n"
				"\talways @(posedge i_clk)\n"
				"\t\tassume((i_ce)||($past(i_ce))||($past(i_ce,2)));\n"
				"`endif // BITREVERSE\n\n");

		fprintf(fp,
"\t// Verilator lint_off UNDRIVEN\n"
"\t(* anyconst *) reg	[LGSIZE:0]\tf_const_addr;\n"
"\t// Verilator lint_on  UNDRIVEN\n"
"\twire\t[LGSIZE:0]\tf_reversed_addr;\n"
"\treg\t		f_addr_loaded;\n"
"\treg\t[(2*WIDTH-1):0]\tf_addr_value;\n"
"\n"
"\t// f_reversed_addr\n"
"\t// {{{\n"
"\tgenerate for(k=0; k<LGSIZE; k=k+1)\n"
"\t\tassign\tf_reversed_addr[k] = f_const_addr[LGSIZE-1-k];\n"
"\tendgenerate\n"
"\tassign\tf_reversed_addr[LGSIZE] = f_const_addr[LGSIZE];\n"
"\t// }}}\n"
"\n"
"\t// f_addr_loaded\n"
"\t// {{{\n"
"\tinitial\tf_addr_loaded = 1'b0;\n"
"\talways @(posedge i_clk)\n"
"\tif (i_reset)\n"
"\t\tf_addr_loaded <= 1'b0;\n"
"\telse if (i_ce)\n"
"\tbegin\n"
"\t\tif (wraddr == f_const_addr)\n"
"\t\t\tf_addr_loaded <= 1'b1;\n"
"\t\telse if (rdaddr == f_const_addr)\n"
"\t\t\tf_addr_loaded <= 1'b0;\n"
"\tend\n"
"\t// }}}\n"
"\n"
"\t// f_addr_value\n"
"\t// {{{\n"
"\talways @(posedge i_clk)\n"
"\tif ((i_ce)&&(wraddr == f_const_addr))\n"
"\tbegin\n"
"\t\tf_addr_value <= i_in;\n"
"\t\t`ASSERT(!f_addr_loaded);\n"
"\tend\n"
"\t// }}}\n"
"\n"
"\talways @(posedge i_clk)\n"
"\tif ((f_past_valid)&&(!$past(i_reset))\n"
"\t\t\t&&($past(f_addr_loaded))&&(!f_addr_loaded))\n"
"\t\tassert(o_out == f_addr_value);\n"
"\n"
		"\talways @(*)\n"
		"\tif (o_sync)\n"
			"\t\tassert(wraddr[LGSIZE-1:0] == 1);\n"
"\n"
"\talways @(*)\n"
"\tif ((wraddr[LGSIZE]==f_const_addr[LGSIZE])\n"
"\t\t\t&&(wraddr[LGSIZE-1:0]\n"
"\t\t\t\t\t<= f_const_addr[LGSIZE-1:0]))\n"
"\t\t`ASSERT(!f_addr_loaded);\n"
"\n"
"\talways @(*)\n"
"\tif ((rdaddr[LGSIZE]==f_const_addr[LGSIZE])&&(f_addr_loaded))\n"
"\t\t`ASSERT(wraddr[LGSIZE-1:0]\n"
"\t\t\t\t<= f_reversed_addr[LGSIZE-1:0]+1);\n"
"\n"
"\talways @(*)\n"
"\tif (f_addr_loaded)\n"
"\t\t`ASSERT(brmem[f_const_addr] == f_addr_value);\n"
"\n\n");

		fprintf(fp,
"\t// Make Verilator happy\n"
"\t// {{{\n"
"\t// Verilator lint_off UNUSED\n"
"\twire\tunused_formal;\n"
"\tassign\tunused_formal = &{ 1\'b0, f_reversed_addr[LGSIZE] };\n"
"\t// Verilator lint_on  UNUSED\n"
"\t// }}}\n");

		fprintf(fp,
"`endif\t// FORMAL\n");
	}

	fprintf(fp,
"// }}}\n"
"endmodule\n");

	fclose(fp);
	free(modulename);
}
// }}}

// build_dblreverse(fname, async_reset)
// {{{
void	build_dblreverse(const char *fname, const bool async_reset) {
	FILE	*fp = fopen(fname, "w");
	if (NULL == fp) {
		fprintf(stderr, "Could not open \'%s\' for writing\n", fname);
		perror("O/S Err was:");
		return;
	}

	std::string	resetw("i_reset");
	if (async_reset)
		resetw = std::string("i_areset_n");

	char	*modulename = strdup(fname), *pslash;
	modulename[strlen(modulename)-2] = '\0';
	pslash = strrchr(modulename, '/');
	if (pslash != NULL)
		strcpy(modulename, pslash+1);

	fprintf(fp,
SLASHLINE
"//\n"
"// Filename:\t%s.v\n"
"// {{{\n" // "}}}"
"// Project:\t%s\n"
"//\n"
"// Purpose:\tThis module bitreverses a pipelined FFT input.  Operation is\n"
"//		expected as follows:\n"
"//\n"
"//		i_clk	A running clock at whatever system speed is offered.\n",
	modulename, prjname);

	if (async_reset)
		fprintf(fp,
"//		i_areset_n	An active low asynchronous reset signal,\n"
"//				that resets all internals\n");
	else
		fprintf(fp,
"//		i_reset	A synchronous reset signal, that resets all internals\n");

	fprintf(fp,
"//		i_ce	If this is one, one input is consumed and an output\n"
"//			is produced.\n"
"//		i_in_0, i_in_1\n"
"//			Two inputs to be consumed, each of width WIDTH.\n"
"//		o_out_0, o_out_1\n"
"//			Two of the bitreversed outputs, also of the same\n"
"//			width, WIDTH.  Of course, there is a delay from the\n"
"//			first input to the first output.  For this purpose,\n"
"//			o_sync is present.\n"
"//		o_sync	This will be a 1\'b1 for the first value in any block.\n"
"//			Following a reset, this will only become 1\'b1 once\n"
"//			the data has been loaded and is now valid.  After that,\n"
"//			all outputs will be valid.\n"
"//\n"
"// How do we do bit reversing at two smples per clock?  Can we separate out\n"
"// our work into eight memory banks, writing two banks at once and reading\n"
"// another two banks in the same clock?\n"
"//\n"
"//	mem[00xxx0] = s_0[n]\n"
"//	mem[00xxx1] = s_1[n]\n"
"//	o_0[n] = mem[10xxx0]\n"
"//	o_1[n] = mem[11xxx0]\n"
"//	...\n"
"//	mem[01xxx0] = s_0[m]\n"
"//	mem[01xxx1] = s_1[m]\n"
"//	o_0[m] = mem[10xxx1]\n"
"//	o_1[m] = mem[11xxx1]\n"
"//	...\n"
"//	mem[10xxx0] = s_0[n]\n"
"//	mem[10xxx1] = s_1[n]\n"
"//	o_0[n] = mem[00xxx0]\n"
"//	o_1[n] = mem[01xxx0]\n"
"//	...\n"
"//	mem[11xxx0] = s_0[m]\n"
"//	mem[11xxx1] = s_1[m]\n"
"//	o_0[m] = mem[00xxx1]\n"
"//	o_1[m] = mem[01xxx1]\n"
"//	...\n"
"//\n"
"//	The answer is that, yes we can but: we need to use four memory banks\n"
"//	to do it properly.  These four banks are defined by the two bits\n"
"//	that determine the top and bottom of the correct address.  Larger\n"
"//	FFT\'s would require more memories.\n"
"//\n"
"//\n"
"//	20150602 -- This module has undergone massive rework in order to\n"
"//		ensure that it uses resources efficiently.  As a result,\n"
"//		it now optimizes nicely into block RAMs.  As an unfortunately\n"
"//		side effect, it now passes it\'s bench test (dblrev_tb) but\n"
"//		fails the integration bench test (fft_tb).\n"
"//\n"
"//\n%s"
"//\n", creator);
	fprintf(fp, "%s", cpyleft);
	fprintf(fp, "//\n//\n`default_nettype\tnone\n//\n");
	fprintf(fp,
"module	%s #(\n"
	// (i_clk, %s, i_ce, i_in_0, i_in_1,o_out_0, o_out_1, o_sync);\n
	"\t\t// {{{\n"
	"\t\tparameter\t\t\tLGSIZE=%d, WIDTH=24\n"
	"\t\t// }}}\n"
	"\t) (\n"
	"\t\t// {{{\n"
	"\t\tinput\twire\t\t\ti_clk, %s, i_ce,\n"
	"\t\tinput\twire\t[(2*WIDTH-1):0]\ti_in_0, i_in_1,\n"
	"\t\toutput\twire\t[(2*WIDTH-1):0]\to_out_0, o_out_1,\n"
	"\t\toutput\treg\t\t\to_sync\n"
	"\t\t// }}}\n"
	"\t);\n\n", modulename,
		TST_DBLREVERSE_LGSIZE, resetw.c_str());

	fprintf(fp,
	"\t// Local declarations\n"
	"\t// {{{\n"
	"\treg\t\t\tin_reset;\n"
	"\treg\t[(LGSIZE-1):0]\tiaddr;\n"
	"\twire\t[(LGSIZE-3):0]\tbraddr;\n"
"\n"
	"\treg\t[(2*WIDTH-1):0]\tmem_e [0:((1<<(LGSIZE))-1)];\n"
	"\treg\t[(2*WIDTH-1):0]\tmem_o [0:((1<<(LGSIZE))-1)];\n"
"\n"
	"\treg [(2*WIDTH-1):0] evn_out_0, evn_out_1, odd_out_0, odd_out_1;\n"
	"\treg\tadrz;\n"
	"\t// }}}\n"
"\n"
	"\t// braddr\n"
	"\t// {{{\n"
	"\tgenvar\tk;\n"
	"\tgenerate for(k=0; k<LGSIZE-2; k=k+1)\n"
	"\tbegin : gen_a_bit_reversed_value\n"
		"\t\tassign braddr[k] = iaddr[LGSIZE-3-k];\n"
	"\tend endgenerate\n"
	"\t// }}}\n"
"\n"
	"\t// iaddr, in_reset, o_sync\n"
	"\t// {{{\n"
	"\tinitial iaddr = 0;\n"
	"\tinitial in_reset = 1\'b1;\n"
	"\tinitial o_sync = 1\'b0;\n");

	if (async_reset)
		fprintf(fp, "\talways @(posedge i_clk, negedge i_areset_n)\n\tif (!i_areset_n)\n");
	else
		fprintf(fp, "\talways @(posedge i_clk)\n\tif (i_reset)\n");
	fprintf(fp,
		"\tbegin\n"
			"\t\tiaddr <= 0;\n"
			"\t\tin_reset <= 1\'b1;\n"
			"\t\to_sync <= 1\'b0;\n"
		"\tend else if (i_ce)\n"
		"\tbegin\n"
			"\t\tiaddr <= iaddr + { {(LGSIZE-1){1\'b0}}, 1\'b1 };\n"
			"\t\tif (&iaddr[(LGSIZE-2):0])\n"
				"\t\t\tin_reset <= 1\'b0;\n"
			"\t\tif (in_reset)\n"
				"\t\t\to_sync <= 1\'b0;\n"
			"\t\telse\n"
				"\t\t\to_sync <= ~(|iaddr[(LGSIZE-2):0]);\n"
		"\tend\n"
	"\t// }}}\n"
"\n"
	"\t// Write to memories mem_e and mem_o\n"
	"\t// {{{\n"
	"\talways @(posedge i_clk)\n"
	"\tif (i_ce)\n"
		"\t\tmem_e[iaddr] <= i_in_0;\n\n"
	"\talways @(posedge i_clk)\n"
	"\tif (i_ce)\n"
		"\t\tmem_o[iaddr] <= i_in_1;\n"
	"\t// }}}\n"
"\n"
	"\t// Read from memories into: [evn|odd]_out_[0|1]\n"
	"\t// {{{\n"
	"\talways @(posedge i_clk)\n"
	"\tif (i_ce)\n"
	"\t\tevn_out_0 <= mem_e[{!iaddr[LGSIZE-1],1\'b0,braddr}];\n"
"\n"
	"\talways @(posedge i_clk)\n"
	"\tif (i_ce)\n"
	"\t\tevn_out_1 <= mem_e[{!iaddr[LGSIZE-1],1\'b1,braddr}];\n"
"\n"
	"\talways @(posedge i_clk)\n"
	"\tif (i_ce)\n"
	"\t\todd_out_0 <= mem_o[{!iaddr[LGSIZE-1],1\'b0,braddr}];\n"
"\n"
	"\talways @(posedge i_clk)\n"
	"\tif (i_ce)\n"
	"\t\todd_out_1 <= mem_o[{!iaddr[LGSIZE-1],1\'b1,braddr}];\n"
	"\t// }}}\n"
"\n"
	"\t// adrz\n"
	"\t// {{{\n"
	"\talways @(posedge i_clk)\n"
	"\tif (i_ce)\n"
		"\t\tadrz <= iaddr[LGSIZE-2];\n"
	"\t// }}}\n"
"\n"
	"\tassign\to_out_0 = (adrz)?odd_out_0:evn_out_0;\n"
	"\tassign\to_out_1 = (adrz)?odd_out_1:evn_out_1;\n"
"\n");

	fprintf(fp,
SLASHLINE
SLASHLINE
SLASHLINE
"//\n"
"// Formal property section\n"
"// {{{\n"
SLASHLINE
SLASHLINE
SLASHLINE );
	if (formal_property_flag) {
		fprintf(fp,
"`ifdef\tFORMAL\n"
	"\t// Formal declarations\n"
	"\t// {{{\n"
"`define\tASSERT	assert\n"
"`ifdef	BITREVERSE\n"
"`define\tASSUME	assume\n");
		if (async_reset)
			fprintf(fp,
		"\n\talways @($global_clock)\n"
			"\t\tassume(i_clk != $past(i_clk));\n\n");

		fprintf(fp,
"`else\n"
"`define\tASSUME	assert\n"
"`endif\n"
"\n"
	"\treg	f_past_valid;\n"
	"\t(* anyconst *) reg	[LGSIZE-1:0]	f_const_addr;\n"
	"\twire	[LGSIZE-3:0]	f_reversed_addr;\n"
	"\t// reg	[LGSIZE:0]	f_now;\n"
	"\treg			f_addr_loaded_0, f_addr_loaded_1;\n"
	"\treg	[(2*WIDTH-1):0]	f_data_0, f_data_1;\n"
	"\twire			f_writing, f_reading;\n"
	"\t// }}}\n"
"\n");

	fprintf(fp,
	"\tinitial	f_past_valid = 1'b0;\n"
	"\talways @(posedge i_clk)\n"
		"\t\tf_past_valid <= 1'b1;\n\n");
	
	if (async_reset)
		fprintf(fp,
	"\tinitial	`ASSUME(!i_areset_n);\n"
	"\talways @($global_clock)\n"
	"\tif (!$rose(i_clk)))\n"
		"\t\t`ASSERT(!$rose(i_areset_n));\n\n"
	"\talways @($global_clock)\n"
	"\tif (!$rose(i_clk))\n"
	"\tbegin\n"
		"\t\t`ASSUME($stable(i_ce));\n"
		"\t\t`ASSUME($stable(i_in_0));\n"
		"\t\t`ASSUME($stable(i_in_1));\n"
		"\t\t//\n"
		"\t\tif (i_areset_n)\n"
		"\t\tbegin\n"
		"\t\t\t`ASSERT($stable(o_out_0));\n"
		"\t\t\t`ASSERT($stable(o_out_1));\n"
		"\t\t\t`ASSERT($stable(o_sync));\n"
		"\t\tend\n"
	"\tend\n"
"\n"
	"\talways @(posedge i_clk)\n"
	"\tif ((!f_past_valid)||(!i_areset_n))\n"
	"\tbegin\n");
	else
		fprintf(fp,
	"\tinitial	`ASSUME(i_reset);\n"
	"\talways @(posedge i_clk)\n"
		"\tif ((!f_past_valid)||($past(i_reset)))\n"
		"\tbegin\n");

	fprintf(fp,
		"\t\t`ASSERT(iaddr == 0);\n"
		"\t\t`ASSERT(in_reset);\n"
		"\t\t`ASSERT(!o_sync);\n");
	fprintf(fp, "\tend\n");


		fprintf(fp, "`ifdef	BITREVERSE\n"
				"\talways @(posedge i_clk)\n"
				"\t\tassume((i_ce)||($past(i_ce))||($past(i_ce,2)));\n"
				"`endif // BITREVERSE\n\n");


		fprintf(fp,
	"\t// f_reversed_addr\n"
	"\t// {{{\n"
	"\tgenerate for(k=0; k<LGSIZE-2; k=k+1)\n"
	"\t	assign	f_reversed_addr[k] = f_const_addr[LGSIZE-3-k];\n"
	"\tendgenerate\n"
	"\t// }}}\n"
"\n"
	"\t\tassign	f_writing=(f_const_addr[LGSIZE-1]==iaddr[LGSIZE-1]);\n"
	"\t\tassign	f_reading=(f_const_addr[LGSIZE-1]!=iaddr[LGSIZE-1]);\n");


	fprintf(fp,
	"\t// f_addr_loaded_[0|1]\n"
	"\t// {{{\n"
	"\tinitial	f_addr_loaded_0 = 1'b0;\n"
	"\tinitial	f_addr_loaded_1 = 1'b0;\n"
	"\talways @(posedge i_clk)\n"
	"\tif (i_reset)\n"
	"\tbegin\n"
	"\t	f_addr_loaded_0 <= 1'b0;\n"
	"\t	f_addr_loaded_1 <= 1'b0;\n"
	"\tend else if (i_ce)\n"
	"\tbegin\n"
	"\t	if (iaddr == f_const_addr)\n"
	"\t	begin\n"
	"\t		f_addr_loaded_0 <= 1'b1;\n"
	"\t		f_addr_loaded_1 <= 1'b1;\n"
	"\t	end\n"
"\n"
	"\t	if (f_reading)\n"
	"\t	begin\n"
	"\t		if ((braddr == f_const_addr[LGSIZE-3:0])\n"
	"\t			&&(iaddr[LGSIZE-2] == 1'b0))\n"
	"\t			f_addr_loaded_0 <= 1'b0;\n"
"\n"
	"\t		if ((braddr == f_const_addr[LGSIZE-3:0])\n"
	"\t			&&(iaddr[LGSIZE-2] == 1'b1))\n"
	"\t			f_addr_loaded_1 <= 1'b0;\n"
	"\t	end\n"
	"\tend\n"
	"\t// }}}\n"
"\n"
	"\t// f_data_0, f_data_1\n"
	"\t// {{{\n"
	"\talways @(posedge i_clk)\n"
	"\tif ((i_ce)&&(iaddr == f_const_addr))\n"
	"\tbegin\n"
	"\t	f_data_0 <= i_in_0;\n"
	"\t	f_data_1 <= i_in_1;\n"
	"\t	`ASSERT(!f_addr_loaded_0);\n"
	"\t	`ASSERT(!f_addr_loaded_1);\n"
	"\tend\n"
	"\t// }}}\n"
"\n"
	"\talways @(posedge i_clk)\n"
	"\tif ((f_past_valid)&&(!$past(i_reset))\n"
	"\t		&&($past(f_addr_loaded_0))&&(!f_addr_loaded_0))\n"
	"\tbegin\n"
	"\t	assert(!$past(iaddr[LGSIZE-2]));\n"
	"\t	if (f_const_addr[LGSIZE-2])\n"
	"\t		assert(o_out_1 == f_data_0);\n"
	"\t	else\n"
	"\t		assert(o_out_0 == f_data_0);\n"
	"\tend\n"
"\n"
	"\talways @(posedge i_clk)\n"
	"\tif ((f_past_valid)&&(!$past(i_reset))\n"
	"\t		&&($past(f_addr_loaded_1))&&(!f_addr_loaded_1))\n"
	"\tbegin\n"
	"\t	assert($past(iaddr[LGSIZE-2]));\n"
	"\t	if (f_const_addr[LGSIZE-2])\n"
	"\t		assert(o_out_1 == f_data_1);\n"
	"\t	else\n"
	"\t		assert(o_out_0 == f_data_1);\n"
	"\tend\n"
"\n"
	"\talways @(*)\n"
	"\t	`ASSERT(o_sync == ((iaddr[LGSIZE-2:0] == 1)&&(!in_reset)));\n"
"\n"
	"\t// Before writing to a section, the loaded flags should be\n"
	"\t// zero\n"
	"\talways @(*)\n"
	"\tif (f_writing)\n"
	"\tbegin\n"
	"\t	`ASSERT(f_addr_loaded_0 == (iaddr[LGSIZE-2:0]\n"
	"\t				> f_const_addr[LGSIZE-2:0]));\n"
	"\t	`ASSERT(f_addr_loaded_1 == (iaddr[LGSIZE-2:0]\n"
	"\t				> f_const_addr[LGSIZE-2:0]));\n"
	"\tend\n"
"\n"
	"\t// If we were writing, and now we are reading, then both\n"
	"\t// f_addr_loaded flags must be set\n"
	"\talways @(posedge i_clk)\n"
	"\tif ((f_past_valid)&&(!$past(i_reset))\n"
	"\t		&&($past(f_writing))&&(f_reading))\n"
	"\tbegin\n"
	"\t	`ASSERT(f_addr_loaded_0);\n"
	"\t	`ASSERT(f_addr_loaded_1);\n"
	"\tend\n"
"\n"
	"\talways @(*)\n"
	"\tif (f_writing)\n"
	"\t	`ASSERT(f_addr_loaded_0 == f_addr_loaded_1);\n"
"\n"
	"\t// When reading, and the loaded flag is zero, our pointer\n"
	"\t// must not have hit the address of interest yet\n"
	"\talways @(*)\n"
	"\tif ((!in_reset)&&(f_reading))\n"
	"\t	`ASSERT(f_addr_loaded_0 ==\n"
	"\t		((!iaddr[LGSIZE-2])&&(iaddr[LGSIZE-3:0]\n"
	"\t			<= f_reversed_addr[LGSIZE-3:0])));\n"
"\n"
	"\talways @(*)\n"
	"\tif ((!in_reset)&&(f_reading))\n"
	"\t	`ASSERT(f_addr_loaded_1 ==\n"
	"\t		((!iaddr[LGSIZE-2])||(iaddr[LGSIZE-3:0]\n"
	"\t			<= f_reversed_addr[LGSIZE-3:0])));\n"
"\n"
	"\talways @(*)\n"
	"\tif ((in_reset)&&(f_reading))\n"
	"\tbegin\n"
	"\t	`ASSERT(!f_addr_loaded_0);\n"
	"\t	`ASSERT(!f_addr_loaded_1);\n"
	"\tend\n"
"\n"
	"\talways @(*)\n"
	"\tif(iaddr[LGSIZE-1])\n"
	"\t	`ASSERT(!in_reset);\n"
"\n"
	"\talways @(*)\n"
	"\tif (f_addr_loaded_0)\n"
	"\t	`ASSERT(mem_e[f_const_addr] == f_data_0);\n"
	"\talways @(*)\n"
	"\tif (f_addr_loaded_1)\n"
	"\t	`ASSERT(mem_o[f_const_addr] == f_data_1);\n"
"\n");


		fprintf(fp,
"`endif\t// FORMAL\n");
	} else
	fprintf(fp, "// Formal properties have not included in this build\n");


	fprintf(fp,
"// }}}\n"
"endmodule\n");

	fclose(fp);
	free(modulename);
}
// }}}
