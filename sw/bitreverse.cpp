////////////////////////////////////////////////////////////////////////////////
//
// Filename: 	bitreverse.cpp
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
// with this program.  (It's in the $(ROOT)/doc directory.  Run make with no
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
"//\n"
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
"module	%s(i_clk, %s, i_ce, i_in, o_out, o_sync);\n"
	"\tparameter\t\t\tLGSIZE=%d, WIDTH=24;\n"
	"\tinput\twire\t\t\ti_clk, %s, i_ce;\n"
	"\tinput\twire\t[(2*WIDTH-1):0]\ti_in;\n"
	"\toutput\treg\t[(2*WIDTH-1):0]\to_out;\n"
	"\toutput\treg\t\t\to_sync;\n", modulename, resetw.c_str(),
		TST_DBLREVERSE_LGSIZE,
		resetw.c_str());

	fprintf(fp,
"	reg	[(LGSIZE):0]	wraddr;\n"
"	wire	[(LGSIZE):0]	rdaddr;\n"
"\n"
"	reg	[(2*WIDTH-1):0]	brmem	[0:((1<<(LGSIZE+1))-1)];\n"
"\n"
"	genvar	k;\n"
"	generate for(k=0; k<LGSIZE; k=k+1)\n"
"		assign rdaddr[k] = wraddr[LGSIZE-1-k];\n"
"	endgenerate\n"
"	assign	rdaddr[LGSIZE] = !wraddr[LGSIZE];\n"
"\n"
"	reg	in_reset;\n"
"\n"
"	initial	in_reset = 1'b1;\n");

	if (async_reset)
		fprintf(fp, "\talways @(posedge i_clk, negedge i_areset_n)\n\t\tif (!i_areset_n)\n");
	else
		fprintf(fp, "\talways @(posedge i_clk)\n\t\tif (i_reset)\n");
	fprintf(fp,
"			in_reset <= 1'b1;\n"
"		else if ((i_ce)&&(&wraddr[(LGSIZE-1):0]))\n"
"			in_reset <= 1'b0;\n"
"\n"
"	initial	wraddr = 0;\n");

	if (async_reset)
		fprintf(fp, "\talways @(posedge i_clk, negedge i_areset_n)\n\t\tif (!i_areset_n)\n");
	else
		fprintf(fp, "\talways @(posedge i_clk)\n\t\tif (i_reset)\n");
	fprintf(fp,
"			wraddr <= 0;\n"
"		else if (i_ce)\n"
"		begin\n"
"			brmem[wraddr] <= i_in;\n"
"			wraddr <= wraddr + 1;\n"
"		end\n"
"\n"
"	always @(posedge i_clk)\n"
"		if (i_ce) // If (i_reset) we just output junk ... not a problem\n"
"			o_out <= brmem[rdaddr]; // w/o a sync pulse\n"
"\n"
"	initial	o_sync = 1'b0;\n");

	if (async_reset)
		fprintf(fp, "\talways @(posedge i_clk, negedge i_areset_n)\n\t\tif (!i_areset_n)\n");
	else
		fprintf(fp, "\talways @(posedge i_clk)\n\t\tif (i_reset)\n");
	fprintf(fp,
"			o_sync <= 1'b0;\n"
"		else if ((i_ce)&&(!in_reset))\n"
"			o_sync <= (wraddr[(LGSIZE-1):0] == 0);\n"
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
"\t\t(* anyconst *) reg	[LGSIZE:0]\tf_const_addr;\n"
"\t\twire\t[LGSIZE:0]\tf_reversed_addr;\n"
"\t\treg\t		f_addr_loaded;\n"
"\t\treg\t[(2*WIDTH-1):0]\tf_addr_value;\n"
"\n"
"\t\tgenerate for(k=0; k<LGSIZE; k=k+1)\n"
"\t\t\tassign\tf_reversed_addr[k] = f_const_addr[LGSIZE-1-k];\n"
"\t\tendgenerate\n"
"\t\tassign\tf_reversed_addr[LGSIZE] = f_const_addr[LGSIZE];\n"
"\n"
"\t\tinitial\tf_addr_loaded = 1'b0;\n"
"\t\talways @(posedge i_clk)\n"
"\t\tif (i_reset)\n"
"\t\t\tf_addr_loaded <= 1'b0;\n"
"\t\telse if (i_ce)\n"
"\t\tbegin\n"
"\t\t\tif (wraddr == f_const_addr)\n"
"\t\t\t\tf_addr_loaded <= 1'b1;\n"
"\t\t\telse if (rdaddr == f_const_addr)\n"
"\t\t\t\tf_addr_loaded <= 1'b0;\n"
"\t\tend\n"
"\n"
"\t\talways @(posedge i_clk)\n"
"\t\tif ((i_ce)&&(wraddr == f_const_addr))\n"
"\t\tbegin\n"
"\t\t\tf_addr_value <= i_in;\n"
"\t\t\t`ASSERT(!f_addr_loaded);\n"
"\t\tend\n"
"\n"
"\t\talways @(posedge i_clk)\n"
"\t\tif ((f_past_valid)&&(!$past(i_reset))\n"
"\t\t\t\t&&($past(f_addr_loaded))&&(!f_addr_loaded))\n"
"\t\t\tassert(o_out == f_addr_value);\n"
"\n"
		"\t\talways @(*)\n"
		"\t\tif (o_sync)\n"
			"\t\t\tassert(wraddr[LGSIZE-1:0] == 1);\n"
"\n"
"\t\talways @(*)\n"
"\t\tif ((wraddr[LGSIZE]==f_const_addr[LGSIZE])\n"
"\t\t\t\t&&(wraddr[LGSIZE-1:0]\n"
"\t\t\t\t\t\t<= f_const_addr[LGSIZE-1:0]))\n"
"\t\t\t`ASSERT(!f_addr_loaded);\n"
"\n"
"\t\talways @(*)\n"
"\t\tif ((rdaddr[LGSIZE]==f_const_addr[LGSIZE])&&(f_addr_loaded))\n"
"\t\t\t`ASSERT(wraddr[LGSIZE-1:0]\n"
"\t\t\t\t\t<= f_reversed_addr[LGSIZE-1:0]+1);\n"
"\n"
"\t\talways @(*)\n"
"\t\tif (f_addr_loaded)\n"
"\t\t\t`ASSERT(brmem[f_const_addr] == f_addr_value);\n"
"\n"
"\n\n");

		fprintf(fp,
"`endif\t// FORMAL\n");
	}

	fprintf(fp,
"endmodule\n");

	fclose(fp);
	free(modulename);
}

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
"//\n"
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
"\n\n"
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
"//\n");
	fprintf(fp,
"module	%s(i_clk, %s, i_ce, i_in_0, i_in_1,\n"
	"\t\to_out_0, o_out_1, o_sync);\n"
	"\tparameter\t\t\tLGSIZE=%d, WIDTH=24;\n"
	"\tinput\twire\t\t\ti_clk, %s, i_ce;\n"
	"\tinput\twire\t[(2*WIDTH-1):0]\ti_in_0, i_in_1;\n"
	"\toutput\twire\t[(2*WIDTH-1):0]\to_out_0, o_out_1;\n"
	"\toutput\treg\t\t\to_sync;\n", modulename,
		resetw.c_str(), TST_DBLREVERSE_LGSIZE, resetw.c_str());

	fprintf(fp,
"\n"
	"\treg\t\t\tin_reset;\n"
	"\treg\t[(LGSIZE-1):0]\tiaddr;\n"
	"\twire\t[(LGSIZE-3):0]\tbraddr;\n"
"\n"
	"\tgenvar\tk;\n"
	"\tgenerate for(k=0; k<LGSIZE-2; k=k+1)\n"
	"\tbegin : gen_a_bit_reversed_value\n"
		"\t\tassign braddr[k] = iaddr[LGSIZE-3-k];\n"
	"\tend endgenerate\n"
"\n"
	"\tinitial iaddr = 0;\n"
	"\tinitial in_reset = 1\'b1;\n"
	"\tinitial o_sync = 1\'b0;\n");

	if (async_reset)
		fprintf(fp, "\talways @(posedge i_clk, negedge i_areset_n)\n\t\tif (!i_areset_n)\n");
	else
		fprintf(fp, "\talways @(posedge i_clk)\n\t\tif (i_reset)\n");
	fprintf(fp,
		"\t\tbegin\n"
			"\t\t\tiaddr <= 0;\n"
			"\t\t\tin_reset <= 1\'b1;\n"
			"\t\t\to_sync <= 1\'b0;\n"
		"\t\tend else if (i_ce)\n"
		"\t\tbegin\n"
			"\t\t\tiaddr <= iaddr + { {(LGSIZE-1){1\'b0}}, 1\'b1 };\n"
			"\t\t\tif (&iaddr[(LGSIZE-2):0])\n"
				"\t\t\t\tin_reset <= 1\'b0;\n"
			"\t\t\tif (in_reset)\n"
				"\t\t\t\to_sync <= 1\'b0;\n"
			"\t\t\telse\n"
				"\t\t\t\to_sync <= ~(|iaddr[(LGSIZE-2):0]);\n"
		"\t\tend\n"
"\n"
	"\treg\t[(2*WIDTH-1):0]\tmem_e [0:((1<<(LGSIZE))-1)];\n"
	"\treg\t[(2*WIDTH-1):0]\tmem_o [0:((1<<(LGSIZE))-1)];\n"
"\n"
	"\talways @(posedge i_clk)\n"
		"\t\tif (i_ce)\tmem_e[iaddr] <= i_in_0;\n"
	"\talways @(posedge i_clk)\n"
		"\t\tif (i_ce)\tmem_o[iaddr] <= i_in_1;\n"
"\n"
"\n"
	"\treg [(2*WIDTH-1):0] evn_out_0, evn_out_1, odd_out_0, odd_out_1;\n"
"\n"
	"\talways @(posedge i_clk)\n"
		"\t\tif (i_ce)\n\t\t\tevn_out_0 <= mem_e[{!iaddr[LGSIZE-1],1\'b0,braddr}];\n"
	"\talways @(posedge i_clk)\n"
		"\t\tif (i_ce)\n\t\t\tevn_out_1 <= mem_e[{!iaddr[LGSIZE-1],1\'b1,braddr}];\n"
	"\talways @(posedge i_clk)\n"
		"\t\tif (i_ce)\n\t\t\todd_out_0 <= mem_o[{!iaddr[LGSIZE-1],1\'b0,braddr}];\n"
	"\talways @(posedge i_clk)\n"
		"\t\tif (i_ce)\n\t\t\todd_out_1 <= mem_o[{!iaddr[LGSIZE-1],1\'b1,braddr}];\n"
"\n"
	"\treg\tadrz;\n"
	"\talways @(posedge i_clk)\n"
		"\t\tif (i_ce) adrz <= iaddr[LGSIZE-2];\n"
"\n"
	"\tassign\to_out_0 = (adrz)?odd_out_0:evn_out_0;\n"
	"\tassign\to_out_1 = (adrz)?odd_out_1:evn_out_1;\n"
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
	"\t\t(* anyconst *) reg	[LGSIZE-1:0]	f_const_addr;\n"
	"\t\twire	[LGSIZE-3:0]	f_reversed_addr;\n"
	"\t\t// reg	[LGSIZE:0]	f_now;\n"
	"\t\treg			f_addr_loaded_0, f_addr_loaded_1;\n"
	"\t\treg	[(2*WIDTH-1):0]	f_data_0, f_data_1;\n"
	"\t\twire			f_writing, f_reading;\n"
"\n"
	"\t\tgenerate for(k=0; k<LGSIZE-2; k=k+1)\n"
	"\t\t	assign	f_reversed_addr[k] = f_const_addr[LGSIZE-3-k];\n"
	"\t\tendgenerate\n"
"\n"
	"\t\tassign	f_writing=(f_const_addr[LGSIZE-1]==iaddr[LGSIZE-1]);\n"
	"\t\tassign	f_reading=(f_const_addr[LGSIZE-1]!=iaddr[LGSIZE-1]);\n"
	"\t\tinitial	f_addr_loaded_0 = 1'b0;\n"
	"\t\tinitial	f_addr_loaded_1 = 1'b0;\n"
	"\t\talways @(posedge i_clk)\n"
	"\t\tif (i_reset)\n"
	"\t\tbegin\n"
	"\t		f_addr_loaded_0 <= 1'b0;\n"
	"\t		f_addr_loaded_1 <= 1'b0;\n"
	"\t\tend else if (i_ce)\n"
	"\t\tbegin\n"
	"\t\t	if (iaddr == f_const_addr)\n"
	"\t\t	begin\n"
	"\t\t		f_addr_loaded_0 <= 1'b1;\n"
	"\t\t		f_addr_loaded_1 <= 1'b1;\n"
	"\t\t	end\n"
"\n"
	"\t\t	if (f_reading)\n"
	"\t\t	begin\n"
	"\t\t		if ((braddr == f_const_addr[LGSIZE-3:0])\n"
	"\t\t			&&(iaddr[LGSIZE-2] == 1'b0))\n"
	"\t\t			f_addr_loaded_0 <= 1'b0;\n"
"\n"
	"\t\t		if ((braddr == f_const_addr[LGSIZE-3:0])\n"
	"\t\t			&&(iaddr[LGSIZE-2] == 1'b1))\n"
	"\t\t			f_addr_loaded_1 <= 1'b0;\n"
	"\t\t	end\n"
	"\t\tend\n"
"\n"
	"\t\talways @(posedge i_clk)\n"
	"\t\tif ((i_ce)&&(iaddr == f_const_addr))\n"
	"\t\tbegin\n"
	"\t\t	f_data_0 <= i_in_0;\n"
	"\t\t	f_data_1 <= i_in_1;\n"
	"\t\t	`ASSERT(!f_addr_loaded_0);\n"
	"\t\t	`ASSERT(!f_addr_loaded_1);\n"
	"\t\tend\n"
"\n"
	"\t\talways @(posedge i_clk)\n"
	"\t\tif ((f_past_valid)&&(!$past(i_reset))\n"
	"\t\t		&&($past(f_addr_loaded_0))&&(!f_addr_loaded_0))\n"
	"\t\tbegin\n"
	"\t\t	assert(!$past(iaddr[LGSIZE-2]));\n"
	"\t\t	if (f_const_addr[LGSIZE-2])\n"
	"\t\t		assert(o_out_1 == f_data_0);\n"
	"\t\t	else\n"
	"\t\t		assert(o_out_0 == f_data_0);\n"
	"\t\tend\n"
"\n"
	"\t\talways @(posedge i_clk)\n"
	"\t\tif ((f_past_valid)&&(!$past(i_reset))\n"
	"\t\t		&&($past(f_addr_loaded_1))&&(!f_addr_loaded_1))\n"
	"\t\tbegin\n"
	"\t\t	assert($past(iaddr[LGSIZE-2]));\n"
	"\t\t	if (f_const_addr[LGSIZE-2])\n"
	"\t\t		assert(o_out_1 == f_data_1);\n"
	"\t\t	else\n"
	"\t\t		assert(o_out_0 == f_data_1);\n"
	"\t\tend\n"
"\n"
	"\t\talways @(*)\n"
	"\t\t	`ASSERT(o_sync == ((iaddr[LGSIZE-2:0] == 1)&&(!in_reset)));\n"
"\n"
	"\t\t// Before writing to a section, the loaded flags should be\n"
	"\t\t// zero\n"
	"\t\talways @(*)\n"
	"\t\tif (f_writing)\n"
	"\t\tbegin\n"
	"\t\t	`ASSERT(f_addr_loaded_0 == (iaddr[LGSIZE-2:0]\n"
	"\t\t				> f_const_addr[LGSIZE-2:0]));\n"
	"\t\t	`ASSERT(f_addr_loaded_1 == (iaddr[LGSIZE-2:0]\n"
	"\t\t				> f_const_addr[LGSIZE-2:0]));\n"
	"\t\tend\n"
"\n"
	"\t\t// If we were writing, and now we are reading, then both\n"
	"\t\t// f_addr_loaded flags must be set\n"
	"\t\talways @(posedge i_clk)\n"
	"\t\tif ((f_past_valid)&&(!$past(i_reset))\n"
	"\t\t		&&($past(f_writing))&&(f_reading))\n"
	"\t\tbegin\n"
	"\t\t	`ASSERT(f_addr_loaded_0);\n"
	"\t\t	`ASSERT(f_addr_loaded_1);\n"
	"\t\tend\n"
"\n"
	"\t\talways @(*)\n"
	"\t\tif (f_writing)\n"
	"\t\t	`ASSERT(f_addr_loaded_0 == f_addr_loaded_1);\n"
"\n"
	"\t\t// When reading, and the loaded flag is zero, our pointer\n"
	"\t\t// must not have hit the address of interest yet\n"
	"\t\talways @(*)\n"
	"\t\tif ((!in_reset)&&(f_reading))\n"
	"\t\t	`ASSERT(f_addr_loaded_0 ==\n"
	"\t\t		((!iaddr[LGSIZE-2])&&(iaddr[LGSIZE-3:0]\n"
	"\t\t			<= f_reversed_addr[LGSIZE-3:0])));\n"
	"\t\talways @(*)\n"
	"\t\tif ((!in_reset)&&(f_reading))\n"
	"\t\t	`ASSERT(f_addr_loaded_1 ==\n"
	"\t\t		((!iaddr[LGSIZE-2])||(iaddr[LGSIZE-3:0]\n"
	"\t\t			<= f_reversed_addr[LGSIZE-3:0])));\n"
	"\t\talways @(*)\n"
	"\t\tif ((in_reset)&&(f_reading))\n"
	"\t\tbegin\n"
	"\t\t	`ASSERT(!f_addr_loaded_0);\n"
	"\t\t	`ASSERT(!f_addr_loaded_1);\n"
	"\t\tend\n"
"\n"
	"\t\talways @(*)\n"
	"\t\tif(iaddr[LGSIZE-1])\n"
	"\t\t	`ASSERT(!in_reset);\n"
"\n"
	"\t\talways @(*)\n"
	"\t\tif (f_addr_loaded_0)\n"
	"\t\t	`ASSERT(mem_e[f_const_addr] == f_data_0);\n"
	"\t\talways @(*)\n"
	"\t\tif (f_addr_loaded_1)\n"
	"\t\t	`ASSERT(mem_o[f_const_addr] == f_data_1);\n"
"\n\n");


		fprintf(fp,
"`endif\t// FORMAL\n");
	}

	fprintf(fp,
"endmodule\n");

	fclose(fp);
	free(modulename);
}
