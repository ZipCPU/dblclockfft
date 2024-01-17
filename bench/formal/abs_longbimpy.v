////////////////////////////////////////////////////////////////////////////////
//
// Filename: 	bench/formal/abs_longbimpy.v
// {{{
// Project:	A General Purpose Pipelined FFT Implementation
//
// Purpose:	A portable shift and add multiply, built with the knowledge
//	of the existence of a six bit LUT and carry chain.  That knowledge
//	allows us to multiply two bits from one value at a time against all
//	of the bits of the other value.  This sub multiply is called the
//	bimpy.
//
//	For minimal processing delay, make the first parameter the one with
//	the least bits, so that AWIDTH <= BWIDTH.
//
//
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
// with this program.  (It's in the $(ROOT)/doc directory, run make with no
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
//
`default_nettype	none
// }}}
module	longbimpy #(
		// {{{
		parameter	IAW=8,	// The width of i_a, min width is 5
				IBW=12,	// The width of i_b, can be anything
				// The following three parameters should not be
				// changed by any implementation, but are based
				// upon hardware and the above values:
				// OW=IAW+IBW;	// The output width
		localparam	AW = (IAW<IBW) ? IAW : IBW,
				BW = (IAW<IBW) ? IBW : IAW,
				LUTB=2,	// How many bits we can mpy by at once
				TLEN=(AW+(LUTB-1))/LUTB // Nbr rows in tableau
		// }}}
	) (
		// {{{
		input	wire			i_clk, i_ce,
		input	wire	[(IAW-1):0]	i_a_unsorted,
		input	wire	[(IBW-1):0]	i_b_unsorted,
		output	reg	[(AW+BW-1):0]	o_r
		// }}}
	);

	// Signal declarations
	// {{{
	reg	f_past_valid;
	wire	[AW-1:0]	i_a;
	wire	[BW-1:0]	i_b;
	reg	[AW-1:0]	f_past_a	[0:TLEN+1];
	reg	[BW-1:0]	f_past_b	[0:TLEN+1];

	// }}}

	// Swap parameter order, so that AW <= BW -- for performance reasons
	// {{{
	generate if (IAW <= IBW)
	begin : NO_PARAM_CHANGE
		assign i_a = i_a_unsorted;
		assign i_b = i_b_unsorted;
	end else begin : SWAP_PARAMETERS
		assign i_a = i_b_unsorted;
		assign i_b = i_a_unsorted;
	end endgenerate
	// }}}

`ifndef	FORMAL
	// This file should only be used in a formal context.
	// The following line should therefore yield a syntax error
	assert(0);
`endif

	initial	f_past_valid = 1'b0;
	always @(posedge i_clk)
		f_past_valid <= 1'b1;

	// f_past_a, f_past_b
	// {{{
	initial	f_past_a[0] = 0;
	initial	f_past_b[0] = 0;
	always @(posedge i_clk)
	if (i_ce)
	begin
		f_past_a[0] <= i_a;
		f_past_b[0] <= i_b;
	end

	genvar	k;

	generate for(k=0; k<TLEN+1; k=k+1)
	begin
		initial	f_past_a[k+1] = 0;
		initial	f_past_b[k+1] = 0;
		always @(posedge i_clk)
		if (i_ce)
		begin
			f_past_a[k+1] <= f_past_a[k];
			f_past_b[k+1] <= f_past_b[k];
		end
	end endgenerate
	// }}}

	// abs_mpy #(.AW(AW), .BW(BW)) thempy(f_past_a[TLEN+1], f_past_b[TLEN+1], o_r);
	(* anyseq *) reg [AW+BW-1:0]	result;
	wire	[AW+BW-1:0]	f_neg_a, f_neg_b;

	// Negate f_past_a and f_past_b
	// {{{
	assign	f_neg_a = - {{(BW){f_past_a[TLEN+1][AW-1]}}, f_past_a[TLEN+1]};
	assign	f_neg_b = - {{(AW){f_past_b[TLEN+1][BW-1]}}, f_past_b[TLEN+1]};
	// }}}

	// Assumptions about the result
	// {{{
	always @(*)
	if (f_past_a[TLEN+1] == 0)
		assume(result == 0);
	else if (f_past_b[TLEN+1] == 0)
		assume(result == 0);
	else if (f_past_a[TLEN+1] == 1)
	begin
		assume(result[BW-1:0] == f_past_b[TLEN+1]);
		assume(result[AW+BW-1:BW] == {(AW){f_past_b[TLEN+1][BW-1]}});
	end else if (f_past_b[TLEN+1] == 1)
	begin
		assume(result[AW-1:0] == f_past_a[TLEN+1]);
		assume(result[AW+BW-1:AW] == {(BW){f_past_a[TLEN+1][AW-1]}});
	end else if (&f_past_a[TLEN+1])
		assume(result == f_neg_b);
	else if (&f_past_b[TLEN+1])
		assume(result == f_neg_a);
	else
		assume(result[AW+BW-1] == (f_past_a[TLEN+1][AW-1]
					^f_past_b[TLEN+1][BW-1]));
	// }}}

	always @(*)
		o_r = result;
endmodule
