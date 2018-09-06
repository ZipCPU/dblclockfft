////////////////////////////////////////////////////////////////////////////////
//
// Filename: 	../rtl/abs_longbimpy.v
//
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
`default_nettype	none
//
module	longbimpy(i_clk, i_ce, i_a_unsorted, i_b_unsorted, o_r);
	parameter	IAW=8,	// The width of i_a, min width is 5
			IBW=12,	// The width of i_b, can be anything
			// The following three parameters should not be changed
			// by any implementation, but are based upon hardware
			// and the above values:
			OW=IAW+IBW;	// The output width
	localparam	AW = (IAW<IBW) ? IAW : IBW,
			BW = (IAW<IBW) ? IBW : IAW,
			IW=(AW+1)&(-2),	// Internal width of A
			LUTB=2,	// How many bits we can multiply by at once
			TLEN=(AW+(LUTB-1))/LUTB; // Nmbr of rows in our tableau
	input	wire			i_clk, i_ce;
	input	wire	[(IAW-1):0]	i_a_unsorted;
	input	wire	[(IBW-1):0]	i_b_unsorted;
	output	reg	[(AW+BW-1):0]	o_r;

	//
	// Swap parameter order, so that AW <= BW -- for performance
	// reasons
	wire	[AW-1:0]	i_a;
	wire	[BW-1:0]	i_b;
	generate if (IAW <= IBW)
	begin : NO_PARAM_CHANGE
		assign i_a = i_a_unsorted;
		assign i_b = i_b_unsorted;
	end else begin : SWAP_PARAMETERS
		assign i_a = i_b_unsorted;
		assign i_b = i_a_unsorted;
	end endgenerate

`ifndef	FORMAL
	// This file should only be used in a formal context.
	// The following line should therefore yield a syntax error
	assert(0);
`endif

	reg	[(IW-1):0]	u_a;
	reg	[(BW-1):0]	u_b;
	reg			sgn;

	genvar k;

	// First step:
	// Switch to unsigned arithmetic for our multiply, keeping track
	// of the along the way.  We'll then add the sign again later at
	// the end.
	//
	// If we were forced to stay within two's complement arithmetic,
	// taking the absolute value here would require an additional bit.
	// However, because our results are now unsigned, we can stay
	// within the number of bits given (for now).
	initial u_a = 0;
	generate if (IW > AW)
	begin
		always @(posedge i_clk)
			if (i_ce)
				u_a <= { 1'b0, (i_a[AW-1])?(-i_a):(i_a) };
	end else begin
		always @(posedge i_clk)
			if (i_ce)
				u_a <= (i_a[AW-1])?(-i_a):(i_a);
	end endgenerate

	initial sgn = 0;
	initial u_b = 0;
	always @(posedge i_clk)
	if (i_ce)
	begin
		u_b <= (i_b[BW-1])?(-i_b):(i_b);
		sgn <= i_a[AW-1] ^ i_b[BW-1];
	end

`ifdef	FORMAL
	reg	f_past_valid;
	initial	f_past_valid = 1'b0;
	always @(posedge i_clk)
		f_past_valid <= 1'b1;

`define	ASSERT	assume

	reg	[AW-1:0]	f_past_a	[0:TLEN];
	reg	[BW-1:0]	f_past_b	[0:TLEN];
	reg	[TLEN+1:0]	f_sgn_a, f_sgn_b;

	initial	f_past_a[0] = 0;
	initial	f_past_b[0] = 0;
	initial	f_sgn_a = 0;
	initial	f_sgn_b = 0;
	always @(posedge i_clk)
	if (i_ce)
	begin
		f_past_a[0] <= u_a;
		f_past_b[0] <= u_b;
		f_sgn_a[0] <= i_a[AW-1];
		f_sgn_b[0] <= i_b[BW-1];
	end

	generate for(k=0; k<TLEN; k=k+1)
	begin
		initial	f_past_a[k+1] = 0;
		initial	f_past_b[k+1] = 0;
		initial	f_sgn_a[k+1] = 0;
		initial	f_sgn_b[k+1] = 0;
		always @(posedge i_clk)
		if (i_ce)
		begin
			f_past_a[k+1] <= f_past_a[k];
			f_past_b[k+1] <= f_past_b[k];

			f_sgn_a[k+1]  <= f_sgn_a[k];
			f_sgn_b[k+1]  <= f_sgn_b[k];
		end
	end endgenerate

	always @(posedge i_clk)
	if (i_ce)
	begin
		f_sgn_a[TLEN+1] <= f_sgn_a[TLEN];
		f_sgn_b[TLEN+1] <= f_sgn_b[TLEN];
	end

	always @(posedge i_clk)
		assert(sgn == (f_sgn_a[0] ^ f_sgn_b[0]));

	wire	[BW-1:0]	f_past_b_neg = - f_past_b[TLEN];

	always @(posedge i_clk)
	if (f_past_a[TLEN]==0)
		`ASSERT(o_r == 0);
	else if ((f_past_valid)&&($past(i_ce))
		&&(f_past_a[TLEN]==1)&&(!f_sgn_a[TLEN+1]))
	begin
		if (!f_sgn_b[TLEN+1])
		begin
			`ASSERT(o_r[BW-1:0] == f_past_b[TLEN][BW-1:0]);
			`ASSERT(o_r[AW+BW-1:BW] == 0);
		end else begin
			`ASSERT(o_r[BW-1:0] == f_past_b_neg);
			`ASSERT(&o_r[AW+BW-1:BW]);
		end
	end

	wire	[AW-1:0]	f_past_a_neg = - f_past_a[TLEN];

	always @(posedge i_clk)
	if (f_past_b[TLEN]==0)
		`ASSERT(o_r == 0);
	else if ((f_past_valid)&&($past(i_ce))
		&&(f_past_b[TLEN]==1)&&(!f_sgn_b[TLEN+1]))
	begin
		if (!f_sgn_a[TLEN+1])
		begin
			`ASSERT(o_r[AW-1:0] == f_past_a[TLEN][AW-1:0]);
			`ASSERT(o_r[AW+BW-1:AW] == 0);
		end else begin
			`ASSERT(o_r[AW-1:0] == f_past_a_neg);
			`ASSERT(&o_r[AW+BW-1:AW]);
		end
	end
`endif
endmodule
