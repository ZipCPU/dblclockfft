////////////////////////////////////////////////////////////////////////////////
//
// Filename:	shiftaddmpy.v
//
// Project:	A General Purpose Pipelined FFT Implementation
//
// Purpose:	A portable shift and add multiply.
//
//	While both Xilinx and Altera will offer single clock multiplies, this
//	simple approach will multiply two numbers on any architecture.  The
//	result maintains the full width of the multiply, there are no extra
//	stuff bits, no rounding, no shifted bits, etc.
//
//	Further, for those applications that can support it, this multiply
//	is pipelined and will produce one answer per clock.
//
//	For minimal processing delay, make the first parameter the one with
//	the least bits, so that AWIDTH <= BWIDTH.
//
//	The processing delay in this multiply is (AWIDTH+1) cycles.  That is,
//	if the data is present on the input at clock t=0, the result will be
//	present on the output at time t=AWIDTH+1;
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
module	shiftaddmpy(i_clk, i_ce, i_a, i_b, o_r);
	parameter	AWIDTH=16,BWIDTH=20;
	input					i_clk, i_ce;
	input		[(AWIDTH-1):0]		i_a;
	input		[(BWIDTH-1):0]		i_b;
	output	reg	[(AWIDTH+BWIDTH-1):0]	o_r;

	reg	[(AWIDTH-1):0]	u_a;
	reg	[(BWIDTH-1):0]	u_b;
	reg			sgn;

	reg	[(AWIDTH-2):0]		r_a[0:(AWIDTH-1)];
	reg	[(AWIDTH+BWIDTH-2):0]	r_b[0:(AWIDTH-1)];
	reg				r_s[0:(AWIDTH-1)];
	reg	[(AWIDTH+BWIDTH-1):0]	acc[0:(AWIDTH-1)];
	genvar k;

	// If we were forced to stay within two's complement arithmetic,
	// taking the absolute value here would require an additional bit.
	// However, because our results are now unsigned, we can stay
	// within the number of bits given (for now).
	always @(posedge i_clk)
		if (i_ce)
		begin
			u_a <= (i_a[AWIDTH-1])?(-i_a):(i_a);
			u_b <= (i_b[BWIDTH-1])?(-i_b):(i_b);
			sgn <= i_a[AWIDTH-1] ^ i_b[BWIDTH-1];
		end

	always @(posedge i_clk)
		if (i_ce)
		begin
			acc[0] <= (u_a[0]) ? { {(AWIDTH){1'b0}}, u_b }
					: {(AWIDTH+BWIDTH){1'b0}};
			r_a[0] <= { u_a[(AWIDTH-1):1] };
			r_b[0] <= { {(AWIDTH-1){1'b0}}, u_b };
			r_s[0] <= sgn; // The final sign, needs to be preserved
		end

	generate
	for(k=0; k<AWIDTH-1; k=k+1)
	begin : genstages
		always @(posedge i_clk)
		if (i_ce)
		begin
			acc[k+1] <= acc[k] + ((r_a[k][0]) ? {r_b[k],1'b0}:0);
			r_a[k+1] <= { 1'b0, r_a[k][(AWIDTH-2):1] };
			r_b[k+1] <= { r_b[k][(AWIDTH+BWIDTH-3):0], 1'b0};
			r_s[k+1] <= r_s[k];
		end
	end
	endgenerate

	always @(posedge i_clk)
		if (i_ce)
			o_r <= (r_s[AWIDTH-1]) ? (-acc[AWIDTH-1]) : acc[AWIDTH-1];

endmodule
