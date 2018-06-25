////////////////////////////////////////////////////////////////////////////////
//
// Filename:	bitreverse.v
//
// Project:	A General Purpose Pipelined FFT Implementation
//
// Purpose:	This module bitreverses a pipelined FFT input.  It differes
//		from the dblreverse module in that this is just a simple and
//	straightforward bitreverse, rather than one written to handle two
//	words at once.
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
module	bitreverse(i_clk, i_reset, i_ce, i_in, o_out, o_sync);
	parameter			LGSIZE=5, WIDTH=24;
	input				i_clk, i_reset, i_ce;
	input		[(2*WIDTH-1):0]	i_in;
	output	wire	[(2*WIDTH-1):0]	o_out;
	output	reg			o_sync;
	reg	[(LGSIZE):0]	wraddr;
	wire	[(LGSIZE):0]	rdaddr;

	reg	[(2*WIDTH-1):0]	brmem	[0:((1<<(LGSIZE+1))-1)];

	genvar	k;
	generate for(k=0; k<LGSIZE; k=k+1)
		assign rdaddr[k] = wraddr[LGSIZE-1-k];
	endgenerate
	assign	rdaddr[LGSIZE] = ~wraddr[LGSIZE];

	reg	in_reset;

	initial	in_reset = 1'b1;
	always @(posedge i_clk)
		if (i_reset)
			in_reset <= 1'b1;
		else if ((i_ce)&&(&wraddr[(LGSIZE-1):0]))
			in_reset <= 1'b0;

	initial	wraddr = 0;
	always @(posedge i_clk)
		if (i_reset)
			wraddr <= 0;
		else if (i_ce)
		begin
			brmem[wraddr] <= i_in;
			wraddr <= wraddr + 1;
		end

	always @(posedge i_clk)
		if (i_ce) // If (i_reset) we just output junk ... not a problem
			o_out <= brmem[rdaddr]; // w/o a sync pulse

	initial	o_sync = 1'b0;
	always @(posedge i_clk)
		if (i_reset)
			o_sync <= 1'b0;
		else if ((i_ce)&&(!in_reset))
			o_sync <= (wraddr[(LGSIZE-1):0] == 0);

endmodule
