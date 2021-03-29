////////////////////////////////////////////////////////////////////////////////
//
// Filename:	laststage.v
// {{{
// Project:	A General Purpose Pipelined FFT Implementation
//
// Purpose:	This is part of an FPGA implementation that will process
//		the final stage of a decimate-in-frequency FFT, running
//	through the data at two samples per clock.  If you notice from the
//	derivation of an FFT, the only time both even and odd samples are
//	used at the same time is in this stage.  Therefore, other than this
//	stage and these twiddles, all of the other stages can run two stages
//	at a time at one sample per clock.
//
// Operation:
// 	Given a stream of values, operate upon them as though they were
// 	value pairs, x[2n] and x[2n+1].  The stream begins when n=0, and ends
// 	when n=1.  When the first x[0] value enters, the synchronization
//	input, i_sync, must be true as well.
//
// 	For this stream, produce outputs
// 	y[2n  ] = x[2n] + x[2n+1], and
// 	y[2n+1] = x[2n] - x[2n+1]
//
// 	When y[0] is output, a synchronization bit o_sync will be true as
// 	well, otherwise it will be zero.
//
//
//	In this implementation, the output is valid one clock after the input
//	is valid.  The output also accumulates one bit above and beyond the
//	number of bits in the input.
//
//		i_clk	A system clock
//		i_reset	A synchronous reset
//		i_ce	Circuit enable--nothing happens unless this line is high
//		i_sync	A synchronization signal, high once per FFT at the start
//		i_left	The first (even) complex sample input.  The higher order
//			bits contain the real portion, low order bits the
//			imaginary portion, all in two's complement.
//		i_right	The next (odd) complex sample input, same format as
//			i_left.
//		o_left	The first (even) complex output.
//		o_right	The next (odd) complex output.
//		o_sync	Output synchronization signal.
//
// Creator:	Dan Gisselquist, Ph.D.
//		Gisselquist Technology, LLC
//
////////////////////////////////////////////////////////////////////////////////
// }}}
// Copyright (C) 2015-2021, Gisselquist Technology, LLC
// {{{
// This file is part of the general purpose pipelined FFT project.
//
// The pipelined FFT project is free software (firmware): you can redistribute
// it and/or modify it under the terms of the GNU Lesser General Public License
// as published by the Free Software Foundation, either version 3 of the
// License, or (at your option) any later version.
//
// The pipelined FFT project is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTIBILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
// General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program.  (It's in the $(ROOT)/doc directory.  Run make
// with no target there if the PDF file isn't present.)  If not, see
// <http://www.gnu.org/licenses/> for a copy.
// }}}
// License:	LGPL, v3, as defined and found on www.gnu.org,
// {{{
//		http://www.gnu.org/licenses/lgpl.html
//
// }}}
////////////////////////////////////////////////////////////////////////////////
//
//
`default_nettype	none
//
module	laststage #(
		// {{{
		parameter	IWIDTH=16,OWIDTH=IWIDTH+1, SHIFT=0
		// }}}
	) (
		// {{{
		input	wire	i_clk, i_reset, i_ce, i_sync,
		input	wire	[(2*IWIDTH-1):0]	i_left, i_right,
		output	reg	[(2*OWIDTH-1):0]	o_left, o_right,
		output	reg			o_sync

		// }}}
	);

	// Local declarations
	// {{{
	wire	signed	[(IWIDTH-1):0]	i_in_0r, i_in_0i, i_in_1r, i_in_1i;
	wire	[(OWIDTH-1):0]		o_out_0r, o_out_0i,
					o_out_1r, o_out_1i;
	reg	rnd_sync, r_sync;
	reg	signed	[(IWIDTH):0]	rnd_in_0r, rnd_in_0i;
	reg	signed	[(IWIDTH):0]	rnd_in_1r, rnd_in_1i;

	// }}}

	assign	i_in_0r = i_left[(2*IWIDTH-1):(IWIDTH)];
	assign	i_in_0i = i_left[(IWIDTH-1):0];
	assign	i_in_1r = i_right[(2*IWIDTH-1):(IWIDTH)];
	assign	i_in_1i = i_right[(IWIDTH-1):0];

	// rnd_sync, r_sync
	// {{{
	// As with any register connected to the sync pulse, these must
	// have initial values and be reset on the i_reset signal.
	// Other data values need only restrict their updates to i_ce
	// enabled clocks, but sync's must obey resets and initial
	// conditions as well.

	initial	rnd_sync      = 1'b0; // Sync into rounding
	initial	r_sync        = 1'b0; // Sync coming out
	always @(posedge i_clk)
	if (i_reset)
		begin
			rnd_sync <= 1'b0;
			r_sync <= 1'b0;
		end else if (i_ce)
		begin
			rnd_sync <= i_sync;
			r_sync <= rnd_sync;
		end
	// }}}

	// rnd_in_0r, rnd_in_0i, rnd_in_1r, rnd_in_1i
	// {{{
	// As with other variables, these are really only updated when in
	// the processing pipeline, after the first i_sync.  However, to
	// eliminate as much unnecessary logic as possible, we toggle
	// these any time the i_ce line is enabled, and don't reset.
	// them on i_reset.
	// Don't forget that we accumulate a bit by adding two values
	// together. Therefore our intermediate value must have one more
	// bit than the two originals.
	always @(posedge i_clk)
	if (i_ce)
	begin
		//
		rnd_in_0r <= i_in_0r + i_in_1r;
		rnd_in_0i <= i_in_0i + i_in_1i;
		//
		rnd_in_1r <= i_in_0r - i_in_1r;
		rnd_in_1i <= i_in_0i - i_in_1i;
		//
	end
	// }}}

	convround #(IWIDTH+1,OWIDTH,SHIFT) do_rnd_0r(i_clk, i_ce,
							rnd_in_0r, o_out_0r);

	convround #(IWIDTH+1,OWIDTH,SHIFT) do_rnd_0i(i_clk, i_ce,
							rnd_in_0i, o_out_0i);

	convround #(IWIDTH+1,OWIDTH,SHIFT) do_rnd_1r(i_clk, i_ce,
							rnd_in_1r, o_out_1r);

	convround #(IWIDTH+1,OWIDTH,SHIFT) do_rnd_1i(i_clk, i_ce,
							rnd_in_1i, o_out_1i);


	// o_left, o_right
	// {{{
	// Prior versions of this routine did not include the extra
	// clock and register/flip-flops that this routine requires.
	// These are placed in here to correct a bug in Verilator, that
	// otherwise struggles.  (Hopefully this will fix the problem ...)
	always @(posedge i_clk)
	if (i_ce)
	begin
		o_left  <= { o_out_0r, o_out_0i };
		o_right <= { o_out_1r, o_out_1i };
	end
	// }}}

	// o_sync
	// {{{
	initial	o_sync = 1'b0; // Final sync coming out of module
	always @(posedge i_clk)
	if (i_reset)
		o_sync <= 1'b0;
	else if (i_ce)
		o_sync <= r_sync;
	// }}}

endmodule
