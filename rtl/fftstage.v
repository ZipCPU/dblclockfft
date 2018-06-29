////////////////////////////////////////////////////////////////////////////////
//
// Filename:	fftstage.v
//
// Project:	A General Purpose Pipelined FFT Implementation
//
// Purpose:	This file is (almost) a Verilog source file.  It is meant to
//		be used by a FFT core compiler to generate FFTs which may be
//	used as part of an FFT core.  Specifically, this file encapsulates
//	the options of an FFT-stage.  For any 2^N length FFT, there shall be
//	(N-1) of these stages.
//
//
// Operation:
// 	Given a stream of values, operate upon them as though they were
// 	value pairs, x[n] and x[n+N/2].  The stream begins when n=0, and ends
// 	when n=N/2-1 (i.e. there's a full set of N values).  When the value
// 	x[0] enters, the synchronization input, i_sync, must be true as well.
//
// 	For this stream, produce outputs
// 	y[n    ] = x[n] + x[n+N/2], and
// 	y[n+N/2] = (x[n] - x[n+N/2]) * c[n],
// 			where c[n] is a complex coefficient found in the
// 			external memory file COEFFILE.
// 	When y[0] is output, a synchronization bit o_sync will be true as
// 	well, otherwise it will be zero.
//
// 	Most of the work to do this is done within the butterfly, whether the
// 	hardware accelerated butterfly (uses a DSP) or not.
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
module	fftstage(i_clk, i_reset, i_ce, i_sync, i_data, o_data, o_sync);
	parameter	IWIDTH=16,CWIDTH=20,OWIDTH=17;
	// Parameters specific to the core that should be changed when this
	// core is built ... Note that the minimum LGSPAN (the base two log
	// of the span, or the base two log of the current FFT size) is 3.
	// Smaller spans (i.e. the span of 2) must use the dbl laststage module.
	parameter	LGWIDTH=11, LGSPAN=10, LGBDLY=4, BFLYSHIFT=0;
	parameter	[0:0]	OPT_HWMPY = 1'b1;
	// Clocks per CE.  If your incoming data rate is less than 50% of your
	// clock speed, you can set CKPCE to 2'b10, make sure there's at least
	// one clock between cycles when i_ce is high, and then use two
	// multiplies instead of three.  Setting CKPCE to 2'b11, and insisting
	// on at least two clocks with i_ce low between cycles with i_ce high,
	// then the hardware optimized butterfly code will used one multiply
	// instead of two.
	parameter	[1:0]	CKPCE = 2'h2;
	// The COEFFILE parameter contains the name of the file containing the
	// FFT twiddle factors
	parameter	COEFFILE="cmem_2048.hex";
	input					i_clk, i_reset, i_ce, i_sync;
	input		[(2*IWIDTH-1):0]	i_data;
	output	reg	[(2*OWIDTH-1):0]	o_data;
	output	reg				o_sync;

	reg	wait_for_sync;
	reg	[(2*IWIDTH-1):0]	ib_a, ib_b;
	reg	[(2*CWIDTH-1):0]	ib_c;
	reg	ib_sync;

	reg	b_started;
	wire	ob_sync;
	wire	[(2*OWIDTH-1):0]	ob_a, ob_b;

	// cmem is defined as an array of real and complex values,
	// where the top CWIDTH bits are the real value and the bottom
	// CWIDTH bits are the imaginary value.
	//
	// cmem[i] = { (2^(CWIDTH-2)) * cos(2*pi*i/(2^LGWIDTH)),
	//		(2^(CWIDTH-2)) * sin(2*pi*i/(2^LGWIDTH)) };
	//
	reg	[(2*CWIDTH-1):0]	cmem [0:((1<<LGSPAN)-1)];
	initial	$readmemh(COEFFILE,cmem);

	reg	[(LGSPAN):0]		iaddr;
	reg	[(2*IWIDTH-1):0]	imem	[0:((1<<LGSPAN)-1)];

	reg	[LGSPAN:0]		oB;
	reg	[(2*OWIDTH-1):0]	omem	[0:((1<<LGSPAN)-1)];

	initial wait_for_sync = 1'b1;
	initial iaddr = 0;
	always @(posedge i_clk)
		if (i_reset)
		begin
			wait_for_sync <= 1'b1;
			iaddr <= 0;
		end
		else if ((i_ce)&&((!wait_for_sync)||(i_sync)))
		begin
			//
			// First step: Record what we're not ready to use yet
			//
			iaddr <= iaddr + { {(LGSPAN){1'b0}}, 1'b1 };
			wait_for_sync <= 1'b0;
		end
	always @(posedge i_clk) // Need to make certain here that we don't read
		if ((i_ce)&&(!iaddr[LGSPAN])) // and write the same address on
			imem[iaddr[(LGSPAN-1):0]] <= i_data; // the same clk

	//
	// Now, we have all the inputs, so let's feed the butterfly
	//
	initial ib_sync = 1'b0;
	always @(posedge i_clk)
		if (i_reset)
			ib_sync <= 1'b0;
		else if ((i_ce)&&(iaddr[LGSPAN]))
			begin
				// Set the sync to true on the very first
				// valid input in, and hence on the very
				// first valid data out per FFT.
				ib_sync <= (iaddr==(1<<(LGSPAN)));
			end
	always	@(posedge i_clk)
		if ((i_ce)&&(iaddr[LGSPAN]))
			begin
				// One input from memory, ...
				ib_a <= imem[iaddr[(LGSPAN-1):0]];
				// One input clocked in from the top
				ib_b <= i_data;
				// and the coefficient or twiddle factor
				ib_c <= cmem[iaddr[(LGSPAN-1):0]];
			end

	generate if (OPT_HWMPY)
	begin : HWBFLY
		hwbfly #(.IWIDTH(IWIDTH),.CWIDTH(CWIDTH),.OWIDTH(OWIDTH),
				.CKPCE(CKPCE), .SHIFT(BFLYSHIFT))
			bfly(i_clk, i_reset, i_ce, ib_c,
				ib_a, ib_b, ib_sync, ob_a, ob_b, ob_sync);
	end else begin : FWBFLY
		butterfly #(.IWIDTH(IWIDTH),.CWIDTH(CWIDTH),.OWIDTH(OWIDTH),
				.MPYDELAY(4'd11),.LGDELAY(LGBDLY),
				.CKPCE(CKPCE),.SHIFT(BFLYSHIFT))
			bfly(i_clk, i_reset, i_ce, ib_c,
				ib_a, ib_b, ib_sync, ob_a, ob_b, ob_sync);
	end endgenerate

	//
	// Next step: recover the outputs from the butterfly
	//
	initial oB        = 0;
	initial o_sync    = 0;
	initial b_started = 0;
	always @(posedge i_clk)
		if (i_reset)
		begin
			oB <= 0;
			o_sync <= 0;
			b_started <= 0;
		end else if (i_ce)
		begin
			o_sync <= (!oB[LGSPAN])?ob_sync : 1'b0;
			if (ob_sync||b_started)
				oB <= oB + { {(LGSPAN){1'b0}}, 1'b1 };
			if ((ob_sync)&&(!oB[LGSPAN]))
			// A butterfly output is available
				b_started <= 1'b1;
		end

	reg	[(LGSPAN-1):0]		dly_addr;
	reg	[(2*OWIDTH-1):0]	dly_value;
	always @(posedge i_clk)
		if (i_ce)
		begin
			dly_addr <= oB[(LGSPAN-1):0];
			dly_value <= ob_b;
		end
	always @(posedge i_clk)
		if (i_ce)
			omem[dly_addr] <= dly_value;

	always @(posedge i_clk)
		if (i_ce)
			o_data <= (!oB[LGSPAN])?ob_a : omem[oB[(LGSPAN-1):0]];

endmodule
