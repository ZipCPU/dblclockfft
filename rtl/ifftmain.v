////////////////////////////////////////////////////////////////////////////////
//
// Filename:	ifftmain.v
// {{{
// Project:	A General Purpose Pipelined FFT Implementation
//
// Purpose:	This is the main module in the General Purpose FPGA FFT
//		implementation.  As such, all other modules are subordinate
//	to this one.  This module accomplish a fixed size Complex FFT on
//	2048 data points.
//	The FFT is fully pipelined, and accepts as inputs one complex two's
//	complement sample per clock.
//
// Parameters:
//	i_clk	The clock.  All operations are synchronous with this clock.
//	i_reset	Synchronous reset, active high.  Setting this line will
//			force the reset of all of the internals to this routine.
//			Further, following a reset, the o_sync line will go
//			high the same time the first output sample is valid.
//	i_ce	A clock enable line.  If this line is set, this module
//			will accept one complex input value, and produce
//			one (possibly empty) complex output value.
//	i_sample	The complex input sample.  This value is split
//			into two two's complement numbers, 15 bits each, with
//			the real portion in the high order bits, and the
//			imaginary portion taking the bottom 15 bits.
//	o_result	The output result, of the same format as i_sample,
//			only having 21 bits for each of the real and imaginary
//			components, leading to 42 bits total.
//	o_sync	A one bit output indicating the first sample of the FFT frame.
//			It also indicates the first valid sample out of the FFT
//			on the first frame.
//
// Arguments:	This file was computer generated using the following command
//		line:
//
//		% ./fftgen -i -d ../rtl -f 2048 -1 -k 1 -p 0 -n 15 -a ../bench/cpp/ifftsize.h
//
//	This core will use hardware accelerated multiplies (DSPs)
//	for 0 of the 11 stages
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
//
//
module ifftmain(i_clk, i_reset, i_ce,
		i_sample, o_result, o_sync);
	// The bit-width of the input, IWIDTH, output, OWIDTH, and the log
	// of the FFT size.  These are localparams, rather than parameters,
	// because once the core has been generated, they can no longer be
	// changed.  (These values can be adjusted by running the core
	// generator again.)  The reason is simply that these values have
	// been hardwired into the core at several places.
	localparam	IWIDTH=15, OWIDTH=21; // LGWIDTH=11;
	//
	input	wire				i_clk, i_reset, i_ce;
	//
	input	wire	[(2*IWIDTH-1):0]	i_sample;
	output	reg	[(2*OWIDTH-1):0]	o_result;
	output	reg				o_sync;


	// Outputs of the FFT, ready for bit reversal.
	wire				br_sync;
	wire	[(2*OWIDTH-1):0]	br_result;


	wire		w_s2048;
	wire	[31:0]	w_d2048;
	fftstage	#(
		// {{{
		.IWIDTH(IWIDTH),
		.CWIDTH(IWIDTH+4),
		.OWIDTH(16),
		.LGSPAN(10),
		.BFLYSHIFT(0),
		.OPT_HWMPY(0),
		.CKPCE(1),
		.COEFFILE("icmem_2048.hex")
		// }}}
	) stage_2048(
		// {{{
		.i_clk(i_clk),
		.i_reset(i_reset),
		.i_ce(i_ce),
		.i_sync(!i_reset),
		.i_data(i_sample),
		.o_data(w_d2048),
		.o_sync(w_s2048)
		// }}}
	);


	wire		w_s1024;
	wire	[33:0]	w_d1024;
	fftstage	#(
		// {{{
		.IWIDTH(16),
		.CWIDTH(20),
		.OWIDTH(17),
		.LGSPAN(9),
		.BFLYSHIFT(0),
		.OPT_HWMPY(0),
		.CKPCE(1),
		.COEFFILE("icmem_1024.hex")
		// }}}
	) stage_1024(
		// {{{
		.i_clk(i_clk),
		.i_reset(i_reset),
		.i_ce(i_ce),
		.i_sync(w_s2048),
		.i_data(w_d2048),
		.o_data(w_d1024),
		.o_sync(w_s1024)
		// }}}
	);

	wire		w_s512;
	wire	[33:0]	w_d512;
	fftstage	#(
		// {{{
		.IWIDTH(17),
		.CWIDTH(21),
		.OWIDTH(17),
		.LGSPAN(8),
		.BFLYSHIFT(0),
		.OPT_HWMPY(0),
		.CKPCE(1),
		.COEFFILE("icmem_512.hex")
		// }}}
	) stage_512(
		// {{{
		.i_clk(i_clk),
		.i_reset(i_reset),
		.i_ce(i_ce),
		.i_sync(w_s1024),
		.i_data(w_d1024),
		.o_data(w_d512),
		.o_sync(w_s512)
		// }}}
	);

	wire		w_s256;
	wire	[35:0]	w_d256;
	fftstage	#(
		// {{{
		.IWIDTH(17),
		.CWIDTH(21),
		.OWIDTH(18),
		.LGSPAN(7),
		.BFLYSHIFT(0),
		.OPT_HWMPY(0),
		.CKPCE(1),
		.COEFFILE("icmem_256.hex")
		// }}}
	) stage_256(
		// {{{
		.i_clk(i_clk),
		.i_reset(i_reset),
		.i_ce(i_ce),
		.i_sync(w_s512),
		.i_data(w_d512),
		.o_data(w_d256),
		.o_sync(w_s256)
		// }}}
	);

	wire		w_s128;
	wire	[35:0]	w_d128;
	fftstage	#(
		// {{{
		.IWIDTH(18),
		.CWIDTH(22),
		.OWIDTH(18),
		.LGSPAN(6),
		.BFLYSHIFT(0),
		.OPT_HWMPY(0),
		.CKPCE(1),
		.COEFFILE("icmem_128.hex")
		// }}}
	) stage_128(
		// {{{
		.i_clk(i_clk),
		.i_reset(i_reset),
		.i_ce(i_ce),
		.i_sync(w_s256),
		.i_data(w_d256),
		.o_data(w_d128),
		.o_sync(w_s128)
		// }}}
	);

	wire		w_s64;
	wire	[37:0]	w_d64;
	fftstage	#(
		// {{{
		.IWIDTH(18),
		.CWIDTH(22),
		.OWIDTH(19),
		.LGSPAN(5),
		.BFLYSHIFT(0),
		.OPT_HWMPY(0),
		.CKPCE(1),
		.COEFFILE("icmem_64.hex")
		// }}}
	) stage_64(
		// {{{
		.i_clk(i_clk),
		.i_reset(i_reset),
		.i_ce(i_ce),
		.i_sync(w_s128),
		.i_data(w_d128),
		.o_data(w_d64),
		.o_sync(w_s64)
		// }}}
	);

	wire		w_s32;
	wire	[37:0]	w_d32;
	fftstage	#(
		// {{{
		.IWIDTH(19),
		.CWIDTH(23),
		.OWIDTH(19),
		.LGSPAN(4),
		.BFLYSHIFT(0),
		.OPT_HWMPY(0),
		.CKPCE(1),
		.COEFFILE("icmem_32.hex")
		// }}}
	) stage_32(
		// {{{
		.i_clk(i_clk),
		.i_reset(i_reset),
		.i_ce(i_ce),
		.i_sync(w_s64),
		.i_data(w_d64),
		.o_data(w_d32),
		.o_sync(w_s32)
		// }}}
	);

	wire		w_s16;
	wire	[39:0]	w_d16;
	fftstage	#(
		// {{{
		.IWIDTH(19),
		.CWIDTH(23),
		.OWIDTH(20),
		.LGSPAN(3),
		.BFLYSHIFT(0),
		.OPT_HWMPY(0),
		.CKPCE(1),
		.COEFFILE("icmem_16.hex")
		// }}}
	) stage_16(
		// {{{
		.i_clk(i_clk),
		.i_reset(i_reset),
		.i_ce(i_ce),
		.i_sync(w_s32),
		.i_data(w_d32),
		.o_data(w_d16),
		.o_sync(w_s16)
		// }}}
	);

	wire		w_s8;
	wire	[39:0]	w_d8;
	fftstage	#(
		// {{{
		.IWIDTH(20),
		.CWIDTH(24),
		.OWIDTH(20),
		.LGSPAN(2),
		.BFLYSHIFT(0),
		.OPT_HWMPY(0),
		.CKPCE(1),
		.COEFFILE("icmem_8.hex")
		// }}}
	) stage_8(
		// {{{
		.i_clk(i_clk),
		.i_reset(i_reset),
		.i_ce(i_ce),
		.i_sync(w_s16),
		.i_data(w_d16),
		.o_data(w_d8),
		.o_sync(w_s8)
		// }}}
	);

	wire		w_s4;
	wire	[41:0]	w_d4;
	qtrstage	#(
		// {{{
		.IWIDTH(20),
		.OWIDTH(21),
		.LGWIDTH(11),
		.INVERSE(1),
		.SHIFT(0)
		// }}}
	) stage_4(
		// {{{
		.i_clk(i_clk),
		.i_reset(i_reset),
		.i_ce(i_ce),
		.i_sync(w_s8),
		.i_data(w_d8),
		.o_data(w_d4),
		.o_sync(w_s4)
		// }}}
	);
	// verilator lint_off UNUSED
	wire		w_s2;
	// verilator lint_on  UNUSED
	wire	[41:0]	w_d2;
	laststage	#(
		// {{{
		.IWIDTH(21),
		.OWIDTH(21),
		.SHIFT(0)
		// }}}
	) stage_2(
		// {{{
		.i_clk(i_clk),
		.i_reset(i_reset),
		.i_ce(i_ce),
		.i_sync(w_s4),
		.i_val(w_d4),
		.o_val(w_d2),
		.o_sync(w_s2)
		// }}}
	);


	wire	br_start;
	reg	r_br_started;
	initial	r_br_started = 1'b0;
	always @(posedge i_clk)
	if (i_reset)
		r_br_started <= 1'b0;
	else if (i_ce)
		r_br_started <= r_br_started || w_s2;
	assign	br_start = r_br_started || w_s2;

	// Now for the bit-reversal stage.
	bitreverse	#(
		// {{{
		.LGSIZE(11), .WIDTH(21)
		// }}}
	) revstage (
		// {{{
		.i_clk(i_clk),
		.i_reset(i_reset),
		.i_ce(i_ce & br_start),
		.i_in(w_d2),
		.o_out(br_result),
		.o_sync(br_sync)
		// }}}
	);


	// Last clock: Register our outputs, we're done.
	initial	o_sync  = 1'b0;
	always @(posedge i_clk)
	if (i_reset)
		o_sync  <= 1'b0;
	else if (i_ce)
		o_sync  <= br_sync;

	always @(posedge i_clk)
	if (i_ce)
		o_result  <= br_result;


endmodule
