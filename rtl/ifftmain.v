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
//	The FFT is fully pipelined, and accepts as inputs two complex two's
//	complement samples per clock.
//
// Parameters:
//	i_clk	The clock.  All operations are synchronous with this clock.
//	i_reset	Synchronous reset, active high.  Setting this line will
//			force the reset of all of the internals to this routine.
//			Further, following a reset, the o_sync line will go
//			high the same time the first output sample is valid.
//	i_ce	A clock enable line.  If this line is set, this module
//			will accept two complex values as inputs, and produce
//			two (possibly empty) complex values as outputs.
//	i_left	The first of two complex input samples.  This value is split
//			into two two's complement numbers, 15 bits each, with
//			the real portion in the high order bits, and the
//			imaginary portion taking the bottom 15 bits.
//	i_right	This is the same thing as i_left, only this is the second of
//			two such samples.  Hence, i_left would contain input
//			sample zero, i_right would contain sample one.  On the
//			next clock i_left would contain input sample two,
//			i_right number three and so forth.
//	o_left	The first of two output samples, of the same format as i_left,
//			only having 21 bits for each of the real and imaginary
//			components, leading to 42 bits total.
//	o_right	The second of two output samples produced each clock.  This has
//			the same format as o_left.
//	o_sync	A one bit output indicating the first valid sample produced by
//			this FFT following a reset.  Ever after, this will
//			indicate the first sample of an FFT frame.
//
// Arguments:	This file was computer generated using the following command
//		line:
//
//		% ./fftgen -i -d ../rtl -f 2048 -2 -p 0 -n 15 -a ../bench/cpp/ifftsize.h
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
		i_left, i_right,
		o_left, o_right, o_sync);
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
	input	wire	[(2*IWIDTH-1):0]	i_left, i_right;
	output	reg	[(2*OWIDTH-1):0]	o_left, o_right;
	output	reg				o_sync;


	// Outputs of the FFT, ready for bit reversal.
	wire				br_sync;
	wire	[(2*OWIDTH-1):0]	br_left, br_right;


	wire		w_s2048;
	// verilator lint_off UNUSED
	wire		w_os2048;
	// verilator lint_on  UNUSED
	wire	[31:0]	w_e2048, w_o2048;
	fftstage	#(IWIDTH,IWIDTH+4,16,9,0,
			0, 1, "icmem_e4096.hex")
		stage_e2048(i_clk, i_reset, i_ce,
			(!i_reset), i_left, w_e2048, w_s2048);
	fftstage	#(IWIDTH,IWIDTH+4,16,9,0,
			0, 1, "icmem_o4096.hex")
		stage_o2048(i_clk, i_reset, i_ce,
			(!i_reset), i_right, w_o2048, w_os2048);


	wire		w_s1024;
	// verilator lint_off UNUSED
	wire		w_os1024;
	// verilator lint_on  UNUSED
	wire	[33:0]	w_e1024, w_o1024;
	fftstage	#(16,20,17,8,0,
			0, 1, "icmem_e2048.hex")
		stage_e1024(i_clk, i_reset, i_ce,
			w_s2048, w_e2048, w_e1024, w_s1024);
	fftstage	#(16,20,17,8,0,
			0, 1, "icmem_o2048.hex")
		stage_o1024(i_clk, i_reset, i_ce,
			w_s2048, w_o2048, w_o1024, w_os1024);

	wire		w_s512;
	// verilator lint_off UNUSED
	wire		w_os512;
	// verilator lint_on  UNUSED
	wire	[33:0]	w_e512, w_o512;
	fftstage	#(17,21,17,7,0,
			0, 1, "icmem_e1024.hex")
		stage_e512(i_clk, i_reset, i_ce,
			w_s1024, w_e1024, w_e512, w_s512);
	fftstage	#(17,21,17,7,0,
			0, 1, "icmem_o1024.hex")
		stage_o512(i_clk, i_reset, i_ce,
			w_s1024, w_o1024, w_o512, w_os512);

	wire		w_s256;
	// verilator lint_off UNUSED
	wire		w_os256;
	// verilator lint_on  UNUSED
	wire	[35:0]	w_e256, w_o256;
	fftstage	#(17,21,18,6,0,
			0, 1, "icmem_e512.hex")
		stage_e256(i_clk, i_reset, i_ce,
			w_s512, w_e512, w_e256, w_s256);
	fftstage	#(17,21,18,6,0,
			0, 1, "icmem_o512.hex")
		stage_o256(i_clk, i_reset, i_ce,
			w_s512, w_o512, w_o256, w_os256);

	wire		w_s128;
	// verilator lint_off UNUSED
	wire		w_os128;
	// verilator lint_on  UNUSED
	wire	[35:0]	w_e128, w_o128;
	fftstage	#(18,22,18,5,0,
			0, 1, "icmem_e256.hex")
		stage_e128(i_clk, i_reset, i_ce,
			w_s256, w_e256, w_e128, w_s128);
	fftstage	#(18,22,18,5,0,
			0, 1, "icmem_o256.hex")
		stage_o128(i_clk, i_reset, i_ce,
			w_s256, w_o256, w_o128, w_os128);

	wire		w_s64;
	// verilator lint_off UNUSED
	wire		w_os64;
	// verilator lint_on  UNUSED
	wire	[37:0]	w_e64, w_o64;
	fftstage	#(18,22,19,4,0,
			0, 1, "icmem_e128.hex")
		stage_e64(i_clk, i_reset, i_ce,
			w_s128, w_e128, w_e64, w_s64);
	fftstage	#(18,22,19,4,0,
			0, 1, "icmem_o128.hex")
		stage_o64(i_clk, i_reset, i_ce,
			w_s128, w_o128, w_o64, w_os64);

	wire		w_s32;
	// verilator lint_off UNUSED
	wire		w_os32;
	// verilator lint_on  UNUSED
	wire	[37:0]	w_e32, w_o32;
	fftstage	#(19,23,19,3,0,
			0, 1, "icmem_e64.hex")
		stage_e32(i_clk, i_reset, i_ce,
			w_s64, w_e64, w_e32, w_s32);
	fftstage	#(19,23,19,3,0,
			0, 1, "icmem_o64.hex")
		stage_o32(i_clk, i_reset, i_ce,
			w_s64, w_o64, w_o32, w_os32);

	wire		w_s16;
	// verilator lint_off UNUSED
	wire		w_os16;
	// verilator lint_on  UNUSED
	wire	[39:0]	w_e16, w_o16;
	fftstage	#(19,23,20,2,0,
			0, 1, "icmem_e32.hex")
		stage_e16(i_clk, i_reset, i_ce,
			w_s32, w_e32, w_e16, w_s16);
	fftstage	#(19,23,20,2,0,
			0, 1, "icmem_o32.hex")
		stage_o16(i_clk, i_reset, i_ce,
			w_s32, w_o32, w_o16, w_os16);

	wire		w_s8;
	// verilator lint_off UNUSED
	wire		w_os8;
	// verilator lint_on  UNUSED
	wire	[39:0]	w_e8, w_o8;
	fftstage	#(20,24,20,1,0,
			0, 1, "icmem_e16.hex")
		stage_e8(i_clk, i_reset, i_ce,
			w_s16, w_e16, w_e8, w_s8);
	fftstage	#(20,24,20,1,0,
			0, 1, "icmem_o16.hex")
		stage_o8(i_clk, i_reset, i_ce,
			w_s16, w_o16, w_o8, w_os8);

	wire		w_s4;
	// verilator lint_off UNUSED
	wire		w_os4;
	// verilator lint_on  UNUSED
	wire	[41:0]	w_e4, w_o4;
	qtrstage	#(20,21,11,0,1,0)	stage_e4(i_clk, i_reset, i_ce,
						w_s8, w_e8, w_e4, w_s4);
	qtrstage	#(20,21,11,1,1,0)	stage_o4(i_clk, i_reset, i_ce,
						w_s8, w_o8, w_o4, w_os4);
	wire		w_s2;
	wire	[41:0]	w_e2, w_o2;
	laststage	#(21,21,0)	stage_2(i_clk, i_reset, i_ce,
					w_s4, w_e4, w_o4, w_e2, w_o2, w_s2);


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
	bitreverse	#(11,21)
	revstage(
		// {{{
		.i_clk(i_clk),
		.i_reset(i_reset),
		.i_ce(i_ce & br_start),
		.i_in_0(w_e2),
		.i_in_1(w_o2),
		.o_out_0(br_left),
		.o_out_1(br_right),
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
	begin
		o_left  <= br_left;
		o_right <= br_right;
	end


endmodule
