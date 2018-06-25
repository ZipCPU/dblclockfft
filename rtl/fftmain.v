////////////////////////////////////////////////////////////////////////////////
//
// Filename:	fftmain.v
//
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
//			into two two's complement numbers, 16 bits each, with
//			the real portion in the high order bits, and the
//			imaginary portion taking the bottom 16 bits.
//	o_result	The output result, of the same format as i_sample,
//			only having 22 bits for each of the real and imaginary
//			components, leading to 44 bits total.
//	o_sync	A one bit output indicating the first sample of the FFT frame.
//			It also indicates the first valid sample out of the FFT
//			on the first frame.
//
// Arguments:	This file was computer generated using the following command
//		line:
//
//		% ./fftgen -v -d ../rtl -k 2 -f 2048 -n 16 -p 6 -a ../bench/cpp/fftsize.h
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
//
//
module fftmain(i_clk, i_reset, i_ce,
		i_sample, o_result, o_sync);
	parameter	IWIDTH=16, OWIDTH=22, LGWIDTH=11;
	//
	input					i_clk, i_reset, i_ce;
	//
	input		[(2*IWIDTH-1):0]	i_sample;
	output	reg	[(2*OWIDTH-1):0]	o_result;
	output	reg				o_sync;


	// Outputs of the FFT, ready for bit reversal.
	wire	[(2*OWIDTH-1):0]	br_sample;




	wire		w_s2048;
	wire	[33:0]	w_d2048;
	fftstage	#(IWIDTH,IWIDTH+4,17,11,10,4,0,
			1'b0, 2, "cmem_2048.hex")
		stage_e2048(i_clk, i_reset, i_ce,
			(!i_reset), i_sample, w_d2048, w_s2048);




	wire		w_s1024;
	wire	[35:0]	w_d1024;
	fftstage	#(17,21,18,11,9,4,0,
			1'b0, 2, "cmem_1024.hex")
		stage_1024(i_clk, i_reset, i_ce,
			w_s2048, w_d2048, w_d1024, w_s1024);


	wire		w_s512;
	wire	[35:0]	w_d512;
	fftstage	#(18,22,18,11,8,4,0,
			1'b0, 2, "cmem_512.hex")
		stage_512(i_clk, i_reset, i_ce,
			w_s1024, w_d1024, w_d512, w_s512);


	wire		w_s256;
	wire	[37:0]	w_d256;
	fftstage	#(18,22,19,11,7,4,0,
			1'b0, 2, "cmem_256.hex")
		stage_256(i_clk, i_reset, i_ce,
			w_s512, w_d512, w_d256, w_s256);


	wire		w_s128;
	wire	[37:0]	w_d128;
	fftstage	#(19,23,19,11,6,4,0,
			1'b0, 2, "cmem_128.hex")
		stage_128(i_clk, i_reset, i_ce,
			w_s256, w_d256, w_d128, w_s128);


	wire		w_s64;
	wire	[39:0]	w_d64;
	fftstage	#(19,23,20,11,5,4,0,
			1'b0, 2, "cmem_64.hex")
		stage_64(i_clk, i_reset, i_ce,
			w_s128, w_d128, w_d64, w_s64);


	// A hardware optimized FFT stage
	wire		w_s32;
	wire	[39:0]	w_d32;
	fftstage	#(20,24,20,11,4,4,0,
			1'b1, 2, "cmem_32.hex")
		stage_32(i_clk, i_reset, i_ce,
			w_s64, w_d64, w_d32, w_s32);


	// A hardware optimized FFT stage
	wire		w_s16;
	wire	[41:0]	w_d16;
	fftstage	#(20,24,21,11,3,4,0,
			1'b1, 2, "cmem_16.hex")
		stage_16(i_clk, i_reset, i_ce,
			w_s32, w_d32, w_d16, w_s16);


	// A hardware optimized FFT stage
	wire		w_s8;
	wire	[41:0]	w_d8;
	fftstage	#(21,25,21,11,2,5,0,
			1'b1, 2, "cmem_8.hex")
		stage_8(i_clk, i_reset, i_ce,
			w_s16, w_d16, w_d8, w_s8);


	wire		w_s4;
	wire	[43:0]	w_d4;
	qtrstage	#(21,22,11,0,0,0)	stage_4(i_clk, i_reset, i_ce,
						w_s8, w_d8, w_d4, w_s4);
	wire		w_s2;
	wire	[43:0]	w_d2;
	sngllast	#(22,22,0)	stage_2(i_clk, i_reset, i_ce,
					w_s4, w_d4, w_d2, w_s2);


	// Prepare for a (potential) bit-reverse stage.
	assign	br_sample= w_d2;

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
	wire	br_sync;
	wire	[(2*OWIDTH-1):0]	br_o_result;
	bitreverse	#(11,22)
		revstage(i_clk, i_reset,
			(i_ce & br_start), br_sample,
			br_o_result, br_sync);


	// Last clock: Register our outputs, we're done.
	initial	o_sync  = 1'b0;
	always @(posedge i_clk)
		if (i_reset)
			o_sync  <= 1'b0;
		else if (i_ce)
			o_sync  <= br_sync;

	always @(posedge i_clk)
		if (i_ce)
			o_result  <= br_o_result;


endmodule
