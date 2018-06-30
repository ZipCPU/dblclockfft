////////////////////////////////////////////////////////////////////////////////
//
// Filename:	fftmain.v
//
// Project:	A General Purpose Pipelined FFT Implementation
//
// Purpose:	This is the main module in the General Purpose FPGA FFT
//		implementation.  As such, all other modules are subordinate
//	to this one.  This module accomplish a fixed size Complex FFT on
//	16 data points.
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
//			only having 19 bits for each of the real and imaginary
//			components, leading to 38 bits total.
//	o_sync	A one bit output indicating the first sample of the FFT frame.
//			It also indicates the first valid sample out of the FFT
//			on the first frame.
//
// Arguments:	This file was computer generated using the following command
//		line:
//
//		% ./fftgen -v -d ../rtl -k 1 -f 16 -n 16 -p 0 -a ../bench/cpp/fftsize.h
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
	parameter	IWIDTH=16, OWIDTH=19, LGWIDTH=4;
	//
	input					i_clk, i_reset, i_ce;
	//
	input		[(2*IWIDTH-1):0]	i_sample;
	output	reg	[(2*OWIDTH-1):0]	o_result;
	output	reg				o_sync;


	// Outputs of the FFT, ready for bit reversal.
	wire	[(2*OWIDTH-1):0]	br_sample;




	wire		w_s16;
	wire	[33:0]	w_d16;
	fftstage	#(IWIDTH,IWIDTH+4,17,4,3,4,0,
			0, 11, 1, "cmem_16.hex")
		stage_16(i_clk, i_reset, i_ce,
			(!i_reset), i_sample, w_d16, w_s16);




	wire		w_s8;
	wire	[35:0]	w_d8;
	fftstage	#(17,21,18,4,2,4,0,
			0, 12, 1, "cmem_8.hex")
		stage_8(i_clk, i_reset, i_ce,
			w_s16, w_d16, w_d8, w_s8);


	wire		w_s4;
	wire	[35:0]	w_d4;
	qtrstage	#(18,18,4,0,0)	stage_4(i_clk, i_reset, i_ce,
						w_s8, w_d8, w_d4, w_s4);
	wire		w_s2;
	wire	[37:0]	w_d2;
	laststage	#(18,19,1)	stage_2(i_clk, i_reset, i_ce,
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
	bitreverse	#(4,19)
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
