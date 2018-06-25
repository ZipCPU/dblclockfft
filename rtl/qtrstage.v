////////////////////////////////////////////////////////////////////////////////
//
// Filename:	qtrstage.v
//
// Project:	A General Purpose Pipelined FFT Implementation
//
// Purpose:	This file encapsulates the 4 point stage of a decimation in
//		frequency FFT.  This particular implementation is optimized
//	so that all of the multiplies are accomplished by additions and
//	multiplexers only.
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
module	qtrstage(i_clk, i_reset, i_ce, i_sync, i_data, o_data, o_sync);
	parameter	IWIDTH=16, OWIDTH=IWIDTH+1;
	// Parameters specific to the core that should be changed when this
	// core is built ... Note that the minimum LGSPAN is 2.  Smaller
	// spans must use the fftdoubles stage.
	parameter	LGWIDTH=8, ODD=0, INVERSE=0,SHIFT=0;
	input					i_clk, i_reset, i_ce, i_sync;
	input		[(2*IWIDTH-1):0]	i_data;
	output	reg	[(2*OWIDTH-1):0]	o_data;
	output	reg				o_sync;
	
	reg		wait_for_sync;
	reg	[3:0]	pipeline;

	reg	[(IWIDTH):0]	sum_r, sum_i, diff_r, diff_i;

	reg	[(2*OWIDTH-1):0]	ob_a;
	wire	[(2*OWIDTH-1):0]	ob_b;
	reg	[(OWIDTH-1):0]		ob_b_r, ob_b_i;
	assign	ob_b = { ob_b_r, ob_b_i };

	reg	[(LGWIDTH-1):0]		iaddr;
	reg	[(2*IWIDTH-1):0]	imem;

	wire	signed	[(IWIDTH-1):0]	imem_r, imem_i;
	assign	imem_r = imem[(2*IWIDTH-1):(IWIDTH)];
	assign	imem_i = imem[(IWIDTH-1):0];

	wire	signed	[(IWIDTH-1):0]	i_data_r, i_data_i;
	assign	i_data_r = i_data[(2*IWIDTH-1):(IWIDTH)];
	assign	i_data_i = i_data[(IWIDTH-1):0];

	reg	[(2*OWIDTH-1):0]	omem;

	wire	signed	[(OWIDTH-1):0]	rnd_sum_r, rnd_sum_i, rnd_diff_r, rnd_diff_i,
					n_rnd_diff_r, n_rnd_diff_i;
	convround #(IWIDTH+1,OWIDTH,SHIFT)	do_rnd_sum_r(i_clk, i_ce,
				sum_r, rnd_sum_r);

	convround #(IWIDTH+1,OWIDTH,SHIFT)	do_rnd_sum_i(i_clk, i_ce,
				sum_i, rnd_sum_i);

	convround #(IWIDTH+1,OWIDTH,SHIFT)	do_rnd_diff_r(i_clk, i_ce,
				diff_r, rnd_diff_r);

	convround #(IWIDTH+1,OWIDTH,SHIFT)	do_rnd_diff_i(i_clk, i_ce,
				diff_i, rnd_diff_i);

	assign n_rnd_diff_r = - rnd_diff_r;
	assign n_rnd_diff_i = - rnd_diff_i;
	initial wait_for_sync = 1'b1;
	initial iaddr = 0;
	always @(posedge i_clk)
		if (i_reset)
		begin
			wait_for_sync <= 1'b1;
			iaddr <= 0;
		end else if ((i_ce)&&((!wait_for_sync)||(i_sync)))
		begin
			iaddr <= iaddr + { {(LGWIDTH-1){1'b0}}, 1'b1 };
			wait_for_sync <= 1'b0;
		end
	always @(posedge i_clk)
		if (i_ce)
			imem <= i_data;


	// Note that we don't check on wait_for_sync or i_sync here.
	// Why not?  Because iaddr will always be zero until after the
	// first i_ce, so we are safe.
	initial pipeline = 4'h0;
	always	@(posedge i_clk)
		if (i_reset)
			pipeline <= 4'h0;
		else if (i_ce) // is our pipeline process full?  Which stages?
			pipeline <= { pipeline[2:0], iaddr[0] };

	// This is the pipeline[-1] stage, pipeline[0] will be set next.
	always	@(posedge i_clk)
		if ((i_ce)&&(iaddr[0]))
		begin
			sum_r  <= imem_r + i_data_r;
			sum_i  <= imem_i + i_data_i;
			diff_r <= imem_r - i_data_r;
			diff_i <= imem_i - i_data_i;
		end

	// pipeline[1] takes sum_x and diff_x and produces rnd_x

	// Now for pipeline[2].  We can actually do this at all i_ce
	// clock times, since nothing will listen unless pipeline[3]
	// on the next clock.  Thus, we simplify this logic and do
	// it independent of pipeline[2].
	always	@(posedge i_clk)
		if (i_ce)
		begin
			ob_a <= { rnd_sum_r, rnd_sum_i };
			// on Even, W = e^{-j2pi 1/4 0} = 1
			if (ODD == 0)
			begin
				ob_b_r <= rnd_diff_r;
				ob_b_i <= rnd_diff_i;
			end else if (INVERSE==0) begin
				// on Odd, W = e^{-j2pi 1/4} = -j
				ob_b_r <=   rnd_diff_i;
				ob_b_i <= n_rnd_diff_r;
			end else begin
				// on Odd, W = e^{j2pi 1/4} = j
				ob_b_r <= n_rnd_diff_i;
				ob_b_i <=   rnd_diff_r;
			end
		end

	always	@(posedge i_clk)
		if (i_ce)
		begin // In sequence, clock = 3
			if (pipeline[3])
			begin
				omem <= ob_b;
				o_data <= ob_a;
			end else
				o_data <= omem;
		end

	// Don't forget in the sync check that we are running
	// at two clocks per sample.  Thus we need to
	// produce a sync every 2^(LGWIDTH-1) clocks.
	initial	o_sync = 1'b0;
	always	@(posedge i_clk)
		if (i_reset)
			o_sync <= 1'b0;
		else if (i_ce)
			o_sync <= &(~iaddr[(LGWIDTH-2):3]) && (iaddr[2:0] == 3'b101);
endmodule
