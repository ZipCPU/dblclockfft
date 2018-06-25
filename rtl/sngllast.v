////////////////////////////////////////////////////////////////////////////////
//
// Filename:	sngllast.v
//
// Project:	A General Purpose Pipelined FFT Implementation
//
// Purpose:	This is part of an FPGA implementation that will process
//		the final stage of a decimate-in-frequency FFT, running
//	through the data at one sample per clock.
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
module	sngllast(i_clk, i_reset, i_ce, i_sync, i_val, o_val, o_sync);
	parameter	IWIDTH=16,OWIDTH=IWIDTH+1, SHIFT=0;
	input					i_clk, i_reset, i_ce, i_sync;
	input		[(2*IWIDTH-1):0]	i_val;
	output	wire	[(2*OWIDTH-1):0]	o_val;
	output	reg				o_sync;

	reg	signed	[(IWIDTH-1):0]	m_r, m_i;
	wire	signed	[(IWIDTH-1):0]	i_r, i_i;

	assign	i_r = i_val[(2*IWIDTH-1):(IWIDTH)]; 
	assign	i_i = i_val[(IWIDTH-1):0]; 

	// Don't forget that we accumulate a bit by adding two values
	// together. Therefore our intermediate value must have one more
	// bit than the two originals.
	reg	signed	[(IWIDTH):0]	rnd_r, rnd_i, sto_r, sto_i;
	reg				wait_for_sync, rnd_sync, stage, pre_sync;

	initial	rnd_sync      = 1'b0;
	initial	o_sync        = 1'b0;
	initial	wait_for_sync = 1'b1;
	initial	stage         = 1'b0;
	always @(posedge i_clk)
		if (i_reset)
		begin
			rnd_sync      <= 1'b0;
			o_sync        <= 1'b0;
			wait_for_sync <= 1'b1;
			stage         <= 1'b0;
		end else if ((i_ce)&&((!wait_for_sync)||(i_sync))&&(!stage))
		begin
			wait_for_sync <= 1'b0;
			//
			stage <= 1'b1;
			//
			o_sync <= rnd_sync;
		end else if (i_ce)
		begin
			rnd_sync <= pre_sync;
			//
			stage <= 1'b0;
			o_sync <= 1'b0;
		end
	always @(posedge i_clk)
	if (i_reset)
		pre_sync <= 1'b0;
	else if (i_ce)
		pre_sync <= i_sync;


	always @(posedge i_clk)
	if (i_ce)
	begin
		if (!stage)
		begin
			// Clock 1
			m_r <= i_r;
			m_i <= i_i;
			// Clock 3
			rnd_r <= sto_r;
			rnd_i <= sto_i;
			//
		end else begin
			// Clock 2
			rnd_r <= m_r + i_r;
			rnd_i <= m_i + i_i;
			//
			sto_r <= m_r - i_r;
			sto_i <= m_i - i_i;
			//
		end
	end

	// Now that we have our results, let's round them and report them
	wire	signed	[(OWIDTH-1):0]	o_r, o_i;

	convround #(IWIDTH+1,OWIDTH,SHIFT) do_rnd_r(i_clk, i_ce, rnd_r, o_r);
	convround #(IWIDTH+1,OWIDTH,SHIFT) do_rnd_i(i_clk, i_ce, rnd_i, o_i);

	assign	o_val  = { o_r, o_i };

endmodule
