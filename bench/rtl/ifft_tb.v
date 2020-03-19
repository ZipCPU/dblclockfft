////////////////////////////////////////////////////////////////////////////////
//
// Filename: 	ifft_tb.v
//
// Project:	A General Purpose Pipelined FFT Implementation
//
// Purpose:	This file is used to test whether an FFT followed by an
//		immediate IFFT produces the input again.  It is a test
//		bench for the composition of the two programs, and specifically
//		for the IFFT.  As designed and built, this test bench would
//		be difficult to implement in an FPGA (although I wouldn't
//		put it past a smart individual), it was designed to be
//		used within Verilator as part of a C++ simulation.  The
//		simulation source may be found in this project inside
//		ifft_tb.cpp.
//
// Creator:	Dan Gisselquist, Ph.D.
//		Gisselquist Technology, LLC
//
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2015-2020, Gisselquist Technology, LLC
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
// with this program.  (It's in the $(ROOT)/doc directory.  Run make with no
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
`default_nettype none
//
// `define	DBLCLK
//
module	ifft_tb(i_clk, i_rst, i_ce,
`ifdef	DBLCLK
		i_left, i_right, o_left, o_right,
`else
		i_sample, o_sample,
`endif
		o_sync);
	parameter	IWIDTH=15, MIDWIDTH=21, OWIDTH=21;
	input					i_clk, i_rst, i_ce;
`ifdef	DBLCLK
	input		[(2*IWIDTH-1):0]	i_left, i_right;
	output	wire	[(2*OWIDTH-1):0]	o_left, o_right;
`else
	input		[(2*IWIDTH-1):0]	i_sample;
	output	wire	[(2*OWIDTH-1):0]	o_sample;
`endif
	output	wire				o_sync;

	wire				m_sync;
`ifdef	DBLCLK
	wire	[(2*MIDWIDTH-1):0]	m_left, m_right;
	fftmain	fft(i_clk, i_rst, i_ce, i_left, i_right,
				m_left, m_right, m_sync);
`else
	wire	[(2*MIDWIDTH-1):0]	m_sample;
	fftmain	fft(i_clk, i_rst, i_ce, i_sample,
				m_sample, m_sync);
`endif

	wire	w_syncd;
	reg	r_syncd;
	always @(posedge i_clk)
	if (i_rst)
		r_syncd <= 1'b0;
	else
		r_syncd <= r_syncd || m_sync;

	assign	w_syncd = r_syncd || m_sync;

`ifdef	DBLCLK
	ifftmain	ifft(i_clk, i_rst, (i_ce)&&(w_syncd), m_left, m_right,
				o_left, o_right, o_sync);
`else
	ifftmain	ifft(i_clk, i_rst, (i_ce)&&(w_syncd),
				{ m_sample[2*MIDWIDTH-1:2*MIDWIDTH-IWIDTH],
					m_sample[MIDWIDTH-1:MIDWIDTH-IWIDTH] },
				o_sample, o_sync);

	wire unused;
	assign unused = &{ 1'b0, m_sample[2*MIDWIDTH-IWIDTH-1:MIDWIDTH],
			m_sample[MIDWIDTH-IWIDTH-1:0] };
`endif
endmodule
