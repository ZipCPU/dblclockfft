////////////////////////////////////////////////////////////////////////////////
//
// Filename:	bitreverse.v
// {{{
// Project:	A General Purpose Pipelined FFT Implementation
//
// Purpose:	This module bitreverses a pipelined FFT input.  Operation is
//		expected as follows:
//
//		i_clk	A running clock at whatever system speed is offered.
//		i_reset	A synchronous reset signal, that resets all internals
//		i_ce	If this is one, one input is consumed and an output
//			is produced.
//		i_in_0, i_in_1
//			Two inputs to be consumed, each of width WIDTH.
//		o_out_0, o_out_1
//			Two of the bitreversed outputs, also of the same
//			width, WIDTH.  Of course, there is a delay from the
//			first input to the first output.  For this purpose,
//			o_sync is present.
//		o_sync	This will be a 1'b1 for the first value in any block.
//			Following a reset, this will only become 1'b1 once
//			the data has been loaded and is now valid.  After that,
//			all outputs will be valid.
//
// How do we do bit reversing at two smples per clock?  Can we separate out
// our work into eight memory banks, writing two banks at once and reading
// another two banks in the same clock?
//
//	mem[00xxx0] = s_0[n]
//	mem[00xxx1] = s_1[n]
//	o_0[n] = mem[10xxx0]
//	o_1[n] = mem[11xxx0]
//	...
//	mem[01xxx0] = s_0[m]
//	mem[01xxx1] = s_1[m]
//	o_0[m] = mem[10xxx1]
//	o_1[m] = mem[11xxx1]
//	...
//	mem[10xxx0] = s_0[n]
//	mem[10xxx1] = s_1[n]
//	o_0[n] = mem[00xxx0]
//	o_1[n] = mem[01xxx0]
//	...
//	mem[11xxx0] = s_0[m]
//	mem[11xxx1] = s_1[m]
//	o_0[m] = mem[00xxx1]
//	o_1[m] = mem[01xxx1]
//	...
//
//	The answer is that, yes we can but: we need to use four memory banks
//	to do it properly.  These four banks are defined by the two bits
//	that determine the top and bottom of the correct address.  Larger
//	FFT's would require more memories.
//
//
//	20150602 -- This module has undergone massive rework in order to
//		ensure that it uses resources efficiently.  As a result,
//		it now optimizes nicely into block RAMs.  As an unfortunately
//		side effect, it now passes it's bench test (dblrev_tb) but
//		fails the integration bench test (fft_tb).
//
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
module	bitreverse #(
		// {{{
		parameter			LGSIZE=5, WIDTH=24
		// }}}
	) (
		// {{{
		input	wire			i_clk, i_reset, i_ce,
		input	wire	[(2*WIDTH-1):0]	i_in_0, i_in_1,
		output	wire	[(2*WIDTH-1):0]	o_out_0, o_out_1,
		output	reg			o_sync
		// }}}
	);

	// Local declarations
	// {{{
	reg			in_reset;
	reg	[(LGSIZE-1):0]	iaddr;
	wire	[(LGSIZE-3):0]	braddr;

	reg	[(2*WIDTH-1):0]	mem_e [0:((1<<(LGSIZE))-1)];
	reg	[(2*WIDTH-1):0]	mem_o [0:((1<<(LGSIZE))-1)];

	reg [(2*WIDTH-1):0] evn_out_0, evn_out_1, odd_out_0, odd_out_1;
	reg	adrz;
	// }}}

	// braddr
	// {{{
	genvar	k;
	generate for(k=0; k<LGSIZE-2; k=k+1)
	begin : gen_a_bit_reversed_value
		assign braddr[k] = iaddr[LGSIZE-3-k];
	end endgenerate
	// }}}

	// iaddr, in_reset, o_sync
	// {{{
	initial iaddr = 0;
	initial in_reset = 1'b1;
	initial o_sync = 1'b0;
	always @(posedge i_clk)
	if (i_reset)
	begin
		iaddr <= 0;
		in_reset <= 1'b1;
		o_sync <= 1'b0;
	end else if (i_ce)
	begin
		iaddr <= iaddr + { {(LGSIZE-1){1'b0}}, 1'b1 };
		if (&iaddr[(LGSIZE-2):0])
			in_reset <= 1'b0;
		if (in_reset)
			o_sync <= 1'b0;
		else
			o_sync <= ~(|iaddr[(LGSIZE-2):0]);
	end
	// }}}

	// Write to memories mem_e and mem_o
	// {{{
	always @(posedge i_clk)
	if (i_ce)
		mem_e[iaddr] <= i_in_0;

	always @(posedge i_clk)
	if (i_ce)
		mem_o[iaddr] <= i_in_1;
	// }}}

	// Read from memories into: [evn|odd]_out_[0|1]
	// {{{
	always @(posedge i_clk)
	if (i_ce)
		evn_out_0 <= mem_e[{!iaddr[LGSIZE-1],1'b0,braddr}];

	always @(posedge i_clk)
	if (i_ce)
		evn_out_1 <= mem_e[{!iaddr[LGSIZE-1],1'b1,braddr}];

	always @(posedge i_clk)
	if (i_ce)
		odd_out_0 <= mem_o[{!iaddr[LGSIZE-1],1'b0,braddr}];

	always @(posedge i_clk)
	if (i_ce)
		odd_out_1 <= mem_o[{!iaddr[LGSIZE-1],1'b1,braddr}];
	// }}}

	// adrz
	// {{{
	always @(posedge i_clk)
	if (i_ce)
		adrz <= iaddr[LGSIZE-2];
	// }}}

	assign	o_out_0 = (adrz)?odd_out_0:evn_out_0;
	assign	o_out_1 = (adrz)?odd_out_1:evn_out_1;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
// Formal property section
// {{{
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
`ifdef	FORMAL
	// Formal declarations
	// {{{
`define	ASSERT	assert
`ifdef	BITREVERSE
`define	ASSUME	assume
`else
`define	ASSUME	assert
`endif

	reg	f_past_valid;
	(* anyconst *) reg	[LGSIZE-1:0]	f_const_addr;
	wire	[LGSIZE-3:0]	f_reversed_addr;
	// reg	[LGSIZE:0]	f_now;
	reg			f_addr_loaded_0, f_addr_loaded_1;
	reg	[(2*WIDTH-1):0]	f_data_0, f_data_1;
	wire			f_writing, f_reading;
	// }}}

	initial	f_past_valid = 1'b0;
	always @(posedge i_clk)
		f_past_valid <= 1'b1;

	initial	`ASSUME(i_reset);
	always @(posedge i_clk)
	if ((!f_past_valid)||($past(i_reset)))
	begin
		`ASSERT(iaddr == 0);
		`ASSERT(in_reset);
		`ASSERT(!o_sync);
	end
`ifdef	BITREVERSE
	always @(posedge i_clk)
		assume((i_ce)||($past(i_ce))||($past(i_ce,2)));
`endif // BITREVERSE

	// f_reversed_addr
	// {{{
	generate for(k=0; k<LGSIZE-2; k=k+1)
		assign	f_reversed_addr[k] = f_const_addr[LGSIZE-3-k];
	endgenerate
	// }}}

		assign	f_writing=(f_const_addr[LGSIZE-1]==iaddr[LGSIZE-1]);
		assign	f_reading=(f_const_addr[LGSIZE-1]!=iaddr[LGSIZE-1]);
	// f_addr_loaded_[0|1]
	// {{{
	initial	f_addr_loaded_0 = 1'b0;
	initial	f_addr_loaded_1 = 1'b0;
	always @(posedge i_clk)
	if (i_reset)
	begin
		f_addr_loaded_0 <= 1'b0;
		f_addr_loaded_1 <= 1'b0;
	end else if (i_ce)
	begin
		if (iaddr == f_const_addr)
		begin
			f_addr_loaded_0 <= 1'b1;
			f_addr_loaded_1 <= 1'b1;
		end

		if (f_reading)
		begin
			if ((braddr == f_const_addr[LGSIZE-3:0])
				&&(iaddr[LGSIZE-2] == 1'b0))
				f_addr_loaded_0 <= 1'b0;

			if ((braddr == f_const_addr[LGSIZE-3:0])
				&&(iaddr[LGSIZE-2] == 1'b1))
				f_addr_loaded_1 <= 1'b0;
		end
	end
	// }}}

	// f_data_0, f_data_1
	// {{{
	always @(posedge i_clk)
	if ((i_ce)&&(iaddr == f_const_addr))
	begin
		f_data_0 <= i_in_0;
		f_data_1 <= i_in_1;
		`ASSERT(!f_addr_loaded_0);
		`ASSERT(!f_addr_loaded_1);
	end
	// }}}

	always @(posedge i_clk)
	if ((f_past_valid)&&(!$past(i_reset))
			&&($past(f_addr_loaded_0))&&(!f_addr_loaded_0))
	begin
		assert(!$past(iaddr[LGSIZE-2]));
		if (f_const_addr[LGSIZE-2])
			assert(o_out_1 == f_data_0);
		else
			assert(o_out_0 == f_data_0);
	end

	always @(posedge i_clk)
	if ((f_past_valid)&&(!$past(i_reset))
			&&($past(f_addr_loaded_1))&&(!f_addr_loaded_1))
	begin
		assert($past(iaddr[LGSIZE-2]));
		if (f_const_addr[LGSIZE-2])
			assert(o_out_1 == f_data_1);
		else
			assert(o_out_0 == f_data_1);
	end

	always @(*)
		`ASSERT(o_sync == ((iaddr[LGSIZE-2:0] == 1)&&(!in_reset)));

	// Before writing to a section, the loaded flags should be
	// zero
	always @(*)
	if (f_writing)
	begin
		`ASSERT(f_addr_loaded_0 == (iaddr[LGSIZE-2:0]
					> f_const_addr[LGSIZE-2:0]));
		`ASSERT(f_addr_loaded_1 == (iaddr[LGSIZE-2:0]
					> f_const_addr[LGSIZE-2:0]));
	end

	// If we were writing, and now we are reading, then both
	// f_addr_loaded flags must be set
	always @(posedge i_clk)
	if ((f_past_valid)&&(!$past(i_reset))
			&&($past(f_writing))&&(f_reading))
	begin
		`ASSERT(f_addr_loaded_0);
		`ASSERT(f_addr_loaded_1);
	end

	always @(*)
	if (f_writing)
		`ASSERT(f_addr_loaded_0 == f_addr_loaded_1);

	// When reading, and the loaded flag is zero, our pointer
	// must not have hit the address of interest yet
	always @(*)
	if ((!in_reset)&&(f_reading))
		`ASSERT(f_addr_loaded_0 ==
			((!iaddr[LGSIZE-2])&&(iaddr[LGSIZE-3:0]
				<= f_reversed_addr[LGSIZE-3:0])));

	always @(*)
	if ((!in_reset)&&(f_reading))
		`ASSERT(f_addr_loaded_1 ==
			((!iaddr[LGSIZE-2])||(iaddr[LGSIZE-3:0]
				<= f_reversed_addr[LGSIZE-3:0])));

	always @(*)
	if ((in_reset)&&(f_reading))
	begin
		`ASSERT(!f_addr_loaded_0);
		`ASSERT(!f_addr_loaded_1);
	end

	always @(*)
	if(iaddr[LGSIZE-1])
		`ASSERT(!in_reset);

	always @(*)
	if (f_addr_loaded_0)
		`ASSERT(mem_e[f_const_addr] == f_data_0);
	always @(*)
	if (f_addr_loaded_1)
		`ASSERT(mem_o[f_const_addr] == f_data_1);

`endif	// FORMAL
// }}}
endmodule
