////////////////////////////////////////////////////////////////////////////////
//
// Filename:	windowfn.v
//
// Project:	A General Purpose Pipelined FFT Implementation
//
// Purpose:	Apply a window function to incoming real data points, so as
//		to create an outgoing stream of data samples that can be used
//	in an FFT construct using 50% overlap.  The overlap, coupled with the
//	FFT's requirements, can make for somewhat of a problem.  Hence, there
//	are two 'ce' signals coming into the core.  A primary ce signal when
//	new data is ready, and an alternate that must take place between
//	primary signals.  This allows the second/alternate CE signal to be
//	appropriately spaced between the primary CE signals so that the
//	outgoing signals to the FFT will still meet separation
//	requirements--whatever they would be for the application.
//
//	For this module, the window size is the FFT length.
//
// Ports:
//	i_clk, i_reset	Should be self explanatory.  The reset is assumed to
//			be synchronous.
//
//	i_tap_wr, i_tap	For use when OPT_FIXED_TAPS is zero, i_tap_wr signals
//			that a "tap" or "coefficient" of the filter should be
//		written.  When i_tap_wr is high, i_tap is taken to be
//		a coefficient to the core.  There's an internal address
//		counter, so no address need be given.  However, the counter is
//		reset on an i_reset signal.
//
//	i_ce, i_alt_ce	As discussed above, these signals need to alternate
//		back and forth.  Following a reset, the first signal coming in
//		should be an i_ce signal.
//
//	i_sample	The incoming sample data, valid any time i_ce is true,
//			and accepted into the core on that clock tick.
//
//	o_ce		True when the core has a valid output
//	o_sample	The output calculated by the core, ready to pass to
//			the FFT
//	o_frame		True on the first sample of any frame.  Following a
//			reset, o_ce will remain false until o_frame is also
//		true with it.  From then on out, o_frame will be true once
//		every FFT length.
//
//	For a timing/signaling diagram, please feel free to run the formal
//	tools in cover mode for this module, 'sby -f windowfn.sby cover',
//	and then check out the generated trace.
//
//
// Creator:	Dan Gisselquist, Ph.D.
//		Gisselquist Technology, LLC
//
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2018, Gisselquist Technology, LLC
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
`default_nettype	none
//
module	windowfn(i_clk, i_reset, i_tap_wr, i_tap,
		i_ce, i_sample, i_alt_ce,
		o_frame, o_ce, o_sample);
	parameter		IW=16, OW=16, TW=16, LGNFFT = 4;
	parameter	[0:0]	OPT_FIXED_TAPS = 1'b0;
	parameter		INITIAL_COEFFS = "";
	//
	localparam	AW=IW+TW;
	//
	input	wire			i_clk, i_reset;
	//
	input	wire			i_tap_wr;
	input	wire	[(TW-1):0]	i_tap;
	//
	input	wire			i_ce;
	input	wire	[(IW-1):0]	i_sample;
	input	wire			i_alt_ce;
	//
	output	reg			o_frame, o_ce;
	output	reg	[(OW-1):0]	o_sample;


	reg	[(TW-1):0]	cmem	[0:(1<<LGNFFT)-1];
	reg	[(TW-1):0]	dmem	[0:(1<<LGNFFT)-1];

	//
	// LOAD THE TAPS
	//
	wire	[LGNFFT-1:0]	tapwidx;
	generate if (OPT_FIXED_TAPS)
	begin : SET_FIXED_TAPS

		initial $readmemh(INITIAL_COEFFS, cmem);

		assign	tapwidx = 0;
		// Make Verilators -Wall happy
		// Verilator lint_off UNUSED
		wire	[TW:0]	ignored_inputs;
		assign	ignored_inputs = { i_tap_wr, i_tap };
		// Verilator lint_on  UNUSED

	end else begin : DYNAMICALLY_SET_TAPS

		// Coef memory write index
		reg	[(LGNFFT-1):0]	r_tapwidx;

		initial	r_tapwidx = 0;
		always @(posedge i_clk)
			if(i_reset)
				r_tapwidx <= 0;
			else if (i_tap_wr)
				r_tapwidx <= r_tapwidx + 1'b1;

		if (INITIAL_COEFFS != 0)
			initial $readmemh(INITIAL_COEFFS, cmem);
		always @(posedge i_clk)
			if (i_tap_wr)
				cmem[r_tapwidx] <= i_tap;

		assign	tapwidx = r_tapwidx;
	end endgenerate


	reg		[LGNFFT-1:0]	dwidx, didx;
	reg		[LGNFFT-1:0]	tidx;
	reg				top_of_block, first_block;
	reg		[1:0]		frame_count;
	reg				p_ce, d_ce;
	reg	signed	[IW-1:0]	data;
	reg	signed	[TW-1:0]	tap;



	//
	//
	// Record the incoming data into a local memory
	//
	//

	// Notice how this data writing section is *independent* of the reset,
	// depending only upon new sample data.

	initial	dwidx = 0;
	always @(posedge i_clk)
	if (i_reset)
		dwidx <= 0;
	else if (i_ce)
		dwidx <= dwidx + 1'b1;
	always @(posedge i_clk)
	if (i_ce)
		dmem[dwidx] <= i_sample;

	initial	first_block = 1'b1;
	always @(posedge i_clk)
	if (i_reset)
		first_block <= 1'b1;
	else if ((i_alt_ce)&&(&tidx)&&(dwidx==0))
		first_block <= 1'b0;


	//
	//
	// Keep track of the top of the block.  The top of the block is the
	// first incoming data sample on an FFT or half FFT boundary.  This
	// is where data processing starts from.
	//
	initial	top_of_block = 1'b0;
	always @(posedge i_clk)
	if (i_reset)
		top_of_block <= 1'b0;
	else if (i_alt_ce)
		top_of_block <= (&tidx)&&((!first_block)||(dwidx==0));
	else if (i_ce)
		top_of_block <= 1'b0;

	//
	// Data and coefficient memory indices.
	//
	initial	didx = 0;
	always @(posedge i_clk)
	if (i_reset)
		didx <= 0;
	else if ((i_alt_ce)&&(dwidx[LGNFFT-2:0]==0))
	begin
		didx[LGNFFT-2:0] <= 0;
		didx[LGNFFT-1] <= dwidx[LGNFFT-1];
	end else if ((i_ce)||(i_alt_ce))
		// Process the next point in this FFT
		didx <= didx + 1'b1;

	initial	tidx = 0;
	always @(posedge i_clk)
	if (i_reset)
		tidx <= 0;
	else if ((i_alt_ce)&&(dwidx[LGNFFT-2:0]==0))
	begin
		// // At the beginning of processing for a given FFT
		tidx <= 0;
	end else if ((i_ce)||(i_alt_ce))
		// Process the next point in the window function
		tidx <= tidx + 1'b1;

	initial	frame_count = 0;
	always @(posedge i_clk)
	if (i_reset)
		frame_count <= 0;
	else if ((i_ce)&&(top_of_block)&&(!first_block))
		frame_count <= 3;
	else if (frame_count != 0)
		frame_count <= frame_count - 1'b1;

	initial	o_frame = 1'b0;
	always @(posedge i_clk)
	if (i_reset)
		o_frame <= 1'b0;
	else if (frame_count == 2)
		o_frame <= 1'b1;
	else
		o_frame <= 1'b0;

	//
	// Following any initial i_ce, ...
	// 	d_ce: The data (and coefficient), read from memory,will be valid
	// 	p_ce: The produc of data and coefficient is valid
	//
	initial	{ p_ce, d_ce } = 2'h0;
	always @(posedge i_clk)
	if (i_reset)
		{ p_ce, d_ce } <= 2'h0;
	else
		{ p_ce, d_ce } <= { d_ce, (!first_block)&&((i_ce)||(i_alt_ce)) };

	// Read the data sample point, and the filter coefficient, from block
	// RAM.  Because this is block RAM, we have to be careful not to
	// do anything else here.
	initial	data = 0;
	initial	tap = 0;
	always @(posedge i_clk)
	begin
		data <= dmem[didx];
		tap  <= cmem[tidx];
	end

	//
	// Multiply the two values together
	 reg	signed	[IW+TW-1:0]	product;
`ifdef	FORMAL
	// We'll implement an abstract multiply below--just to make sure the
	// timing is right.
`else
	always @(posedge i_clk)
		product <= data * tap;
`endif

	//
	// Output CE
	//
	initial	o_ce = 0;
	always @(posedge i_clk)
	if (i_reset)
		o_ce <= 0;
	else
		o_ce <= p_ce;

	generate if (OW == AW)
	begin : BIT_ADJUSTMENT_NONE

		initial	o_sample = 0;
		always @(posedge i_clk)
		if (i_reset)
			o_sample <= 0;
		else if (p_ce)
			o_sample <= product;

	end else if (OW < AW)
	begin : BIT_ADJUSTMENT_ROUNDING
		wire	[AW-1:0]	rounded;

		assign	rounded = product + { {(OW){1'b0}}, product[AW-OW],
				{(AW-OW-1){!product[AW-OW]}} };

		initial	o_sample = 0;
		always @(posedge i_clk)
		if (i_reset)
			o_sample <= 0;
		else if (p_ce)
			o_sample <= rounded[(AW-1):(AW-OW)];

		// Make Verilator happy
		// verilator lint_off UNUSED
		wire	[AW-OW-1:0]	unused_rounding_bits;
		assign	unused_rounding_bits = rounded[AW-OW-1:0];
		// verilator lint_on  UNUSED

	end else // if (OW > AW)
	begin : BIT_ADJUSTMENT_EXTENDING

		always @(posedge i_clk)
		if (i_reset)
			o_sample <= 0;
		else if (p_ce)
			o_sample <= { product, {(OW-AW){1'b0}} };

	end endgenerate

	// Make Verilator happy
	// verilator lint_off UNUSED
	wire	[LGNFFT-1:0]	unused;
	assign	unused = tapwidx;
	// verilator lint_on  UNUSED

`ifdef	FORMAL
	reg	f_past_valid;
	initial	f_past_valid = 1'b0;
	always @(posedge i_clk)
		f_past_valid <= 1'b1;

	// Keep track of the phase of this operation
	reg	[LGNFFT:0]	f_phase;

	initial	f_phase = 0;
	always @(posedge i_clk)
	if (i_reset)
	begin
		f_phase <= 0;
	end else if ((i_ce)&&(top_of_block))
		f_phase[LGNFFT-1:0] <= 1;
	else if ((i_ce)||(i_alt_ce))
		f_phase <= f_phase + 1;

	///////
	//
	// Assumptions about the input
	always @(posedge i_clk)
	if ((f_past_valid)&&(!$past(i_reset))&&(!$past(i_ce)))
		restrict($stable(i_sample));

	always @(*)
	if (tapwidx != 0)
		assume(!i_ce);

	always @(posedge i_clk)
	if (f_phase[0])
		assume(!i_ce);

	always @(posedge i_clk)
	if (!f_phase[0])
		assume(!i_alt_ce);

	always @(*)
	if (i_tap_wr)
		assume((!i_ce)&&(!i_alt_ce));

	always @(*)
		assume(!i_tap_wr);

	always @(*)
	if ((i_ce)&&(top_of_block))
	begin
		assert(dwidx[LGNFFT-2:0]==0);
		assert(didx[LGNFFT-2:0]==0);
	end else if (dwidx[LGNFFT-2:0]!=0)
		assert(!top_of_block);

	always @(*)
	if ((dwidx[LGNFFT-2:0]==0)&&(!didx[0]))
		assert(top_of_block||f_first_block||first_block);
	else
		assert(!top_of_block);

	always @(posedge i_clk)
	if ((f_past_valid)&&($past(first_block))&&(!first_block))
		assert(top_of_block);

	/*
	always @(posedge i_clk)
	if (f_past_valid)
	begin
		if ($past(i_reset))
			assert(!top_of_block);
		else if ($past(i_ce))
			assert(top_of_block == ((!first_block)&&(dwidx[LGNFFT-2:0]==0)));
		else

			assert(top_of_block == ((!first_block)&&(&dwidx[LGNFFT-2:1])));
	end
	*/

	always @(posedge i_clk)
	if ((f_past_valid)&&($past(first_block))&&(!$past(first_block)))
		assert(top_of_block);

	//////////////////////////////////////////////////////////////////
	//
	// Assertions about our outputs
	//
	/////////////////////////
	//
	//
//	always @(*)
//	if (o_frame)
//		assert(o_ce);
	always @(*)
	if (first_block)
		assert(!o_ce);

	reg	f_waiting_for_first_frame;

	initial	f_waiting_for_first_frame = 1;
	always @(posedge i_clk)
	if ((i_reset)||(first_block))
		f_waiting_for_first_frame <= 1'b1;
	else if (o_ce)
		f_waiting_for_first_frame <= 1'b0;

	always @(*)
	if ((f_waiting_for_first_frame)&&(o_ce))
		assert(o_frame);
	always @(*)
	if (f_phase == 0)
		assert((top_of_block)||(first_block));

	always @(*)
	if (f_waiting_for_first_frame)
	begin
		if (f_phase[LGNFFT-1:0] > 3)
		begin
			assert((first_block)&&(!o_frame));
			assert({o_ce, p_ce, d_ce} == 0);
			assert(frame_count == 0);
		end else if (f_phase == 0)
		begin
			assert((!o_frame)&&({o_ce, p_ce, d_ce} == 0));
			assert(frame_count == 0);
		end else if ((f_phase > 0)&&(!first_block))
			assert(|{o_ce, p_ce, d_ce });
		else if ((frame_count != 0)||(first_block))
			assert(!o_frame);
		else
			assert(o_frame);
	end

	always @(posedge i_clk)
	if ((f_past_valid)&&(!$past(i_reset)))
		assert(o_ce || $stable(o_sample));

	/*
	always @(*)
		assert(m_ce == (f_phase == 0));
	always @(*)
		assert(d_ce == (f_phase == 1));
	always @(*)
		assert(p_ce == (f_phase == 2));
	always @(*)
		assert(o_ce == (f_phase == 3));
	always @(*)
		assert(o_frame == (f_phase == 3));
	*/

	always @(*)
		assert(didx[LGNFFT-2:0] == tidx[LGNFFT-2:0]);

	always @(posedge i_clk)
	if ((f_past_valid)&&(!$past(i_reset))
			&&($past(i_ce))&&($past(top_of_block)))
		assert(tidx == 1);


	wire	[LGNFFT:0]	f_phase_plus_one;
	always @(*)
		f_phase_plus_one = f_phase + 1;
	/*
	end else if ($past(i_ce))
	begin
// 1 0		// 0	0		Top of block	0
// 2 1		// 1	1	ALT			0
// 3 1		// 2	1				1
// 4 2		// 3	2	ALT			1
// 5 2		// 4	2				2
// 6 3		// 5	3	ALT			2
// 7 3		// 6	3				3
// 0 4		// 7	4	ALT			3
// 1		// 4	4		Top of block	0
// 2		// 5	5				0
//		// 6	5	ALT			1
//		// 7	6				1
//		// 0	6	ALT			2
//		// 1	7				2
//		// 2	7	ALT
//		// 3	0
//		// 0	0	ALT	Top of block
//		// 1	1
//		// 2	1
//		// 3	2
//		// 4	2
	*/
	always @(*)
	       assert(f_phase_plus_one[LGNFFT:1] == dwidx[LGNFFT-1:0]);
	always @(*)
		assert(f_phase[0] == didx[0]);
	always @(*)
	if (f_phase[LGNFFT])
	begin
		assert(f_phase[LGNFFT-1:0] == {!didx[LGNFFT-1],didx[LGNFFT-2:0]});
		assert((dwidx[LGNFFT-1]==1)
				||(dwidx[LGNFFT-2:0]==0));
	end else begin
		assert((dwidx[LGNFFT-1]==0)
			||((dwidx[LGNFFT-1])&&(dwidx[LGNFFT-2:0]==0)&&(&f_phase[LGNFFT-1:0])));
		assert(f_phase[LGNFFT-1:0] == didx[LGNFFT-1:0]);
	end

	always @(*)
		assert(f_phase[LGNFFT-1:0] == tidx[LGNFFT-1:0]);
	always @(*)
		assert(f_phase[LGNFFT-2:0] == didx[LGNFFT-2:0]);
	wire	[LGNFFT-1:0]	f_diff_idx;
	assign	f_diff_idx = didx - dwidx;
	always @(*)
	if (!first_block)
		assert(f_diff_idx < { 1'b1, {(LGNFFT-1){1'b0}} });


	always @(*)
		if (top_of_block)
			assert(!first_block);

	reg	f_first_block;
	initial	f_first_block = 1'b1;
	always @(posedge i_clk)
	if (i_ce)
		f_first_block <= first_block;

	always @(*)
	if ((top_of_block)&&(i_ce)&&(!f_first_block))
		assert(didx[LGNFFT-2:0] == dwidx[LGNFFT-2:0]);

	always @(*)
	if ((!f_first_block)&&(!first_block)&&(dwidx[LGNFFT-2:0]==0)
			&&(didx[LGNFFT-1:0]==0))
		assert(top_of_block);

////////////////////////////////////////////////////////////////////////////////
//
//  Abstract Multiply
//
////////////////////////////////////////////////////////////////////////////////
	// Gin up a really quick abstract multiply for formal testing
	// only.  always @(posedge i_clk)
	(* anyconst *) signed reg [IW+TW-1:0] pre_product;
	always @(posedge i_clk)
	if (data == 0)
		assume(pre_product == 0);
	always @(posedge i_clk)
	if (tap == 0)
		assume(pre_product == 0);
	always @(posedge i_clk)
	if (data == 1)
		assume(pre_product == tap);
	always @(posedge i_clk)
	if (tap == 1)
		assume(pre_product == data);
	always @(posedge i_clk)
	if ((i_ce)||(i_alt_ce))
		product <= pre_product;

////////////////////////////////////////////////////////////////////////////////
//
//  Arbitrary memory test
//
////////////////////////////////////////////////////////////////////////////////
	(* anyconst *)		reg	[LGNFFT-1:0]	f_addr;
			signed	reg	[TW-1:0]	f_tap;
			signed	reg	[IW-1:0]	f_value, f_tap;
				reg			f_this_dce, f_this_pce,
							f_this_oce, f_this_tap;

	initial	assume(f_tap == cmem[f_addr]);
	always @(*)
		assert(f_tap == cmem[f_addr]);
	always @(posedge i_clk)
	if ((i_tap_wr)&&(f_addr == tapwidx))
		f_tap <= i_tap;

	initial	f_value = 0;
	initial	assume(dmem[f_addr] == f_value);
	always @(*)
		assert(f_value == dmem[f_addr]);
	always @(posedge i_clk)
	if ((i_ce)&&(dwidx == f_addr))
		f_value <= i_sample;
	initial	{ f_this_oce, f_this_pce, f_this_dce } = 3'h0;
	always @(posedge i_clk)
	if ((i_reset)||(i_tap_wr))
		{ f_this_oce, f_this_pce, f_this_dce } <= 3'h0;
	else 
		{ f_this_oce, f_this_pce, f_this_dce }
			<= { f_this_pce, f_this_dce,
				(((i_ce)||(i_alt_ce))
					&&(f_past_valid)&&(f_addr == didx)) };
	initial f_this_tap = 0;
	always @(posedge i_clk)
	if (i_reset)
		f_this_tap <= 0;
	else if ((i_ce)||(i_alt_ce))
		f_this_tap <= (f_past_valid)&&(f_addr == tidx);
	else
		f_this_tap <= 0;


	always @(posedge i_clk)
	if (f_this_tap)
		assert(tap == f_tap);

	always @(posedge i_clk)
	if ((f_past_valid)&&(f_this_dce))
		assert(data == $past(f_value));

	reg	signed	[IW-1:0]	f_past_data;
	reg	signed	[TW-1:0]	f_past_tap;

	always @(posedge i_clk)
	begin
		f_past_data <= data;
		f_past_tap <= tap;
	end

	always @(posedge i_clk)
	if ((f_past_valid)&&(f_this_pce))
	begin
		if (f_past_tap == 0)
			assert(product == 0);
		if (f_past_data == 0)
			assert(product == 0);
		if (f_past_tap == 1)
			assert(product=={{(TW){f_past_data[IW-1]}},f_past_data});
		if (f_past_data == 1)
			assert(product=={{{IW}{f_past_tap[TW-1]}},f_past_tap});
	end


////////////////////////////////////////////////////////////////////////////////
//
//  Cover tests
//
////////////////////////////////////////////////////////////////////////////////
	reg	f_second_frame;
	initial	f_second_frame = 1'b0;
	always @(posedge i_clk)
	if ((o_ce)&&(o_frame))
		f_second_frame <= 1'b1;

	always @(posedge i_clk)
		cover((o_ce)&&(o_frame)&&(f_second_frame));
`endif
endmodule
