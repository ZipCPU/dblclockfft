////////////////////////////////////////////////////////////////////////////////
//
// Filename:	windowfn.v
// {{{
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
// }}}
// Copyright (C) 2018-2021, Gisselquist Technology, LLC
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
//
// License:	LGPL, v3, as defined and found on www.gnu.org,
//		http://www.gnu.org/licenses/lgpl.html
//
//
////////////////////////////////////////////////////////////////////////////////
//
//
`default_nettype	none
// }}}
module	windowfn #(
		// {{{
		parameter		IW=16, OW=16, TW=16, LGNFFT = 4,
		parameter	[0:0]	OPT_FIXED_TAPS = 1'b0,
		parameter		INITIAL_COEFFS = "",
		//
		// OPT_TLAST_FRAME
		//
		// This core can be a challenge to use with an AXI stream
		// interface, simply because the o_frame output normally
		// indicates the *first* value of any frame, not the last.
		// Set OPT_TLAST_FRAME high to adjust this behavior and make
		// o_frame a reference to the last data sample in any FFT frame.
		parameter [0:0]	OPT_TLAST_FRAME = 1'b0,
		// AW : Accumulator width (number of bits)
		localparam	AW=IW+TW
		// }}}
	) (
		// {{{
		input	wire			i_clk, i_reset,
		//
		input	wire			i_tap_wr,
		input	wire	[(TW-1):0]	i_tap,
		//
		input	wire			i_ce,
		input	wire	[(IW-1):0]	i_sample,
		input	wire			i_alt_ce,
		//
		output	reg			o_frame, o_ce,
		output	reg	[(OW-1):0]	o_sample
		// }}}
	);

	// Register declarations
	// {{{
	reg	[(TW-1):0]	cmem	[0:(1<<LGNFFT)-1];
	reg	[(IW-1):0]	dmem	[0:(1<<LGNFFT)-1];

	reg		[LGNFFT-1:0]	dwidx, didx; // Data indices: write&read
	wire		[LGNFFT-1:0]	tidx;	// Coefficient index
	reg				r_tidx;	// Coefficient index, MSB
	reg				top_of_block, first_block;
	reg		[1:0]		frame_count;
	reg				p_ce, d_ce;
	reg	signed	[IW-1:0]	data;
	reg	signed	[TW-1:0]	tap;

	reg	signed	[IW+TW-1:0]	product;
	// }}}
	////////////////////////////////////////////////////////////////////////
	//
	// Load the coefficients
	// {{{
	wire	[LGNFFT-1:0]	tapwidx;
	generate if (OPT_FIXED_TAPS)
	begin : SET_FIXED_TAPS
		// {{{
		initial $readmemh(INITIAL_COEFFS, cmem);

		assign	tapwidx = 0;
		// Make Verilators -Wall happy
		// Verilator lint_off UNUSED
		wire	[TW:0]	ignored_inputs;
		assign	ignored_inputs = { i_tap_wr, i_tap };
		// Verilator lint_on  UNUSED
		// }}}
	end else begin : DYNAMICALLY_SET_TAPS
		// {{{
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
		// }}}
	end endgenerate
	// }}}
	////////////////////////////////////////////////////////////////////////
	//
	// Clock #0: Incoming data
	// {{{
	////////////////////////////////////////////////////////////////////////
	//
	//	Write data to memory, so we can come back to it later to get
	//	the other half of our 50% overlap.
	//

	//
	//
	// Record the incoming data into a local memory
	//
	//

	// Notice how this data writing section is *independent* of the reset,
	// depending only upon new sample data.

	// dwidx
	// {{{
	initial	dwidx = 0;
	always @(posedge i_clk)
	if (i_reset)
		dwidx <= 0;
	else if (i_ce)
		dwidx <= dwidx + 1'b1;
	// }}}

	// Write to dmem
	// {{{
	always @(posedge i_clk)
	if (i_ce)
		dmem[dwidx] <= i_sample;
	// }}}

	// first_block
	// {{{
	// first_block is our attempt to make certain that the entire first
	// FFT's worth of data is repressed.  This is the data that wouldn't
	// be valid, due to the data in memory not being valid.  In the case
	// of a first_block, we do everything like we otherwise would--only
	// we don't output either o_ce or o_frame--hence none of the following
	// processing should operate on our outputs until we have a full and
	// valid frame of data.
	//
	initial	first_block = 1'b1;
	always @(posedge i_clk)
	if (i_reset)
		first_block <= 1'b1;
	else if ((i_alt_ce)&&(&tidx)&&(dwidx==0))
		first_block <= 1'b0;
	// }}}
	// }}}
	////////////////////////////////////////////////////////////////////////
	//
	// Clock #0: Sample arrives
	// {{{
	////////////////////////////////////////////////////////////////////////
	//
	//	Indicator signals:	i_ce, or i_alt_ce
	//
	//	Valid signals:
	//		This is the state machine stage.  Signals here are
	//		valid from clock to clock
	//
	//		first_block
	//


	// top_of_block
	// {{{
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
	// }}}

	// didx, tidx
	// {{{
	// Data and coefficient memory indices.
	//
	// The data index goes from 0:(1<<LGNFFT)-1 for the first FFT, and then
	// (1<<(LGNFFT-1)):(1<<LGNFFT)-1,0:(1<<(LGNFFT-1))-1 for the overlapped
	// FFT.  During these same two runs, the tap (coefficient) index goes
	// from 0:(1<<LGNFFT)-1 each time, so the bottom (LGNFFT-1) bits can
	// be shared between the two pointers.
	//
	// Note that each of these pointers is set the clock before the data
	// actually arrives (as it should be)
	//
	initial	didx = 0;
	always @(posedge i_clk)
	if (i_reset)
		didx <= 0;
	else if ((i_alt_ce)&&(dwidx[LGNFFT-2:0]==0))
	begin
		// Restart on the first point of the next FFT
		didx[LGNFFT-2:0] <= 0;
		// Maintain the top bit, so as to keep
		// the overlap working properly
		didx[LGNFFT-1] <= dwidx[LGNFFT-1];
	end else if ((i_ce)||(i_alt_ce))
		// Process the next point in this FFT
		didx <= didx + 1'b1;

	assign tidx = { r_tidx, didx[LGNFFT-2:0] };

	initial	r_tidx = 0;
	always @(posedge i_clk)
	if (i_reset)
		r_tidx <= 0;
	else if ((i_alt_ce)&&(dwidx[LGNFFT-2:0]==0))
	begin
		// Restart the counter for the first point
		// of the next FFT.
		r_tidx <= 0;
	end else if ((i_ce)||(i_alt_ce))
	begin
		// Process the next point in the window function
		//
		// Here we implement the carry only
		if (&tidx[LGNFFT-2:0])
			r_tidx <= !r_tidx;
	end
	// }}}

	// frame_count
	// {{{
	//
	// frame_count is based off of the first index at the top of the
	// block.  It's used to make certain that the o_frame output is
	// properly aligned with the first valid clock of the output.
	//
	initial	frame_count = 0;
	always @(posedge i_clk)
	if (i_reset)
		frame_count <= 0;
	else if (OPT_TLAST_FRAME && i_alt_ce && (&tidx)&& !first_block)
		frame_count <= 3;
	else if (!OPT_TLAST_FRAME && (i_ce)&&(top_of_block)&&(!first_block))
		frame_count <= 3;
	else if (frame_count != 0)
		frame_count <= frame_count - 1'b1;
	// }}}

	// p_ce, d_ce
	// {{{
	// Following any initial i_ce, ...
	// 	d_ce: The data (and coefficient), read from memory,will be valid
	// 	p_ce: The produc of data and coefficient is valid
	//
	initial	{ p_ce, d_ce } = 2'h0;
	always @(posedge i_clk)
	if (i_reset)
		{ p_ce, d_ce } <= 2'h0;
	else
		{ p_ce, d_ce } <= { d_ce, (!first_block)&&((i_ce)||(i_alt_ce))};
	// }}}
	// }}}
	////////////////////////////////////////////////////////////////////////
	//
	// Clock #1
	// {{{
	////////////////////////////////////////////////////////////////////////
	//
	//	Indicator signals: (i_ce || i_alt_ce)
	//	Valid signals:	didx, tidx, top_of_block
	//

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
	// }}}
	////////////////////////////////////////////////////////////////////////
	//
	// Clock #2
	// {{{
	////////////////////////////////////////////////////////////////////////
	//
	//	Indicator signals: d_ce
	//	Valid signals:	data, tap
	//

	//
	// Multiply the two values together
`ifdef	FORMAL
	// We'll implement an abstract multiply below--just to make sure the
	// timing is right.
`else
	always @(posedge i_clk)
		product <= data * tap;
`endif
	// }}}
	////////////////////////////////////////////////////////////////////////
	//
	// Clock #3
	// {{{
	////////////////////////////////////////////////////////////////////////
	//
	//	Indicator signals:	p_ce
	//	Valid signals:		product, (frame_count)
	//

	// o_ce
	// {{{
	// Output CE.  This value will be true exactly two clocks after any
	//	(i_ce || i_alt_ce) input, indicating that we have a new valid
	//	output sample.
	//
	initial	o_ce = 0;
	always @(posedge i_clk)
	if (i_reset)
		o_ce <= 0;
	else
		o_ce <= p_ce;
	// }}}

	// o_frame
	// {{{
	// o_frame is the indicator that this is the first clock cycle of a
	// new FFT frame.  It's *like* TLAST in its function, but acts on the
	// first cycle, rather than the last.  If you want to use this with
	// an AXI stream, you can set OPT_TLAST_FRAME above, and then o_frame
	// will be set on the last o_ce clock cycle of any outgoing frame.
	//
	initial	o_frame = 1'b0;
	always @(posedge i_clk)
	if (i_reset)
		o_frame <= 1'b0;
	else if (frame_count == 2)
		o_frame <= 1'b1;
	else
		o_frame <= 1'b0;
	// }}}

	generate if (OW == AW)
	begin : BIT_ADJUSTMENT_NONE
		// {{{
		initial	o_sample = 0;
		always @(posedge i_clk)
		if (i_reset)
			o_sample <= 0;
		else if (p_ce)
			o_sample <= product;
		// }}}
	end else if (OW < AW)
	begin : BIT_ADJUSTMENT_ROUNDING
		// {{{
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
		// }}}
	end else // if (OW > AW)
	begin : BIT_ADJUSTMENT_EXTENDING
		// {{{
		always @(posedge i_clk)
		if (i_reset)
			o_sample <= 0;
		else if (p_ce)
			o_sample <= { product, {(OW-AW){1'b0}} };
		// }}}
	end endgenerate
	// }}}
	////////////////////////////////////////////////////////////////////////
	//
	// Clock #4: All done
	// {{{
	////////////////////////////////////////////////////////////////////////
	//
	//	Indicator signals:	o_ce
	//	Valid signals:		o_sample
	//

	// }}}
	// Make Verilator happy
	// {{{
	// verilator lint_off UNUSED
	wire	[LGNFFT-1:0]	unused;
	assign	unused = tapwidx;
	// verilator lint_on  UNUSED
	// }}}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
// Formal properties
// {{{
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
`ifdef	FORMAL
	// Declarations
	// {{{
`ifdef	WINDOWFN
`define	ASSUME	assume
`else
`define	ASSUME	assert
`endif
	reg			f_past_valid;
	reg			f_waiting_for_first_frame;
	reg	[LGNFFT-1:0]	f_tidx;
	reg	[LGNFFT:0]	f_phase_plus_one;
	reg	[LGNFFT:0]	f_phase;
	wire	[LGNFFT-1:0]	f_diff_idx;
	reg	f_first_block;
	(* anyconst *)	reg		[LGNFFT-1:0]	f_addr;
			reg	signed	[TW-1:0]	f_tap;
			reg	signed	[IW-1:0]	f_value;
			reg				f_this_dce, f_this_pce,
							f_this_oce, f_this_tap;
	reg	signed	[IW-1:0]	f_past_data;
	reg	signed	[TW-1:0]	f_past_tap;

	initial	f_past_valid = 1'b0;
	always @(posedge i_clk)
		f_past_valid <= 1'b1;
	// }}}
	////////////////////////////////////////////////////////////////////////
	//
	// f_phase -- Keeping track of where we are at in this operation
	// {{{
	////////////////////////////////////////////////////////////////////////
	//
	//

	// f_phase is a one up counter starting from reset, to which
	// we will attach all of our other counters via assertions
	//

	initial	f_phase = 0;
	always @(posedge i_clk)
	if (i_reset)
	begin
		f_phase <= 0;
	end else if ((i_ce)&&(top_of_block))
		f_phase[LGNFFT-1:0] <= 1;
	else if ((i_ce)||(i_alt_ce))
		f_phase <= f_phase + 1;
	// }}}
	////////////////////////////////////////////////////////////////////////
	//
	// Assumptions about the input
	// {{{
	////////////////////////////////////////////////////////////////////////
	//
	//
`ifdef	VERIFIC
	restrict property (@(posedge i_clk)
		i_reset && !i_ce
		|=> $stable(i_sample));
`else
	always @(posedge i_clk)
	if ((f_past_valid)&&(!$past(i_reset))&&(!$past(i_ce)))
		restrict($stable(i_sample));
`endif

	//
	// Only accept data after the coefficients have been written
	//
	always @(*)
	if (tapwidx != 0)
		`ASSUME(!i_ce);

	always @(*)
	if (i_tap_wr)
		`ASSUME((!i_ce)&&(!i_alt_ce));

	//
	// Assume the user isn't writing to our taps at all
	//
	always @(*)
		`ASSUME(!i_tap_wr);

	//
	// Insist that the various CE's alternate: i_ce, i_alt_ce, i_ce, etc.
	//
	always @(posedge i_clk)
	if (f_phase[0])
		`ASSUME(!i_ce);

	always @(posedge i_clk)
	if (!f_phase[0])
		`ASSUME(!i_alt_ce);

	always @(*)
		assert(!i_ce || !i_alt_ce);
	// }}}
	////////////////////////////////////////////////////////////////////////
	//
	// The top_of_block signal
	// {{{
	////////////////////////////////////////////////////////////////////////
	//
	//

	//
	// Our top_of_block signal should only be high when we are ready to
	// accept a new FFT.  That is, when the data write index is zero,
	// otherwise top_of_block should be zero
	//
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
	// }}}
	////////////////////////////////////////////////////////////////////////
	//
	// Assertions about our outputs
	// {{{
	////////////////////////////////////////////////////////////////////////
	//
	//

	always @(*)
	if (first_block)
		assert(!o_ce);

	initial	f_waiting_for_first_frame = 1;
	always @(posedge i_clk)
	if ((i_reset)||(first_block))
		f_waiting_for_first_frame <= 1'b1;
	else if (o_ce)
		f_waiting_for_first_frame <= 1'b0;

	always @(*)
	if (!o_ce)
		assert(!o_frame);
	else if (!OPT_TLAST_FRAME && f_waiting_for_first_frame)
		assert(o_frame);

	always @(posedge i_clk)
	if (OPT_TLAST_FRAME && f_past_valid && $past(f_past_valid) && o_ce)
		assert(o_frame == $past(top_of_block,2));
	// else if (OPT_TLAST_FRAME

	always @(*)
	if (OPT_TLAST_FRAME)
	begin
		if (frame_count == 3)
			assert(top_of_block);
		else if (frame_count != 0)
			assert(top_of_block || d_ce || p_ce);
		if (frame_count <= 2 && top_of_block)
			assert(!d_ce);
		if (frame_count <= 1 && top_of_block)
			assert(!d_ce && !p_ce);
		if (frame_count == 0 && top_of_block)
			assert(!d_ce && !p_ce && !o_ce);
	end

	always @(*)
	if (f_phase[LGNFFT-1:0] == 0)
		assert((top_of_block)||(first_block));
	else
		assert(!top_of_block);

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
			// Never get here
			assert(!o_frame);
		else
			assert(o_frame);
	end

	always @(posedge i_clk)
	if ((f_past_valid)&&(!$past(i_reset)))
		assert(o_ce || $stable(o_sample));

	always @(*)
		assert(didx[LGNFFT-2:0] == tidx[LGNFFT-2:0]);

	always @(*)
		f_tidx = f_phase;

	always @(*)
		assert(f_tidx == tidx);

	always @(posedge i_clk)
	if ((f_past_valid)&&(!$past(i_reset))
			&&($past(i_ce))&&($past(top_of_block)))
		assert(tidx == 1);


	always @(*)
		f_phase_plus_one = f_phase + 1;

	//
	// Assert the relationship between dwidx and f_phase
	//
	always @(*)
	       assert(f_phase_plus_one[LGNFFT:1] == dwidx[LGNFFT-1:0]);

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

	assign	f_diff_idx = didx - dwidx;
	always @(*)
	if (!first_block)
		assert(f_diff_idx < { 1'b1, {(LGNFFT-1){1'b0}} });


	always @(*)
	if (top_of_block)
		assert(!first_block);

	initial	f_first_block = 1'b1;
	always @(posedge i_clk)
	if (i_reset)
		f_first_block <= 1'b1;
	else if (i_ce)
		f_first_block <= first_block;

	always @(*)
	if (first_block)
		assert(f_first_block);

	always @(*)
	if ((top_of_block)&&(i_ce))
		assert(didx[LGNFFT-2:0] == dwidx[LGNFFT-2:0]);

	always @(*)
	if ((!f_first_block)&&(!first_block)&&(dwidx[LGNFFT-2:0]==0)
			&&(didx[LGNFFT-1:0]==0))
		assert(top_of_block);
	// }}}
	////////////////////////////////////////////////////////////////////////
	//
	//  Abstract Multiply
	// {{{
	////////////////////////////////////////////////////////////////////////
	//
	//

	// Gin up a really quick abstract multiply for formal testing
	// only.  always @(posedge i_clk)
	(* anyseq *) reg signed [IW+TW-1:0] f_pre_product;
	always @(posedge i_clk)
	if (data == 0)
		assume(f_pre_product == 0);
	always @(posedge i_clk)
	if (tap == 0)
		assume(f_pre_product == 0);
	always @(posedge i_clk)
	if (data == 1)
		assume(f_pre_product == tap);
	always @(posedge i_clk)
	if (tap == 1)
		assume(f_pre_product == data);
	always @(posedge i_clk)
	if ($stable(data) && $stable(tap))
		assume($stable(f_pre_product));
	always @(posedge i_clk)
		product <= f_pre_product;
	// }}}
	////////////////////////////////////////////////////////////////////////
	//
	//  Arbitrary memory test
	// {{{
	////////////////////////////////////////////////////////////////////////

`ifdef	VERIFIC
	always @(*)
	if (!f_past_valid)
	begin
		assume(f_tap == cmem[f_addr]);
		assume(dmem[f_addr] == f_value);
	end
`else
	initial	assume(f_tap == cmem[f_addr]);
	initial	assume(dmem[f_addr] == f_value);
`endif

	always @(*)
		assert(f_tap == cmem[f_addr]);
	always @(posedge i_clk)
	if ((i_tap_wr)&&(f_addr == tapwidx))
		f_tap <= i_tap;

	initial	f_value = 0;
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
			<= { f_this_pce, f_this_dce, (((i_ce)||(i_alt_ce))
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

	always @(posedge i_clk)
	begin
		f_past_data <= data;
		f_past_tap  <= tap;
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
			assert(product ==
				{ {{IW}{f_past_tap[TW-1]}},f_past_tap});
	end
	// }}}
	////////////////////////////////////////////////////////////////////////
	//
	//  Cover tests
	// {{{
	////////////////////////////////////////////////////////////////////////
	reg	cvr_second_frame;

	initial	cvr_second_frame = 1'b0;
	always @(posedge i_clk)
	if ((o_ce)&&(o_frame))
		cvr_second_frame <= 1'b1;

	always @(posedge i_clk)
		cover((o_ce)&&(o_frame)&&(cvr_second_frame));
	// }}}
`endif
// }}}
endmodule
