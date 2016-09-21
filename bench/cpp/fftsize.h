/////////////////////////////////////////////////////////////////////////////
//
// Filename: 	../bench/cpp/fftsize.h
//
// Project:	A Doubletime Pipelined FFT
//
// Purpose:	This simple header file captures the internal constants
//		within the FFT that were used to build it, for the purpose
//		of making C++ integration (and test bench testing) simpler.  That
//		is, should the FFT change size, this will note that size change
//		and thus any test bench or other C++ program dependent upon
//		either the size of the FFT, the number of bits in or out of
//		it, etc., can pick up the changes in the defines found within
//		this file.
//
// Creator:	Dan Gisselquist, Ph.D.
//		Gisselquist Technology, LLC
//
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2015, Gisselquist Technology, LLC
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
#ifndef FFTHDR_H
#define FFTHDR_H

#define	FFT_IWIDTH	16
#define	FFT_OWIDTH	22
#define	FFT_LGWIDTH	11
#define	FFT_SIZE	(1<<FFT_LGWIDTH)

#define	DBLCLKFFT

// Parameters for testing the longbimpy
#define	TST_LONGBIMPY_AW	16
#define	TST_LONGBIMPY_BW	20

// Parameters for testing the shift add multiply
#define	TST_SHIFTADDMPY_AW	16
#define	TST_SHIFTADDMPY_BW	20

// Parameters for testing the butterfly
#define	TST_BUTTERFLY_IWIDTH	16
#define	TST_BUTTERFLY_CWIDTH	20
#define	TST_BUTTERFLY_OWIDTH	17
#define	TST_BUTTERFLY_MPYDELAY	11

// Parameters for testing the quarter stage
#define	TST_QTRSTAGE_IWIDTH	16
#define	TST_QTRSTAGE_LGWIDTH	8

// Parameters for testing the double stage
#define	TST_DBLSTAGE_IWIDTH	16
#define	TST_DBLSTAGE_SHIFT	0

// Parameters for testing the bit reversal stage
#define	TST_DBLREVERSE_LGSIZE	5


#endif

