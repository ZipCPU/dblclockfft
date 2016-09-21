/////////////////////////////////////////////////////////////////////////////
//
// Filename: 	../bench/cpp/ifftsize.h
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
///////////////////////////////////////////////////////////////////////////
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
///////////////////////////////////////////////////////////////////////////
//
//
#ifndef IFFTHDR_H
#define IFFTHDR_H

#define	IFFT_IWIDTH	22
#define	IFFT_OWIDTH	28
#define	IFFT_LGWIDTH	11
#define	IFFT_SIZE	(1<<IFFT_LGWIDTH)

#define	DBLCLKIFFT


#endif

