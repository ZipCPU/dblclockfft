////////////////////////////////////////////////////////////////////////////////
//
// Filename: 	defaults.h
//
// Project:	A General Purpose Pipelined FFT Implementation
//
// Purpose:
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
#ifndef	DEFAULTS_H
#define	DEFAULTS_H

#define	DEF_NBITSIN	16
#define	DEF_COREDIR	"fft-core"
#define	DEF_XTRACBITS	4
#define	DEF_NMPY	0
#define	DEF_XTRAPBITS	0
#define	USE_OLD_MULTIPLY	false

// To coordinate testing, it helps to have some defines in our header file that
// are common with the default parameters found within the various subroutines.
// We'll define those common parameters here.  These values, however, have no
// effect on anything other than bench testing.  They do, though, allow us to
// bench test exact copies of what is going on within the FFT when necessary
// in order to find problems.
// First, parameters for the new multiply based upon the bi-multiply structure
// (2-bits/2-tableau rows at a time).
#define	TST_LONGBIMPY_AW	16
#define	TST_LONGBIMPY_BW	20	// Leave undefined to match AW

//  We also include parameters for the shift add multiply
#define	TST_SHIFTADDMPY_AW	16
#define	TST_SHIFTADDMPY_BW	20	// Leave undefined to match AW

// Now for parameters matching the butterfly
#define	TST_BUTTERFLY_IWIDTH	16
#define	TST_BUTTERFLY_CWIDTH	20
#define	TST_BUTTERFLY_OWIDTH	17

// Now for parameters matching the qtrstage
#define	TST_QTRSTAGE_IWIDTH	16
#define	TST_QTRSTAGE_LGWIDTH	8

// Parameters for the dblstage
#define	TST_DBLSTAGE_IWIDTH	16
#define	TST_DBLSTAGE_SHIFT	0

// Now for parameters matching the dblreverse stage
#define	TST_DBLREVERSE_LGSIZE	5

static	const	bool	formal_property_flag = false;

#endif