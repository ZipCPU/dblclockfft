////////////////////////////////////////////////////////////////////////////////
//
// Filename: 	rounding.h
// {{{
// Project:	A General Purpose Pipelined FFT Implementation
//
// Purpose:	To create one of a series of modules to handle dropping bits
// 		within the FFT implementation.
//
// Creator:	Dan Gisselquist, Ph.D.
//		Gisselquist Technology, LLC
//
////////////////////////////////////////////////////////////////////////////////
// }}}
// Copyright (C) 2015-2021, Gisselquist Technology, LLC
// {{{
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
// }}}
// License:	GPL, v3, as defined and found on www.gnu.org,
// {{{
//		http://www.gnu.org/licenses/gpl.html
//
//
////////////////////////////////////////////////////////////////////////////////
//
// }}}
#ifndef	ROUNDING_H
#define	ROUNDING_H

typedef	enum {
	RND_TRUNCATE, RND_FROMZERO, RND_HALFUP, RND_CONVERGENT
} ROUND_T;


extern	void	build_truncator(const char *fname);
extern	void	build_roundhalfup(const char *fname);
extern	void	build_roundfromzero(const char *fname);
extern	void	build_convround(const char *fname);

#endif
