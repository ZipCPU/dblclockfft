////////////////////////////////////////////////////////////////////////////////
//
// Filename: 	fftlib.h
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
#ifndef	FFTLIB_H
#define	FFTLIB_H

#define	USE_OLD_MULTIPLY	false

extern	int	lgval(int vl);
extern	int	nextlg(int vl);
extern	int	bflydelay(int nbits, int xtra);
extern	int	lgdelay(int nbits, int xtra);
extern	void	gen_coeffs(FILE *cmem, int stage, int cbits,
			int nwide, int offset, bool inv);
extern	std::string	gen_coeff_fname(const char *coredir,
			int stage, int nwide, int offset, bool inv);
extern	FILE	*gen_coeff_open(const char *fname);
extern	void	gen_coeff_file(const char *coredir, const char *fname,
			int stage, int cbits, int nwide, int offset, bool inv);

#endif	// FFTLIB_H
