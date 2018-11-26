////////////////////////////////////////////////////////////////////////////////
//
// Filename: 	softmpy.h
//
// Project:	A General Purpose Pipelined FFT Implementation
//
// Purpose:	If the chip doesn't have any hardware multiplies, you'll need
// 		a soft-multiply implementation.  This provides that
// 	implementation.
//
// Creator:	Dan Gisselquist, Ph.D.
//		Gisselquist Technology, LLC
//
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2015-2018, Gisselquist Technology, LLC
//
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
#ifndef	SOFTMPY_H
#define	SOFTMPY_H

extern	void	build_multiply(const char *fname);
extern	void	build_bimpy(const char *fname);
extern	void	build_longbimpy(const char *fname);

#endif	// SOFTMPY_H
