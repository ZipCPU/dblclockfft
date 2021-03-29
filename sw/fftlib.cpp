////////////////////////////////////////////////////////////////////////////////
//
// Filename: 	fftlib.cpp
// {{{
// Project:	A General Purpose Pipelined FFT Implementation
//
// Purpose:	
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
#define _CRT_SECURE_NO_WARNINGS   //  ms vs 2012 doesn't like fopen
#include <stdio.h>
#include <stdlib.h>

#ifdef _MSC_VER //  added for ms vs compatibility

#include <io.h>
#include <direct.h>
#define _USE_MATH_DEFINES

#if _MSC_VER <= 1700

long long llround(double d) {
	if (d<0) return -(long long)(-d+0.5);
	else	return (long long)(d+0.5); }

#endif

#else
// And for G++/Linux environment

#include <unistd.h>	// Defines the R_OK/W_OK/etc. macros
#endif

#include <string.h>
#include <string>
#include <math.h>
// #include <ctype.h>
#include <assert.h>

#include "fftlib.h"


// lgval(vl) -- returns ceil(log_2(val))
// {{{
int	lgval(int vl) {
	int	lg;

	for(lg=1; (1<<lg) < vl; lg++)
		;
	return lg;
}
// }}}

// nextlg(vl) -- returns log_2(ceil(log_2(val)))
// {{{
int	nextlg(int vl) {
	int	r;

	for(r=1; r<vl; r<<=1)
		;
	return r;
}
// }}}

// bflydelay -- Clocks used in the butterfly
// {{{
int	bflydelay(int nbits, int xtra) {
	int	cbits = nbits + xtra;
	int	delay;

	if (USE_OLD_MULTIPLY) {
		if (nbits+1<cbits)
			delay = nbits+4;
		else
			delay = cbits+3;
	} else {
		int	na=nbits+2, nb=cbits+1;
		if (nb<na) {
			int tmp = nb;
			nb = na; na = tmp;
		} delay = ((na)/2+(na&1)+2);
	}
	return delay;
}
// }}}

// lgdelay -- log of thebutterfly delay (i.e. bits needed to hold the value)
// {{{
int	lgdelay(int nbits, int xtra) {
	// The butterfly code needs to compare a valid address, of this
	// many bits, with an address two greater.  This guarantees we
	// have enough bits for that comparison.  We'll also end up with
	// more storage space to look for these values, but without a
	// redesign that's just what we'll deal with.
	return lgval(bflydelay(nbits, xtra)+3);
}
// }}}

// gen_coeffs -- generate twiddle factors
// {{{
void	gen_coeffs(FILE *cmem, int stage, int cbits,
			int nwide, int offset, bool inv) {
	unsigned long	ucbits = (unsigned long)cbits;
	//
	// For an FFT stage of 2^n elements, we need 2^(n-1) butterfly
	// coefficients, sometimes called twiddle factors.  Stage captures the
	// width of the FFT at this point.  If this is a 2x at a time FFT,
	// nwide will be equal to 2, and offset will be one or two.
	//
	// assert(nwide > 0);
	// assert(offset < nwide);
	// assert(stage / nwide >  1);
	// assert(stage % nwide == 0);
	// printf("GEN-COEFFS(): stage =%4d, bits =%2d, nwide = %d, offset = %d, nverse = %d\n", stage, cbits, nwide, offset, inv);
	if (ucbits >= 8*sizeof(long long)) {
		fprintf(stderr, "ERROR: CMEM coefficient precision requested (%d / coefficient) overflows long long data type\n", cbits);
		exit(EXIT_FAILURE);
	}

	fprintf(cmem, "// Coefficient memory\n");
	fprintf(cmem, "// ----------------------------------------------\n");
	fprintf(cmem, "//   Stage:               %3d\n", stage);
	fprintf(cmem, "//   Bits per coefficient:%3d\n", cbits);
	fprintf(cmem, "//   NWide:               %3d%s\n", nwide,
			(nwide == 2) ? " Double-wide FFT, two samples per clock"
			:" Normal FFT, one sample per cycle");
	fprintf(cmem, "//   Offset:              %3d%s\n", offset,
			(nwide < 2) ? " (Ignored)"
			: (offset == 0) ? " (Even coefficients)"
			: " (Odd coefficients");
	fprintf(cmem, "//   Inv:               %s\n",
			(inv) ? " True -- FFT is inverted"
			: "False -- This is a forward FFT");
	fprintf(cmem, "//\n//\n");
	fprintf(cmem, "// Each line contains a coefficient.  The real portion\n");
	fprintf(cmem, "// of the coefficient is in the upper %d bits, whereas\n", cbits);
	fprintf(cmem, "// the lower %d bits contain the imaginary portion\n", cbits);
	fprintf(cmem, "//\n//\n");
	int	ncoeffs = stage/nwide/2;
	for(int i=0; i<ncoeffs; i++) {
		int k = nwide*i+offset;
		double	W = ((inv)?1:-1)*2.0*M_PI*k/(double)(stage);
		double	c, s;
		typedef	long long		ll;
		typedef	unsigned long long	ull;
		ll	ic, is;
		ull	vl, uic, uis;

		c = cos(W); s = sin(W);
		ic = (long long)llround((1ll<<(cbits-2)) * c);
		is = (long long)llround((1ll<<(cbits-2)) * s);

		if (ucbits != 8*sizeof(long long)) {
			uic = ic & (~(-1ll << ucbits));
			uis = is & (~(-1ll << ucbits));
		} else {
			uic = ic;
			uis = is;
		}

		if (2*ucbits <= 8 * sizeof(long long)) {
			vl = uic << ucbits;
			vl |= uis;
			fprintf(cmem, "%0*llx\n", ((cbits*2+3)/4), vl);
		} else if (ucbits == 8 * sizeof(long long)) {
			// Each coefficient perfectly fits into a long long.
			fprintf(cmem, "%0*llx\n", ((cbits+3)/4), uic);
			fprintf(cmem, "%0*llx\n", ((cbits+3)/4), uis);
		} else {
			unsigned long long	hi, lo;
			unsigned sright;

			// The "lo" number must always have 8*sizeof(long long)
			// bits in it.  That means that the "hi" number must
			// have (cbits * 2) - (8*sizeof(long long)) valid bits
			// in it.  Since we're in this branch of the "if", we
			// know that there's a non-zero number of bits in hi.
			// Further, we know that this non-zero number is more
			// than zero, and less than 8*sizeof(long long).
			// Therefore, hi must contain a *part* of ic and not
			// all of it, whereas lo must contain the last part of
			// ic.

			vl = (ic & (~(-1ll << (ucbits))));
			sright = (8*sizeof(long long)) - ucbits;
			hi = vl >> sright;
			fprintf(cmem, "%0*llx",
				(int)(((2*cbits-8*sizeof(long long))+3)/4), hi);

			lo = uis;
			vl <<= cbits;
			lo  |= vl;

			fprintf(cmem, "%0*llx\n", (int)(8*sizeof(long long)/4), lo);
		}
		//
	} fclose(cmem);
}
// }}}

// gen_coef_fname -- Generates the hex file name for the twiddle factors
// {{{
std::string	gen_coeff_fname(const char *coredir,
			int stage, int nwide, int offset, bool inv) {
	std::string	result;
	char	*memfile;

	assert((nwide == 1)||(nwide == 2));

	memfile = new char[strlen(coredir)+3+10+strlen(".hex")+64];
	if (nwide == 2) {
		if (coredir[0] == '\0') {
			sprintf(memfile, "%scmem_%c%d.hex",
				(inv)?"i":"", (offset==1)?'o':'e', stage*nwide);
		} else {
			sprintf(memfile, "%s/%scmem_%c%d.hex",
				coredir, (inv)?"i":"",
				(offset==1)?'o':'e', stage*nwide);
		}
	} else if (coredir[0] == '\0') // if (nwide == 1)
		sprintf(memfile, "%scmem_%d.hex",
			(inv)?"i":"", stage);
	else
		sprintf(memfile, "%s/%scmem_%d.hex",
			coredir, (inv)?"i":"", stage);

	result = std::string(memfile);
	delete[] memfile;
	return	result;
}
// }}}

// gen_coeff_open
// {{{
FILE	*gen_coeff_open(const char *fname) {
	FILE	*cmem;

	cmem = fopen(fname, "w");
	if (NULL == cmem) {
		fprintf(stderr, "Could not open FFT coefficient file "
				"\'%s\' for writing\n", fname);
		perror("Err from O/S:");
		exit(EXIT_FAILURE);
	}

	return cmem;
}
// }}}

// gen_coeff_file
// {{{
void	gen_coeff_file(const char *coredir, const char *fname,
			int stage, int cbits, int nwide, int offset, bool inv) {
	std::string	fstr;
	FILE	*cmem;

	fstr= gen_coeff_fname(coredir, stage, nwide, offset, inv);
	cmem = gen_coeff_open(fstr.c_str());
	gen_coeffs(cmem, stage,  cbits, nwide, offset, inv);
}
// }}}
