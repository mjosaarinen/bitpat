//	bitpat.h
//	2021-01-31	Markku-Juhani O. Saarinen <mjos@pqshield.com>
//	Copyright (c) 2021, PQShield Ltd.  All rights reserved.

//	===	main header file

#ifndef _BITPAT_H_
#define _BITPAT_H_

#include <stdint.h>
#include <stddef.h>

//	Gaussian tailcut, in standard deviation "sigmas".
//	From CDF, epsilon ~= erfc( tau / sqrt(2) ); 15 is good for IEEE Float.

#ifndef TAILCUT_TAU
#define TAILCUT_TAU 15.0
#endif

//	===	Compute min-entropy -log2(max p_z)
//	f	fequency [0,1] (peak)
//	d	cutoff (0.5 = no bias)
//	s2	jitter variance
//	n	Zn -- the bit sample size
//	m	FFT size (must be power of 2)

double jitter_fft(double f, double d, double s2, size_t n, size_t m);

//	===	Simulate min-entropy -log2(max p_z)
//	f	fequency [0,1] (peak)
//	d	cutoff  [0,1] (0.5 = no bias)
//	s2	jitter variance
//	n	Zn -- the bit sample size
//	m	number of iterations

double jitter_sim(double f, double d, double s2, size_t n, size_t m);

//	===	Simulate byte strings into *zv. "zvlen" is length in bytes.

void zbytes(uint8_t *zv, size_t zvlen, double f, double d, double s2);

//	===	f_s	step function.
//	vs	result, scaled as vs[i] = f_s(i/m)
//	m	length
//	f	fequency [0,1] (peak)
//	s2	jitter variance

double vec_fs(double *vs, size_t m, double f, double s2);

//	===	Chop vs[m] to "bit", cutoff at "d". return vr and area sum.

double vec_chop(double *vr, const double *vs, size_t m, double d, int bit);

#endif	//	_BITPAT_H_
