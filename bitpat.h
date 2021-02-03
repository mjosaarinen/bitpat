//	bitpat.h
//	2021-01-31	Markku-Juhani O. Saarinen <mjos@pqshield.com>
//	Copyright (c) 2021, PQShield Ltd.  All rights reserved.

//	===	Main header file for jitter bit pattern things.

#ifndef _BITPAT_H_
#define _BITPAT_H_

#include <stdint.h>
#include <stddef.h>

//	Gaussian tailcut, in standard deviation "sigmas".
//	From CDF, epsilon ~= erfc( tau / sqrt(2) ); 15 is good for IEEE Float.

#ifndef TAILCUT_TAU
#define TAILCUT_TAU 15.0
#endif

//	===	Evaluate min-entropy -log2(max p_z) using FFT (pzfft.c)
//	f	fequency [0,1] (peak)
//	d	cutoff  [0,1] (0.5 = no bias)
//	s2	jitter variance
//	n	Zn -- the bit sample size
//	m	FFT size (must be power of 2)
//	v	verbose (0 = print nothing, 1 = distribution to stdout)


double entropy_fft(double f, double d, double s2, size_t n, size_t m, int v);

//	===	Estimate min-entropy -log2(max p_z) using simulation. (simul.c)
//	f	fequency [0,1] (peak)
//	d	cutoff  [0,1] (0.5 = no bias)
//	s2	jitter variance
//	n	Zn -- the bit sample size
//	m	number of iterations
//	v	verbose (0 = print nothing, 1 = distribution to stdout)

double entropy_sim(double f, double d, double s2, size_t n, size_t m, int v);


//	===	Evaluate autocorrelation Ck from a definite integral. (autoc.c)
//	(f, d, s2)	parameters
//	k	distance 

double ck_eval(double f, double d, double s2, int k);


//	===	Estimate autocorrelation vector Ck[] from simulation (autoc.c)
//	*ck	result vector, each entry scaled to -1 <= ck[k] <= 1
//	(f, d, s2)	parameters
//	l	length of ck[] (at most 64)
//	m	number of simulation iterations (at least l)

void ck_sim(double *ck, double f, double d, double s2, size_t l, size_t m);


//	---	helper functions

//	===	computation of s1 integral (autoc.c)
//	(f, d, s2)	parameters
//	x	point where computed

double eval_s1(double f, double d, double s2, double x);


//	===	f_s	step function (pzfft.c)
//	vs	result, scaled as vs[i] = f_s(i/m)
//	m	length
//	f	fequency [0,1] (peak)
//	s2	jitter variance

double vec_fs(double *vs, size_t m, double f, double s2);


//	===	Chop vs[m] to "bit", cutoff at "d". return vr and area sum (pzfft.c)

double vec_chop(double *vr, const double *vs, size_t m, double d, int bit);


//	=== Simulate byte strings into *zv. "zvlen" is length in bytes (simul.c)

void zbytes(uint8_t *zv, size_t zvlen, double f, double d, double s2);


//	---	experimental stuff (testing.c)

//	===	Estimate min-entropy -log2(max p_z) by depth-first search on z
//	(f, d, s2)	parameters
//	n	Zn -- the bit sample size
//	m	FFT size (must be power of 2)
//	how	strategy
//	v	verbose (0 = print nothing, 1 = distribution to stdout)

enum dfs_how {
	DFS_BIG_MASS = 0,
	DFS_FOLLOW_X
};

typedef enum dfs_how dfs_how_t;

double entropy_dfs(double f, double d, double s2, size_t n, size_t m, 
	int v, dfs_how_t how);

#endif	//	_BITPAT_H_

