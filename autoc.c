//	autoc.c
//	2021-01-08	Markku-Juhani O. Saarinen <mjos@pqshield.com>
//	Copyright (c) 2021, PQShield Ltd.  All rights reserved.

//	===	Autocorrelation computations for jitter model.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "bitpat.h"

//	===	computation of s1 integral
//	(f, d, s2)	parameters
//	x	point where computed

double eval_s1(double f, double d, double s2, double x)
{
	double xi, ai, bi;

	double u = sqrt(2.0 * s2);
	double t = 1.0 / u;			//	scale factor
	double tc = ceil(TAILCUT_TAU * sqrt(s2));
	double erfs = 0.0;
	double exps = 0.0;
	
	//	run sum x-f - tc ... x-f + tc
	for (xi = x - f - tc; xi <= x - f + tc; xi += 1.0) {
		ai = t * xi;
		bi = t * (xi - d);
		erfs += ai * erf(ai) - bi * erf(bi);
		exps += exp(-ai*ai) - exp(-bi*bi);
	}
	
	//	1/sqrt(Pi)
	erfs += 0.56418958354775628694807945156077258584 * exps;

	return u * erfs  / 2.0;
}

//	===	Evaluate autocorrelation Ck from a definite integral.
//	(f, d, s2)	parameters
//	k	distance 

double ck_eval(double f, double d, double s2, int k)
{
	if (k == 0)
		return 2.0 * d - 1.0;

	//	scale f and s2 as per theorem 1
	f = fmod( f * k, 1.0 );
	s2 *= (double) k;

	return 4.0 * ( eval_s1(f, d, s2, d) - eval_s1(f, d, s2, 0.0) - d ) + 1.0;
}

