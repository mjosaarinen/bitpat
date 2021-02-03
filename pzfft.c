//	pzfft.c
//	2021-01-31	Markku-Juhani O. Saarinen <mjos@pqshield.com>
//	Copyright (c) 2021, PQShield Ltd.  All rights reserved.

//	===	Compute entropy and bit pattern prob. via FFT convolutions.

#include <math.h>
#include <fftw3.h>
#include <stdio.h>

#include "bitpat.h"

//	===	f_s	step function.
//	vs	result, scaled as vs[i] = f_s(i/m)
//	m	length
//	f	fequency [0,1] (peak)
//	s2	jitter variance

double vec_fs(double *vs, size_t m, double f, double s2)
{
	double x, y, t, s;
	double tc, sum;
	size_t i;

	t = -0.5 / s2;
	tc = ceil(TAILCUT_TAU * sqrt(s2)) + 1.0;

	sum = 0.0;
	for (i = 0; i < m; i++) {
		x = ((double)i) / ((double)m);
		s = 0.0;
		for (y = x - f - tc; y < tc; y += 1.0) {
			s += exp(t * y * y);
		}
		sum += s;
		vs[i] = s;
	}

	//	expected sum: m * sqrt( 6.283185307179586476925286766559 * s2 );
	return sum;
}

//	chop. return sum over chopped area

double vec_chop(double *vr, const double *vs, size_t m, double d, int bit)
{
	double sum, t, u;
	size_t i, l;

	t = modf(d * ((double)m), &u);
	l = lrint(u);

	sum = 0.0;

	if (bit) {
		//	bit 1: x < D
		for (i = 0; i < l; i++) {
			u = vs[i];
			vr[i] = u;
			sum += u;
		}
		if (l >= m) return sum;

		u = vs[l] * t;
		vr[l] = u;
		sum += u;
		for (i = l + 1; i < m; i++) {
			vr[i] = 0.0;
		}

	} else {
		//	bit 0: x >= D
		for (i = 0; i < l; i++) {
			vr[i] = 0.0;
		}
		if (l >= m) return sum;

		u = vs[l] * (1.0 - t);
		vr[l] = u;
		sum += u;

		for (i = l + 1; i < m; i++) {
			u = vs[i];
			vr[i] = u;
			sum += u;
		}
	}

	return sum;
}

//	===	Evaluate min-entropy -log2(max p_z) using FFT (pzfft)
//	f	fequency [0,1] (peak)
//	d	cutoff (0.5 = no bias)
//	s2	jitter variance
//	n	Zn -- the bit sample size
//	m	FFT size (must be power of 2)
//	v	verbose (0 = print nothing, 1 = distribution to stdout)

double entropy_fft(double f, double d, double s2, size_t n, size_t m, int v)
{
	size_t i, j;
	double *vx, *vy;
	fftw_complex *vt, *vu;
	fftw_plan px, py, pi;
	double p, t, h1, pmax;
	uint64_t z;

	if (v) {
		printf("entropy_fft()  n= %zu  m= %zu\n", n, m);
	}
	
	//	plan px: vt = fft(vx)
	vx = (double *)fftw_malloc(sizeof(double) * m);
	vt = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * m);
	px = fftw_plan_dft_r2c_1d(m, vx, vt, FFTW_ESTIMATE);

	//	plan py: vu = fft(vy)
	vy = (double *)fftw_malloc(sizeof(double) * m);
	vu = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * m);
	py = fftw_plan_dft_r2c_1d(m, vy, vu, FFTW_ESTIMATE);

	//	plan pi: vx = fft^-1(vt)
	pi = fftw_plan_dft_c2r_1d(m, vt, vx, FFTW_ESTIMATE);

	//	compute, normalize, and transform the step function
	t = 1.0 / vec_fs(vy, m, f, s2);
	for (i = 0; i < m; i++) {
		vy[i] *= t;
	}
	fftw_execute(py);

	pmax = 0.0;
	h1 = 0.0;
	
	for (z = 0; z < (1UL << n); z++) {

		//	vx = uniform
		t = 1.0 / ((double)m);
		for (i = 0; i < m; i++) {
			vx[i] = t;
		}

		//	select bit
		p = vec_chop(vx, vx, m, d, z & 1);

		//	iterate over bits
		for (j = 1; j < n; j++) {
			//	perform convolution with the step function
			fftw_execute(px);
			for (i = 0; i < m; i++) {
				t = vt[i][0] * vu[i][0] - vt[i][1] * vu[i][1];
				vt[i][1] = vt[i][0] * vu[i][1] + vt[i][1] * vu[i][0];
				vt[i][0] = t;
			}
			fftw_execute(pi);

			//	scale
			t = 1.0 / ((double)m);
			for (i = 0; i < m; i++) {
				vx[i] *= t;
			}

			//	select bit
			p = vec_chop(vx, vx, m, d, (z >> j) & 1);
		}

		//	verbose
		if (v) {
			for (j = 0; j < n; j++) {
				putchar('0' + ((z >> (n - 1 - j)) & 1));
			}
			printf("  %18.16f\n", p);
		}
			
		//	entropies
		if (p > 0.0) {
			h1 -= p * log2(p);
		}
		if (p > pmax) {
			pmax = p;
		}
	}

	//	final stats
	if (v) {
		printf("Hm= %10.8f  H1= %10.8f  f= %8.6f  d= %8.6f  s2= %8.6f\n",
			  -log2(pmax) / ((double) n), h1 / ((double) n), f, d, s2);
	}

	//	free resources
	fftw_destroy_plan(px);
	fftw_destroy_plan(py);
	fftw_destroy_plan(pi);

	fftw_free(vx);
	fftw_free(vt);
	fftw_free(vy);
	fftw_free(vu);
	
	return -log2(pmax) / ((double) n);
}
