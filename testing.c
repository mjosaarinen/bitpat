//	testing.c
//	2021-01-31	Markku-Juhani O. Saarinen <mjos@pqshield.com>
//	Copyright (c) 2021, PQShield Ltd.  All rights reserved.

//	===	things for experiments

#include <math.h>
#include <fftw3.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bitpat.h"
#include "xcrand.h"

//	debug hex string

void hexvec(const uint8_t *x, size_t xlen, const char *lab)
{
	size_t i;

	printf("%s[%zu] = ", lab, xlen);

	for (i = 0; i < xlen; i++) {
		printf("%02X", x[i]);
	}
	printf("\n");
}

void realvec(const double *v, size_t n, const char *lab)
{
	size_t i;

	for (i = 0; i < n; i++) {
		printf("%s[%4zu]  %17.14f\n", lab, i, v[i]);
	}
}

void clxvec(const fftw_complex *v, size_t n, const char *lab)
{
	size_t i;

	for (i = 0; i < n; i++) {
		printf("%s[%4zu] ( %17.14f, %17.14f )\n", lab, i, v[i][0], v[i][1]);
	}
}

//	conventional convolution

void vconv(double *z, const double *x, const double *y, size_t n)
{
	size_t i, j;

	for (i = 0; i < n; i++) {
		z[i] = 0.0;
	}
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			z[(i + j) % n] += x[i] * y[j];
		}
	}
}

//	write a file

size_t wrfile(const char *path, const uint8_t *x, size_t xlen)
{
	FILE *f;
	size_t wlen;

	f = fopen(path, "wb");
	if (f == NULL) {
		perror(path);
		return 0;
	}

	wlen = fwrite(x, 1, xlen, f);
	if (ferror(f)) perror(path);
	fclose(f);

	return wlen;
}

//	write in binary

size_t wrbin(const char *path, const uint8_t *x, size_t xlen)
{
	FILE *f;
	size_t i, j;
	uint8_t buf[8];

	f = fopen(path, "wb");
	if (f == NULL) {
		perror(path);
		return 0;
	}

	for (i = 0; i < xlen; i++) {
		for (j = 0; j < 8; j++) {
			buf[7 - j] = '0' + ((x[i] >> j) & 1);
		}
		if (fwrite(buf, 1, 8, f) != 8) break;
	}

	if (ferror(f)) perror(path);
	fclose(f);

	return 8 * i;
}

//	used to generate the test files for SP 800-90B suites

#ifndef LEN_90B
#define LEN_90B 1000000
#endif

int gen_90b_dat()
{
	double f, d, s;
	xcrand_t xcr;
	char fn[100];
	uint8_t dat[LEN_90B];

	xcrand_init(&xcr);

	//	randomize parameters
	s = xcrand_d(&xcr);
	f = xcrand_d(&xcr);
	d = 0.5;

	snprintf(fn, sizeof(fn), "s%08ld-f%08ld.dat", lrint(100000000.0 * s),
			 lrint(100000000.0 * f));

	zbytes(dat, LEN_90B, f, d, s * s);
	printf("%s\n", fn);
	wrfile(fn, dat, LEN_90B);

	//	strncat(fn, ".bin", sizeof(fn));
	//	printf("%s\n", fn);
	//	wrbin(fn, dat, LEN);

	return 0;
}

//	===	Estimate min-entropy -log2(max p_z) by depth-first search on z
//	f	fequency [0,1] (peak)
//	d	cutoff (0.5 = no bias)
//	s2	jitter variance
//	n	Zn -- the bit sample size
//	m	FFT size (must be power of 2)
//	dfs	strategy
//	v	verbose (0 = print nothing, 1 = distribution to stdout)

double entropy_dfs(double f, double d, double s2, size_t n, size_t m, 
	int v, dfs_how_t how)
{
	int bit;
	size_t i, j;
	double *vx, *vy;
	double *g0, *g1;
	fftw_complex *vt, *vu;
	fftw_plan px, py, pz;

	double h, q0, q1;
	double x, r, t;

	//	distributions  vt = fft(vx)
	g0 = (double *)fftw_malloc(sizeof(double) * m);
	g1 = (double *)fftw_malloc(sizeof(double) * m);
	vx = (double *)fftw_malloc(sizeof(double) * m);
	vt = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * m);
	px = fftw_plan_dft_r2c_1d(m, vx, vt, FFTW_ESTIMATE);

	//	step function  vu = fft(vy)
	vy = (double *)fftw_malloc(sizeof(double) * m);
	vu = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * m);
	py = fftw_plan_dft_r2c_1d(m, vy, vu, FFTW_ESTIMATE);

	//	convolution result  vx = fft^-1(vx)
	pz = fftw_plan_dft_c2r_1d(m, vt, vx, FFTW_ESTIMATE);

	//	compute transformed  step function
	r = vec_fs(vy, m, f, s2);
	fftw_execute(py);

	//	uniform start
	for (i = 0; i < m; i++) {
		vx[i] = 1.0;
	}

	//	for deterministic step, xenter of the larger bias
	x = d >= 0.5 ? d / 2.0 : d + 0.5 * (1.0 - d);

	//	iterate over bits

	h = 0.0;
	for (j = 0; j < n; j++) {
	

		//	convolution with the step function

		if (j > 0) {		
			fftw_execute(px);
			for (i = 0; i < m; i++) {
				t = vt[i][0] * vu[i][0] - vt[i][1] * vu[i][1];
				vt[i][1] = vt[i][0] * vu[i][1] + vt[i][1] * vu[i][0];
				vt[i][0] = t;
			}
			fftw_execute(pz);
		}
		
		//	select bit, normalize

		q0 = vec_chop(g0, vx, m, d, 0);
		q1 = vec_chop(g1, vx, m, d, 1);

		//	track a point
		x += f;
		x -= floor(x);

		//	strategy
		switch (how) {		
				
			case DFS_BIG_MASS:
				bit = q0 > q1 ? 0 : 1;
				break;

			default:
			case DFS_FOLLOW_X:
				bit = x < d ? 1 : 0;
				break;
		}

		//	bit selection
		if (bit == 0) {
			h -= log2(q0 / (q0 + q1));
			r = 1.0 / q0;
			for (i = 0; i < m; i++) {
				vx[i] = r * g0[i];
			}
		} else {
			h -= log2(q1 / (q0 + q1));
			r = 1.0 / q1;
			for (i = 0; i < m; i++) {
				vx[i] = r * g1[i];
			}
		}
		t = (q0 + q1);
		if (v) {
			printf("[%3zu] %d  q0= %16.14f  q1= %16.14f  h= %16.8f\n",
				j, bit, q0 / t, q1 / t, h / ((double)j + 1));
		}
	}

	fftw_destroy_plan(px);
	fftw_destroy_plan(py);
	fftw_destroy_plan(pz);

	fftw_free(g0);
	fftw_free(g1);
	fftw_free(vx);
	fftw_free(vt);
	fftw_free(vy);
	fftw_free(vu);

	return h / ((double) n);
}

//	create data for scatterplot

void gen_estim()
{
	double f, d, s;
	xcrand_t xcr;
	size_t i, n, m;
	
	double h_mass, h_usex;
	
	xcrand_init(&xcr);

	n = 100;
	m = 1 << 10;

	for (i = 0; i < 10000; i++) {

		f = xcrand_d(&xcr);
		d = 0.5;
		s = 0.5 * xcrand_d(&xcr);

		h_mass = entropy_dfs(f, d, s*s, n, m, 0, DFS_BIG_MASS);
		h_usex = entropy_dfs(f, d, s*s, n, m, 0, DFS_FOLLOW_X);
		
		printf("s= %10.8f f= %10.7f  h_mass= %10.8f  h_usex= %10.8f  %s\n",
			s, f, h_mass, h_usex,
			h_mass < h_usex ? "mass" : "usex");
	}
}


