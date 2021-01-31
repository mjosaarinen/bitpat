//	jitter.c
//	2021-01-24	Markku-Juhani O. Saarinen <mjos@pqshield.com>
//	Copyright (c) 2021, PQShield Ltd.  All rights reserved.

//	===	Compute jitter probability distributions via simulation.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "xcrand.h"
#include "bitpat.h"

//	simulation step

#define JITTER_STEP                    \
	{                                  \
		x += f + s * xcrand_std(&xcr); \
		x -= floor(x);                 \
		z <<= 1;                       \
		if (x < d) z ^= 1;             \
	}

//	===	Simulate byte strings into *zv. "zvlen" is length in bytes.

void zbytes(uint8_t *zv, size_t zvlen, double f, double d, double s2)
{
	size_t i, j;
	double x, s;
	uint32_t z;
	xcrand_t xcr;

	xcrand_init(&xcr);

	s = sqrt(s2);
	x = xcrand_d(&xcr);

	z = 0;
	for (i = 0; i < zvlen; i++) {
		for (j = 0; j < 8; j++) {
			JITTER_STEP
		}
		zv[i] = z & 0xFF;
	}
}

//	===	Simulate min-entropy -log2(max p_z)
//	f	fequency [0,1] (peak)
//	d	cutoff  [0,1] (0.5 = no bias)
//	s2	jitter variance
//	n	Zn -- the bit sample size
//	m	number of iterations

//	simulation

double jitter_sim(double f, double d, double s2, size_t n, size_t m)
{
	size_t i, l;
	int j;
	double x, s, p, hm, h1, pmax;
	uint64_t *cnt;
	uint32_t z;

	xcrand_t xcr;

	xcrand_init(&xcr);

	s = sqrt(s2);
	x = xcrand_d(&xcr);

	l = 1L << n;
	cnt = calloc(l, sizeof(uint64_t));
	if (cnt == NULL) {
		perror("calloc()");
		return 0.0;
	}

	//	compute initial z
	z = 0;
	for (i = 0; i < m; i++) {
		JITTER_STEP
	}

	//	main loop (frequency count)
	for (i = 0; i < m; i++) {
		JITTER_STEP
		z &= (l - 1);
		cnt[z]++;
	}

	//	print analysis
	pmax = 0.0;
	h1 = 0.0;

	for (i = 0; i < l; i++) {
		for (j = n; j > 0; j--) {
			putchar('0' + ((i >> (j - 1)) & 1));
		}
		p = ((double)cnt[i]) / ((double)m);

		printf("  %18.16f\n", p);

		if (p > 0.0) {
			h1 -= p * log(p) / log(2.0);
		}
		if (p > pmax) {
			pmax = p;
		}
	}
	hm = -log(pmax) / log(2.0);

	printf("H1= %10.8f  Hm= %10.8f  f= %8.6f  d= %8.6f  s2= %8.6f  n= %zu\n",
		   h1 / ((double)n), hm / ((double)n), f, d, s2, n);
	free(cnt);

	return -log2(pmax) / ((double) n);
}
