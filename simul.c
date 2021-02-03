//	simul.c
//	2021-01-24	Markku-Juhani O. Saarinen <mjos@pqshield.com>
//	Copyright (c) 2021, PQShield Ltd.  All rights reserved.

//	===	Compute jitter probability distributions via simulation.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

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

//	=== Simulate byte strings into *zv. "zvlen" is length in bytes.

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

//	===	Estimate autocorrelation vector Ck[] from simulation.
//	*ck	result vector, each entry scaled to -1 <= ck[k] <= 1
//	(f, d, s2)	parameters
//	l	length of ck[] (at most 64)
//	m	number of simulation iterations (at least l)

void ck_sim(double *ck, double f, double d, double s2, size_t l, size_t m)
{
	size_t i, j, k;;
	double x, s;
	uint64_t z, sum[64];
	xcrand_t xcr;

	assert(l <= 64);
	assert(m > l);
	xcrand_init(&xcr);

	s = sqrt(s2);
	x = xcrand_d(&xcr);
	
	//	initialize
	z = 0;
	memset(sum, 0, sizeof(uint64_t) * l);

	//	m steps
	for (i = 0; i < m; i++) {

		JITTER_STEP

		//	ck[0] is simple bias
		sum[0] += z & 1;
		
		//	autocorrelation counts
		j = i < l ? i : l;
		for (k = 1; k < j; k++) {
			//	0 = different  1 = same  (xnor!)
			sum[k] += ~(z ^ (z >> k)) & 1;
		}
	}
	
	for (k = 0; k < l; k++) {
		ck[k] = 2.0 * (((double) sum[k]) / ((double) (m - k))) - 1.0;
	}
}

//	===	Estimate min-entropy -log2(max p_z) using simulation.
//	f	fequency [0,1] (peak)
//	d	cutoff (0.5 = no bias)
//	s2	jitter variance
//	n	Zn -- the bit sample size
//	m	number of iterations
//	v	verbose (0 = print nothing, 1 = distribution to stdout)

double entropy_sim(double f, double d, double s2, size_t n, size_t m, int v)
{
	size_t i, j, l;
	double x, s, p, h1, pmax;
	uint64_t z, *cnt;

	xcrand_t xcr;

	if (v) {
		printf("entropy_sim()  n= %zu  m= %zu\n", n, m);
	}
	
	xcrand_init(&xcr);

	s = sqrt(s2);
	x = xcrand_d(&xcr);

	l = 1UL << n;
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

	pmax = 0.0;
	h1 = 0.0;

	for (i = 0; i < l; i++) {

		p = ((double) cnt[i]) / ((double) m);

		//	verbose	
		if (v) {
			for (j = 0; j < n; j++) {
				putchar('0' + ((i >> (n - 1 - j)) & 1));
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
	free(cnt);
	
	//	final stats
	if (v) {
		printf("Hm= %10.8f  H1= %10.8f  f= %8.6f  d= %8.6f  s2= %8.6f\n",
			  -log2(pmax) / ((double) n), h1 / ((double) n), f, d, s2);
	}
	
	return -log2(pmax) / ((double) n);
}
