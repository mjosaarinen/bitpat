//	xcrand.c
//	2021-01-29	Markku-Juhani O. Saarinen <mjos@pqshield.com>
//	Copyright (c) 2021, PQShield Ltd.  All rights reserved.

//	===	Cryptographic PRNG for numerical simulations. ===

#include "xcrand.h"

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/random.h>

//	ChaCha "quarter round"

#define ROL32(v, n) (((v) << (n)) | ((v) >> (32 - (n))))

#define CHACHA_QR(a, b, c, d)      \
	x[a] += x[b];                  \
	x[d] = ROL32(x[d] ^ x[a], 16); \
	x[c] += x[d];                  \
	x[b] = ROL32(x[b] ^ x[c], 12); \
	x[a] += x[b];                  \
	x[d] = ROL32(x[d] ^ x[a], 8);  \
	x[c] += x[d];                  \
	x[b] = ROL32(x[b] ^ x[c], 7);

static void chacha_core(uint8_t r[64], const uint32_t seed[16])
{
	uint32_t x[16];
	int i;

	memcpy(x, seed, sizeof(x));
	for (i = 20; i > 0; i -= 2) {
		CHACHA_QR(0, 4, 8, 12)
		CHACHA_QR(1, 5, 9, 13)
		CHACHA_QR(2, 6, 10, 14)
		CHACHA_QR(3, 7, 11, 15)
		CHACHA_QR(0, 5, 10, 15)
		CHACHA_QR(1, 6, 11, 12)
		CHACHA_QR(2, 7, 8, 13)
		CHACHA_QR(3, 4, 9, 14)
	}
	for (i = 0; i < 16; ++i) {
		x[i] += seed[i];
	}
	memcpy(r, x, sizeof(x));
}

//	initialize with given seed

void xcrand_seed(xcrand_t *xcr, uint8_t *seed, size_t len)
{
	size_t i, l;

	l = len;
	if (l > 62) l = 62;
	memcpy(xcr->seed, seed, l);
	memset(((uint8_t *)xcr->seed) + l, 0, 62 - l);
	xcr->seed[15] = l;
	len -= l;
	seed += l;

	//	yes, it's an over-engineered "permutation" hash mode :)
	while (len > 0) {
		chacha_core(xcr->obuf, xcr->seed);
		l = len;
		if (l > 32) l = 32;
		for (i = 0; i < l; i++) xcr->obuf[i] ^= seed[i];
		memcpy(xcr->seed, xcr->obuf, 64);
		xcr->seed[15] ^= l;
		len -= l;
		seed += l;
	}

	xcr->optr = 64;
	xcr->xflg = 0;
}

//	initialize with random seed

int xcrand_init(xcrand_t *xcr)
{
	if (getrandom(xcr->seed, 64, 0) != 64) {
		perror("xcrand_init(): getrandom() fail");
		return -1;
	}
	xcr->optr = 64;
	xcr->xflg = 0;

	return 0;
}

//	users can do big steps manually

static inline void xcrand_step_1(xcrand_t *xcr)
{
	int i;

	for (i = 0; i < 16; i++) {
		xcr->seed[i]++;
		if (xcr->seed[i] != 0) break;
	}
}

//	get some bytes

void xcrand_bytes(xcrand_t *xcr, void *r, size_t n)
{
	size_t i, l;
	uint8_t *p = (uint8_t *)r;

	i = xcr->optr;
	if (i < 64) {
		l = 64 - i;
		//	handle the small-increment case right here
		if (l >= n) {
			memcpy(p, &xcr->obuf[i], n);
			xcr->optr = i + n;
			return;
		}
		//	use up rest of buffer
		memcpy(p, &xcr->obuf[i], l);
		p += l;
		n -= l;
	}
	l = 64;

	//	copy blocks
	while (n > 0) {
		chacha_core(xcr->obuf, xcr->seed);
		xcrand_step_1(xcr);

		if (n < 64) l = n;
		memcpy(p, xcr->obuf, l);
		p += l;
		n -= l;
	}
	xcr->optr = l;
}

//	get a single random byte

uint8_t xcrand_u8(xcrand_t *xcr)
{
	if (xcr->optr < 64) {
		return xcr->obuf[xcr->optr++];
	}
	chacha_core(xcr->obuf, xcr->seed);
	xcrand_step_1(xcr);

	xcr->optr = 1;
	return xcr->obuf[0];
}

//	random 32-bit unsigned integer

uint64_t xcrand_u32(xcrand_t *xcr)
{
	uint32_t x;
	xcrand_bytes(xcr, &x, sizeof(uint32_t));

	return x;
}

//	random 64-bit integer

uint64_t xcrand_u64(xcrand_t *xcr)
{
	uint64_t x;
	xcrand_bytes(xcr, &x, sizeof(uint64_t));

	return x;
}

//	random in interval [0,1) (IEEE 784 double)

double xcrand_d(xcrand_t *xcr)
{
	return ((double)xcrand_u64(xcr)) / 18446744073709551616.0;
}

//	sample normal distribution  N(0,1)

double xcrand_std(xcrand_t *xcr)
{
	double x, y, r2;

	//	saved
	if (xcr->xflg) {
		xcr->xflg = 0;
		return xcr->x;
	}

	do {
		x = 2.0 * xcrand_d(xcr) - 1.0;
		y = 2.0 * xcrand_d(xcr) - 1.0;
		r2 = x * x + y * y;
	} while (r2 > 1.0 || r2 == 0.0);
	r2 = sqrt(-2.0 * log(r2) / r2);

	x *= r2;
	y *= r2;

	xcr->x = x;
	xcr->xflg = 1;

	return y;
}
