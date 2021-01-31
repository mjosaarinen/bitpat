//	xcrand.h
//	2021-01-29	Markku-Juhani O. Saarinen <mjos@pqshield.com>
//	Copyright (c) 2021, PQShield Ltd.  All rights reserved.

//	===	Cryptographic PRNG for numerical simulations. ===

#ifndef _XCRAND_H_
#define _XCRAND_H_

#include <stdint.h>
#include <stddef.h>

//	context structure

typedef struct {
	uint32_t seed[16];
	uint8_t obuf[64];
	int optr, xflg;
	double x;
} xcrand_t;

//	initialize with random seed
int xcrand_init(xcrand_t *xcr);

//	initialize with given seed
void xcrand_seed(xcrand_t *xcr, uint8_t *seed, size_t len);

//	get "n" bytes at r
void xcrand_bytes(xcrand_t *xcr, void *r, size_t n);

//	get a single random byte
uint8_t xcrand_u8(xcrand_t *xcr);

//	random 32-bit unsigned integer
uint64_t xcrand_u32(xcrand_t *xcr);

//	random 64-bit unsigned integer
uint64_t xcrand_u64(xcrand_t *xcr);

//	random in interval [0,1) (IEEE 784 double)
double xcrand_d(xcrand_t *xcr);

//	sample normal distribution  N(0,1)
double xcrand_std(xcrand_t *xcr);

#endif
