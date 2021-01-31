//	main.c
//	2021-01-31	Markku-Juhani O. Saarinen <mjos@pqshield.com>
//	Copyright (c) 2021, PQShield Ltd.  All rights reserved.

//	===	main() stub

#include <stdio.h>
#include <stdlib.h>

#include "bitpat.h"

int main(int argc, char **argv)
{
	(void)argc;
	(void)argv;

	double h, f, d, s2;
	size_t n;

	f = 0.1;
	d = 0.625;
	s2 = 0.04;
	n = 5;

	h = jitter_fft(f, d, s2, n, 1 << 10);
	printf("jitter_fft()  H = %f\n", h);

	h = jitter_sim(f, d, s2, n, 1000000);
	printf("jitter_sim()  H = %f\n", h);

	return 0;
}
