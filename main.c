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

	double f, d, s;

	//	print FFT and simulated values
	f = 0.15;
	s = 0.3;
	d = 0.5;
	
	entropy_fft(f, d, s*s, 3, 1 << 10, 1);
	entropy_sim(f, d, s*s, 3, 10000000, 1);

	entropy_fft(f, d, s*s, 4, 1 << 10, 1);
	entropy_sim(f, d, s*s, 4, 10000000, 1);

	entropy_fft(f, d, s*s, 5, 1 << 10, 1);
	entropy_sim(f, d, s*s, 5, 10000000, 1);

	return 0;
}
