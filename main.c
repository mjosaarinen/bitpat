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

	double f, d, s2;
	double ck[20];
	int i;
	size_t m_fft, m_sim;
	double h_fft, h_sim;
	
	f = 0.15;			//	fractional frequency
	d = 0.5;			//	bias towards 1 (0.5 = no bias)
	s2 = 0.04;			//	jitter variance
	
	m_fft = 1 << 12;	//	fft vector size (power of 2)
	m_sim = 10000000;	//	iteration count for simulations
	
	//	autocorrelation vectors
	
	printf("ck_sim() m= %zu\n", m_sim);
	ck_sim(ck, f, d, s2, 10, m_sim);

	for (i = 0; i < 10; i++) {
		printf("C_%d  %+18.15f  ( sim: %+9.6f )\n",
			i, ck_eval(f, d, s2, i), ck[i]);
	}
	printf("\n");
	
	//	print FFT and simulated distributions, entropy estimates

	for (i = 3; i <= 3; i++) {
		h_fft = entropy_fft(f, d, s2, i, m_fft, 1);
		h_sim = entropy_sim(f, d, s2, i, m_sim, 1);
		printf("Z%d  H_fft= %10.8f  H_sim= %10.8f\n\n", i, h_fft, h_sim);
	}
	
	return 0;
}
