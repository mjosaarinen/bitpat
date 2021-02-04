# bitpat

2021-02-01	Markku-Juhani O. Saarinen  <mjos@pqshield.com>

Supplementary material for the paper Saarinen, M.-J. (2021):
*"On Entropy and Bit Patterns of Ring Oscillator Jitter."*
available as [arXiv:2102.02196](https://arxiv.org/abs/2102.02196).

The code is plain C99, requires just [FFTW3](http://www.fftw.org):
```
sudo apt install libfftw3-dev
```

See [bitpat.h](bitpat.h) for function definitions. Each function takes
in (F, D, s2) ring oscillator model parameters. Some main functions:

*	[pzfft.c](pzfft.c): `entropy_fft()` computes bit pattern probabilities
	and Shannon / min-entropy from it.
*	[autoc.c](autoc.c): `ck_eval()` computes the autocorrelation coefficients
	of Eqn. 6, for which `eval_s1()` computes the indefinite integral of
	Eqn. 12.
*	[simul.c](simul.c) has functions `entropy_sim()` and `ck_sim()` which
	obtain distribution, entropy, and autocorrelation estimates via 
	(much slower) Monte Carlo simulation.

There's a stub [main.c](main.c) that prints calculated and simulated
autocorrelated vectors for some parameters, and then does the same for bit
pattern distributions and entropy.

Cheers,
- markku
