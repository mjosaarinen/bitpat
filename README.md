# bitpat

2021-02-01	Markku-Juhani O. Saarinen  <mjos@pqshield.com>

Supplementary material to the paper
*"On Entropy and Bit Patterns of Ring Oscillator Jitter."*

See [bitpat.h](bitpat.h) for function definitions.

Currently there's just a stub `main()` that prints 
calculated and simulated autocorrelated vectors for
some parameters, and then the same for bit pattern
distributions and entropy.

The code is plain C99, requires just [FFTW3](http://www.fftw.org):
```
sudo apt install libfftw3-dev
```

Cheers,
- markku
