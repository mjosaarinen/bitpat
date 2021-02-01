# bitpat

2021-02-01	Markku-Juhani O. Saarinen  <mjos@pqshield.com>

Complementary material to the paper
"On Entropy and Bit Patterns of Ring Oscillator Jitter."

See [bitpat.h](bitpat.h) for function definitions.
Currently there's just a stub `main()` that prints both 
calculated (`entropy_fft`) and simulated (`entropy_sim`)
probability distributions for some parameters.

Plain C code, requires FFTW3
```
sudo apt install libfftw3-dev
```

Cheers,
- markku
