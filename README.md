# PowerIterations.jl
Compute the largest eigenvalues of a diagonalisable matrix with non-degenerate eigenvalues

Currently only works for Hermitian matrices, i.e. real eigenvalues.

I noticed that there is another implementation of the power iteration method in IterativeSolvers.jl.

References:
Richard von Mises and H. Pollaczek-Geiringer, Praktische Verfahren der Gleichungsauflösung, ZAMM - Zeitschrift für Angewandte Mathematik und Mechanik 9, p. 152-164 (1929)

An extension for complex eigen values can be found on the [wikipedia page](https://en.wikipedia.org/wiki/Power_iteration) using the Rayleigh quotient.

Note: would be good to merge this with the official IterativeSolvers.jl?

