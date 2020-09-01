# Convergence-divergence game

This repository contains code required to reproduce the results presented in the following paper:

> Kauhanen, H. (2020) Replicator–mutator dynamics of linguistic convergence and divergence. Ms., University of Konstanz.

Please cite this publication if you find this useful.


## Symbolic calculation of Jacobian and eigenvalues

The Jacobian matrix and its eigenvalues can be computed in the case of the fully symmetric two-strategy, two-population game. Since carrying out the calculations by hand would be unwieldy, they have been implemented in the Maxima CAS (Computer Algebra System). A [wxmaxima](https://wxmaxima-developers.github.io/wxmaxima/) worksheet is included in the `wxmaxima` directory.


## Numerical solutions of asymmetric 3-population game

To produce the numerical solutions of the particular instance of the game studied in Section REF of the paper, navigate to the `julia` directory, and, from the Julia REPL, issue the following commands:

``` julia
include("con-di.jl")
Random.seed!(2020)
popthreesweep("../results/sweep.csv")
popthreesweep("../results/combo.csv")
```

Tested with Julia 1.4.2.


## Plots


## Acknowledgements

Preparation of this software was supported by the Federal Ministry of Education and Research (BMBF) and the Baden-Württemberg Ministry of Science as part of the Excellence Strategy of the German Federal and State Governments.
