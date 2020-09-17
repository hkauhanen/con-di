# Convergence–divergence game

This repository contains code required to reproduce the results presented in the following paper:

> Kauhanen, H. (2020) A model of linguistic identity dynamics. Ms., University of Konstanz.

Please cite this publication if you find this useful.


## Symbolic calculation of equilibria, Jacobian and eigenvalues

The system's equilibria, the Jacobian matrix, and the Jacobian's eigenvalues can be computed in the case of the fully symmetric two-strategy, two-population game. Since carrying out the calculations by hand would be unwieldy, they have been implemented in the Maxima CAS (Computer Algebra System). A [wxmaxima](https://wxmaxima-developers.github.io/wxmaxima/) worksheet is included in the `wxmaxima` directory.


## Numerical solutions of asymmetric 3-population game

To produce the numerical solutions of the particular instance of the game studied in Section 6 of the paper, navigate to the `julia` directory, and, from the Julia REPL, issue the following commands:

``` julia
include("con-di.jl")
Random.seed!(2020)
popthreesweep("../results/sweep.csv")
popthreesweepone("../results/combo.csv")
```

The Julia packages DelimitedFiles, DifferentialEquations, LinearAlgebra and Random are required. Tested with Julia 1.4.2.


## Plots

The figures in the above-cited paper were produced using a combination of R/ggplot2 and LaTeX/PGF/TikZ. To reproduce these, first navigate into the `Rsession` directory and issue the following commands from R (the packages ggplot2, RColorBrewer, reshape2 and viridis are required):

``` R
source("../R/plots.R")
set.seed(2020)
plot.all()
```

This produces temporary plots which serve as input to the LaTeX/PGF/TikZ routines. The latter may be run using the `do_plots.sh` shell script in the `tex` directory:

``` bash
sh do_plots.sh
```

Output is copied into the `plots` directory. Note that an installation of LaTeX (and pdflatex) is required, as are a number of LaTeX packages. In case of trouble, consult the error messages and install the required packages.


## Acknowledgement

Preparation of this software was made possible by financial support from the Federal Ministry of Education and Research (BMBF) and the Baden-Württemberg Ministry of Science as part of the Excellence Strategy of the German Federal and State Governments.
