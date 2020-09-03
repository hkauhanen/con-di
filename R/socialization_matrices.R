# Socialization matrices for use with the convergence-divergence game.
#
# Henri Kauhanen 2020


# Simple 2D case
#
S2D <- function(s11 = 1,
                s12,
                s21,
                s22 = 1) {
  matrix(c(s11, s12,
           s21, s22),
         byrow=TRUE,
         nrow=2)
}


# Jocks, Burnouts, Adults
#
JBA <- function(sigma,
                tau1,
                tau2,
                upsilon) {
  matrix(c(1, sigma, tau1,
           sigma, 1, tau2,
           upsilon, upsilon, 1),
         byrow=TRUE,
         nrow=3)
}


