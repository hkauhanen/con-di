# Convergence--divergence game: time derivative
#
# Henri Kauhanen 2020


# We also refer to the convergence--divergence model as the "vergence" model.
# The function that supplies the time derivative is thus "verger". The
# derivatives are here calculated in a vectorized fashion, utilizing
# matrix operations. Difficult to read, but short (and efficient) code.
#
verger <- function(x,
                   S,
                   m,
                   tm) {
  dS <- diag(S)
  S0 <- S
  diag(S0) <- 0

  # fitness of variant 1
  f <- x*(dS + S0 %*% (1-x))

  # fitness of variant 2
  tf <- (1-x)*(dS + S0 %*% x)

  # time derivative
  (1-x - tm)*x*f - (x - m)*(1-x)*tf
}
