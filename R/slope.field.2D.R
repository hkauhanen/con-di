slope.field.2D <- function(lower_bounds = c(0, 0),
                           upper_bounds = c(1, 1),
                           resolution = 20,
                           len = 0.1,
                           normalize = TRUE,
                           FUN,
                           ...) {
  do_one <- function(x,
                     y,
                     FUN,
                     ...) {
    df <- FUN(x=c(x, y), ...)
    names(df) <- c("dotx1", "dotx2")
    df$x1 <- x
    df$x2 <- y

    next1 <- x + df$dotx1
    next2 <- y + df$dotx2
    if (normalize) {
      length <- sqrt((next1 - x)^2 + (next2 - y)^2)
    } else {
      length <- 1
    }

    if (length == 0) {
      df$x1end <- x
      df$x2end <- y
    } else {
      df$x1end <- x + (len/length)*df$dotx1
      df$x2end <- y + (len/length)*df$dotx2
    }

    as.data.frame(df)
  }

  mesh <- expand.grid(x1=seq(from=lower_bounds[1], to=upper_bounds[1], length.out=resolution), x2=seq(from=lower_bounds[2], to=upper_bounds[2], length.out=resolution))

  do.call(rbind, mapply(do_one, mesh$x1, mesh$x2, MoreArgs=c(FUN=FUN, list(...)), SIMPLIFY=FALSE))
}
