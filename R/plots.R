# Required packages
require(ggplot2)
require(RColorBrewer)
require(reshape2)
require(viridis)


# Required code
source("../R/slope.field.2D.R")
source("../R/socialization_matrices.R")
source("../R/verger.R")


# Font family
deffam <- "Times"
#deffam <- "Palatino"

# Dimensions for phase space plots
spacesize <- 2.5
spacesizeh <- 2.4

# Palettes
mypal <- brewer.pal(4, "Set3")


# Plot all figures for manuscript
plot.all <- function() {
  # Plot "raw" plots, i.e. the input to beautification by TeX
  pdf("sweep.pdf", width=6, height=1.95, family=deffam)
  print(plot_numerics_sweep())
  dev.off()
  embedFonts("sweep.pdf")

  pdf("indiv.pdf", width=5.3, height=1.95, family=deffam)
  print(plot_numerics_individualpoint())
  dev.off()
  embedFonts("indiv.pdf")

  pdf("phaseI.pdf", height=spacesizeh, width=spacesize, family=deffam)
  print(plot.space.FS(which="I"))
  dev.off()
  embedFonts("phaseI.pdf")

  pdf("phaseII.pdf", height=spacesizeh, width=spacesize, family=deffam)
  print(plot.space.FS(which="II"))
  dev.off()
  embedFonts("phaseII.pdf")

  pdf("phaseIII.pdf", height=spacesizeh, width=spacesize, family=deffam)
  print(plot.space.FS(which="III"))
  dev.off()
  embedFonts("phaseIII.pdf")

  pdf("phaseIV.pdf", height=spacesizeh, width=spacesize, family=deffam)
  print(plot.space.FS(which="IV"))
  dev.off()
  embedFonts("phaseIV.pdf")

  pdf("phases.pdf", width=5, height=2.5, family=deffam)
  print(plot.phases.FS())
  dev.off()
  embedFonts("phases.pdf")

  pdf("eval.pdf", width=4, height=4, family=deffam)
  print(empirical_eval())
  dev.off()
  embedFonts("eval.pdf")
}


# Nice-looking scientific notation; from https://stackoverflow.com/questions/11610377/how-do-i-change-the-formatting-of-numbers-on-an-axis-with-ggplot
fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  l <- gsub("0e\\+00", "0", l)

  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)

  # remove + after exponent
  l <- gsub("e\\+", "e", l)

  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  #l <- gsub("\\'1[\\.0]*\\'\\%\\*\\%", "", l)

  # return this as an expression
  parse(text=l)
}


# Plot results of the numerical solution sweep (grid)
plot_numerics_sweep <- function(df = read.csv("../results/sweep.csv")) {
  df$type <- ifelse(df$x1 < 0.5 & df$x2 > 0.5, 1, 0)
  df <- aggregate(df[, 7], by=list(df$sigma, df$upsilon, df$k), FUN=mean)
  names(df) <- c("sigma", "upsilon", "k", "type")

  g <- ggplot(df, aes(x=sigma, y=k, fill=type)) + geom_tile() + facet_wrap(.~upsilon, nrow=1)
  g <- g + scale_fill_viridis(name="")
  g <- g + theme_bw() + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), panel.spacing = unit(1, "lines"), strip.background = element_blank(), strip.text.x = element_blank(), axis.text=element_text(color="black"), legend.position="right")
  g <- g + scale_y_continuous(limits=c(0,5), expand=c(0,0))
  g <- g + scale_x_continuous(limits=c(0,5), expand=c(0,0))
  g <- g + xlab("") + ylab("")

  df3 <- data.frame(sigma=4, k=4, type=0, upsilon=1)
  g <- g + geom_point(data=df3, inherit.aes=FALSE, aes(x=sigma, y=k), pch=8, size=2.0)

  g
}


# Plots numerical solutions for one parameter combination
plot_numerics_individualpoint <- function(df = read.csv("../results/combo.csv")) {
  df <- reshape2::melt(df, id.vars=c("sigma", "upsilon", "k", "id", "t"))
  df <- df[df$id %in% sample(1:100, size=10, replace=FALSE), ]

  g <- ggplot(df, aes(x=t, y=value, color=factor(id))) + geom_line(alpha=0.8) + facet_wrap(.~variable, ncol=3)
  g <- g + theme_bw() + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), panel.spacing = unit(1, "lines"), strip.background = element_blank(), strip.text.x = element_blank(), axis.text=element_text(color="black"), legend.position="none")
  g <- g + scale_color_brewer(palette="Spectral")
  g <- g + scale_x_log10()
  g <- g + ylim(0,1)
  g <- g + annotation_logticks(sides="bt", size=0.3)
  g <- g + xlab("") + ylab("")

  g
}


# Equilibria on the tilde-D diagonal. Alpha means sigma in current notation.
tDeq <- function(mu, alpha) {
  wtD <- sqrt(-(3*alpha^2 + 4*alpha)*mu^2 - 2*(alpha^2 + 3*alpha + 2)*mu + (alpha+1)^2)/(2*(alpha*mu + alpha + 1))
  c(0.5 + wtD, 0.5 - wtD)
}

# Equilibria on the D-diagonal. Alpha means sigma in current notation.
Deq <- function(mu, alpha) {
  wD <- sqrt((alpha^2 + 4*alpha)*mu^2 - 2*(alpha + 2)*mu + 1)/(2*(alpha*mu - 1))
  c(0.5 + wD, 0.5 - wD)
}

# Off-diagonal equilibria. Taken from Maxima output.
offeq <- function(alpha, mu) {
  a <- alpha

  first <- (mu*(a*(4*sqrt(a^2+4*a+4)*sqrt(-2*sqrt(a^4*mu^4+(2*a^4+8*a^3+8*a^2)*mu^3+(-a^4+2*a^3+18*a^2+32*a+16)*mu^2+(-2*a^3-10*a^2-16*a-8)*mu+a^2+2*a+1)-2*a^2*mu^2+(-2*a^2-8*a-8)*mu+a^2+2*a+2)+4*sqrt(a^4*mu^4+(2*a^4+8*a^3+8*a^2)*mu^3+(-a^4+2*a^3+18*a^2+32*a+16)*mu^2+(-2*a^3-10*a^2-16*a-8)*mu+a^2+2*a+1)-4)+a^2*(sqrt(a^2+4*a+4)*sqrt(-2*sqrt(a^4*mu^4+(2*a^4+8*a^3+8*a^2)*mu^3+(-a^4+2*a^3+18*a^2+32*a+16)*mu^2+(-2*a^3-10*a^2-16*a-8)*mu+a^2+2*a+1)-2*a^2*mu^2+(-2*a^2-8*a-8)*mu+a^2+2*a+2)+2*sqrt(a^4*mu^4+(2*a^4+8*a^3+8*a^2)*mu^3+(-a^4+2*a^3+18*a^2+32*a+16)*mu^2+(-2*a^3-10*a^2-16*a-8)*mu+a^2+2*a+1)-6)+4*sqrt(a^2+4*a+4)*sqrt(-2*sqrt(a^4*mu^4+(2*a^4+8*a^3+8*a^2)*mu^3+(-a^4+2*a^3+18*a^2+32*a+16)*mu^2+(-2*a^3-10*a^2-16*a-8)*mu+a^2+2*a+1)-2*a^2*mu^2+(-2*a^2-8*a-8)*mu+a^2+2*a+2)-a^4-4*a^3)+mu^2*(a^2*(24-sqrt(a^2+4*a+4)*sqrt(-2*sqrt(a^4*mu^4+(2*a^4+8*a^3+8*a^2)*mu^3+(-a^4+2*a^3+18*a^2+32*a+16)*mu^2+(-2*a^3-10*a^2-16*a-8)*mu+a^2+2*a+1)-2*a^2*mu^2+(-2*a^2-8*a-8)*mu+a^2+2*a+2))+3*a^4+14*a^3+16*a)+a*(-sqrt(a^2+4*a+4)*sqrt(-2*sqrt(a^4*mu^4+(2*a^4+8*a^3+8*a^2)*mu^3+(-a^4+2*a^3+18*a^2+32*a+16)*mu^2+(-2*a^3-10*a^2-16*a-8)*mu+a^2+2*a+1)-2*a^2*mu^2+(-2*a^2-8*a-8)*mu+a^2+2*a+2)+2*sqrt(a^4*mu^4+(2*a^4+8*a^3+8*a^2)*mu^3+(-a^4+2*a^3+18*a^2+32*a+16)*mu^2+(-2*a^3-10*a^2-16*a-8)*mu+a^2+2*a+1)+2)-sqrt(a^2+4*a+4)*sqrt(a^4*mu^4+(2*a^4+8*a^3+8*a^2)*mu^3+(-a^4+2*a^3+18*a^2+32*a+16)*mu^2+(-2*a^3-10*a^2-16*a-8)*mu+a^2+2*a+1)*sqrt(-2*sqrt(a^4*mu^4+(2*a^4+8*a^3+8*a^2)*mu^3+(-a^4+2*a^3+18*a^2+32*a+16)*mu^2+(-2*a^3-10*a^2-16*a-8)*mu+a^2+2*a+1)-2*a^2*mu^2+(-2*a^2-8*a-8)*mu+a^2+2*a+2)-sqrt(a^2+4*a+4)*sqrt(-2*sqrt(a^4*mu^4+(2*a^4+8*a^3+8*a^2)*mu^3+(-a^4+2*a^3+18*a^2+32*a+16)*mu^2+(-2*a^3-10*a^2-16*a-8)*mu+a^2+2*a+1)-2*a^2*mu^2+(-2*a^2-8*a-8)*mu+a^2+2*a+2)+a^2*(sqrt(a^4*mu^4+(2*a^4+8*a^3+8*a^2)*mu^3+(-a^4+2*a^3+18*a^2+32*a+16)*mu^2+(-2*a^3-10*a^2-16*a-8)*mu+a^2+2*a+1)+3)+(2*a^4+4*a^3)*mu^3+a^3)/(mu*(a*(8*sqrt(a^4*mu^4+(2*a^4+8*a^3+8*a^2)*mu^3+(-a^4+2*a^3+18*a^2+32*a+16)*mu^2+(-2*a^3-10*a^2-16*a-8)*mu+a^2+2*a+1)-8)+a^2*(4*sqrt(a^4*mu^4+(2*a^4+8*a^3+8*a^2)*mu^3+(-a^4+2*a^3+18*a^2+32*a+16)*mu^2+(-2*a^3-10*a^2-16*a-8)*mu+a^2+2*a+1)-12)-2*a^4-8*a^3)+a*(4*sqrt(a^4*mu^4+(2*a^4+8*a^3+8*a^2)*mu^3+(-a^4+2*a^3+18*a^2+32*a+16)*mu^2+(-2*a^3-10*a^2-16*a-8)*mu+a^2+2*a+1)+4)+a^2*(2*sqrt(a^4*mu^4+(2*a^4+8*a^3+8*a^2)*mu^3+(-a^4+2*a^3+18*a^2+32*a+16)*mu^2+(-2*a^3-10*a^2-16*a-8)*mu+a^2+2*a+1)+6)+(4*a^4+8*a^3)*mu^3+(6*a^4+28*a^3+48*a^2+32*a)*mu^2+2*a^3)
  second <- (-sqrt(a^2+4*a+4)*sqrt(-2*sqrt(a^4*mu^4+(2*a^4+8*a^3+8*a^2)*mu^3+(-a^4+2*a^3+18*a^2+32*a+16)*mu^2+(-2*a^3-10*a^2-16*a-8)*mu+a^2+2*a+1)-2*a^2*mu^2-(2*a^2+8*a+8)*mu+a^2+2*a+2)+a^2+4*a+4)/(2*a^2+8*a+8)

  c(first, second, 1-first, 1-second)
}


# Phase diagram for fully symmetric case
plot.phases.FS <- function(alphalimits = c(0.01, 100),
                           resolution = 1000) {
  df <- expand.grid(alpha=hipster::logseq(from=alphalimits[1], to=alphalimits[2], length.out=resolution),
                    critpoint=c("tD", "D", "D2"),
                    value=NA)

  for (i in 1:nrow(df)) {
    if (df[i,]$critpoint == "tD") {
      df[i,]$value <- (df[i,]$alpha + 1)/(3*df[i,]$alpha + 4)
    }
    if (df[i,]$critpoint == "D") {
      df[i,]$value <- 1/(df[i,]$alpha + 4)
    }
    if (df[i,]$critpoint == "D2") {
      a <- df[i,]$alpha
      df[i,]$value <- (sqrt(2*a^2 + 4*a + 4) - a - 2)/a^2
    }
  }

  # plot
  g <- ggplot(df, aes(x=value, y=alpha, group=critpoint)) + geom_path()
  #g <- g + scale_y_continuous(expand=c(0,0))
  g <- g + scale_y_log10(expand=c(0,0))
  g <- g + scale_x_reverse(limits=c(0.5, 0.0), expand=c(0,0))
  g <- g + theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.text=element_text(color="black"))
  g <- g + xlab(expression(mu)) + ylab("s")
  g <- g + annotation_logticks(sides="l", size=0.3)
  g <- g + xlab("") + ylab("")

  # add coloured polygons
  g <- g + geom_polygon(inherit.aes=FALSE, data=data.frame(x=c(0.5, 1/4, unique(df[df$critpoint=="tD", ]$value), 0.5), y=c(0, 0, unique(df$alpha), alphalimits[2])), aes(x=x, y=y), fill=mypal[1], alpha=0.3)
  g <- g + geom_polygon(inherit.aes=FALSE, data=data.frame(x=c(rev(unique(df[df$critpoint=="tD", ]$value)), 1/4, unique(df[df$critpoint=="D", ]$value), max(unique(df[df$critpoint=="tD", ]$value))), y=c(rev(unique(df$alpha)), 0, unique(df$alpha), alphalimits[2])), aes(x=x, y=y), fill=mypal[2], alpha=0.3)
  g <- g + geom_polygon(inherit.aes=FALSE, data=data.frame(x=c(1/4, unique(df[df$critpoint=="D2", ]$value), max(unique(df[df$critpoint=="D", ]$value)), rev(unique(df[df$critpoint=="D", ]$value)), 1/4), y=c(0, unique(df$alpha), alphalimits[2], rev(unique(df$alpha)), 0)), aes(x=x, y=y), fill=mypal[3], alpha=0.3)
  g <- g + geom_polygon(inherit.aes=FALSE, data=data.frame(x=c(1/4, unique(df[df$critpoint=="D2", ]$value), 0, 0, 1/4), y=c(0, unique(df$alpha), alphalimits[2], 0, 0)), aes(x=x, y=y), fill=mypal[4], alpha=0.3)

  # add phase labels
  #dfl <- data.frame(mu=c(0.4, 0.2, 0.15, 0.1), alpha=c(1.0, 10.0, 2.0, 0.1), text=c("I", "II", "III", "IV"), critpoint=0)
  #g <- g + geom_text(data=dfl, aes(x=mu, y=alpha, label=text), size=5)

  g
}


# Phase a phase space plot for FS case
plot.space.FS <- function(which,
                          alphas = c(1, 1, 1, 1),
                          sloperes = 16,
                          slopelen = 0.025,
                          nullclines = TRUE,
                          ncreso = 10000,
                          mus = c(0.35, 0.25, 0.18, 0.13)) {
  if (which == "I") {
    palwhich <- 1
  } else if (which == "II") {
    palwhich <- 2
  } else if (which == "III") {
    palwhich <- 3
  } else if (which == "IV") {
    palwhich <- 4
  }


  # construct data
  dfI <- data.frame(x=0.5, y=0.5, stability="sink")
  dfI$phase <- "I"

  dfII <- data.frame(x=c(0.5, tDeq(alpha=alphas[2], mu=mus[2])[1], tDeq(alpha=alphas[2], mu=mus[2])[2]), y=c(0.5, 1-tDeq(alpha=alphas[2], mu=mus[2])[1], 1-tDeq(alpha=alphas[2], mu=mus[2])[2]), stability=c("saddle", "sink", "sink"))
  dfII$phase <- "II"

  dfIII <- data.frame(x=0.5, y=0.5, stability="source")
  dfIII <- rbind(dfIII, data.frame(x=tDeq(alpha=alphas[3], mu=mus[3])[1], y=1-tDeq(alpha=alphas[3], mu=mus[3])[1], stability="sink"))
  dfIII <- rbind(dfIII, data.frame(x=tDeq(alpha=alphas[3], mu=mus[3])[2], y=1-tDeq(alpha=alphas[3], mu=mus[3])[2], stability="sink"))
  dfIII <- rbind(dfIII, data.frame(x=Deq(alpha=alphas[3], mu=mus[3])[1], y=Deq(alpha=alphas[3], mu=mus[3])[1], stability="saddle"))
  dfIII <- rbind(dfIII, data.frame(x=Deq(alpha=alphas[3], mu=mus[3])[2], y=Deq(alpha=alphas[3], mu=mus[3])[2], stability="saddle"))
  dfIII$phase <- "III"

  dfIV <- data.frame(x=0.5, y=0.5, stability="source")
  dfIV <- rbind(dfIV, data.frame(x=tDeq(alpha=alphas[4], mu=mus[4])[1], y=1-tDeq(alpha=alphas[4], mu=mus[4])[1], stability="sink"))
  dfIV <- rbind(dfIV, data.frame(x=tDeq(alpha=alphas[4], mu=mus[4])[2], y=1-tDeq(alpha=alphas[4], mu=mus[4])[2], stability="sink"))
  dfIV <- rbind(dfIV, data.frame(x=Deq(alpha=alphas[4], mu=mus[4])[1], y=Deq(alpha=alphas[4], mu=mus[4])[1], stability="sink"))
  dfIV <- rbind(dfIV, data.frame(x=Deq(alpha=alphas[4], mu=mus[4])[2], y=Deq(alpha=alphas[4], mu=mus[4])[2], stability="sink"))
  dfIV <- rbind(dfIV, data.frame(x=offeq(alpha=alphas[4], mu=mus[4])[1], y=offeq(alpha=alphas[4], mu=mus[4])[2], stability="saddle"))
  dfIV <- rbind(dfIV, data.frame(x=1-offeq(alpha=alphas[4], mu=mus[4])[1], y=1-offeq(alpha=alphas[4], mu=mus[4])[2], stability="saddle"))
  #dfIV <- rbind(dfIV, data.frame(x=offeq(alpha=alpha, mu=mus[4])[1], y=1-offeq(alpha=alpha, mu=mus[4])[2], stability="saddle"))
  #dfIV <- rbind(dfIV, data.frame(x=1-offeq(alpha=alpha, mu=mus[4])[1], y=offeq(alpha=alpha, mu=mus[4])[2], stability="saddle"))
  #dfIV <- rbind(dfIV, data.frame(x=offeq(alpha=alpha, mu=mus[4])[2], y=offeq(alpha=alpha, mu=mus[4])[3], stability="saddle"))
  dfIV <- rbind(dfIV, data.frame(x=1-offeq(alpha=alphas[4], mu=mus[4])[2], y=offeq(alpha=alphas[4], mu=mus[4])[3], stability="saddle"))
  #dfIV <- rbind(dfIV, data.frame(x=1-offeq(alpha=alpha, mu=mus[4])[2], y=1-offeq(alpha=alpha, mu=mus[4])[3], stability="saddle"))
  dfIV <- rbind(dfIV, data.frame(x=offeq(alpha=alphas[4], mu=mus[4])[2], y=1-offeq(alpha=alphas[4], mu=mus[4])[3], stability="saddle"))
  dfIV$phase <- "IV"

  df <- rbind(dfI, dfII, dfIII, dfIV)
  df <- df[df$phase == which, ]

  g <- ggplot(df, aes(x=x, y=y)) 
  g <- g + theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.text=element_text(color="black"), strip.background=element_blank(), legend.position="none")
  g <- g + xlab("") + ylab("")
  g <- g + scale_x_continuous(breaks=c(0,0.5,1), labels=c("0", "0.5", "1"), limits=c(0,1), expand=c(0,0)) + scale_y_continuous(breaks=c(0,0.5,1), labels=c("0", "0.5", "1"), limits=c(0,1), expand=c(0,0)) + scale_shape_manual(values=c(19, 10, 1), name="")

  # add coloured polygons
  g <- g + geom_polygon(inherit.aes=FALSE, data=data.frame(phase="I", x=c(0, 1, 1, 0), y=c(0, 0, 1, 1)), aes(x=x, y=y), fill=mypal[palwhich], alpha=0.3)
  #g <- g + geom_polygon(inherit.aes=FALSE, data=data.frame(phase="II", x=c(0, 1, 1, 0), y=c(0, 0, 1, 1)), aes(x=x, y=y), fill=mypal[2], alpha=0.3)
  #g <- g + geom_polygon(inherit.aes=FALSE, data=data.frame(phase="III", x=c(0, 1, 1, 0), y=c(0, 0, 1, 1)), aes(x=x, y=y), fill=mypal[3], alpha=0.3)
  #g <- g + geom_polygon(inherit.aes=FALSE, data=data.frame(phase="IV", x=c(0, 1, 1, 0), y=c(0, 0, 1, 1)), aes(x=x, y=y), fill=mypal[4], alpha=0.3)

  # slope field
  slopes <- slope.field.2D(FUN=verger, resolution=sloperes, len=slopelen, S=S2D(s12=alphas[palwhich], s21=alphas[palwhich]), m=mus[palwhich], tm=mus[palwhich], normalize=TRUE)
  g <- g + geom_segment(data=slopes, aes(x=x1, y=x2, xend=x1end, yend=x2end), arrow=arrow(angle=30, length=unit(0.07, "cm")), alpha=0.6, lwd=0.2)

  if (nullclines) {
    # nullclines
    nullcline <- function(x, alpha, mu) {
      x1 <- x
      a1 <- alpha
      m1 <- mu
      M1 <- 2*mu + 1
      ((2+a1)*x1^3 - ((1-m1)*a1 + 3)*x1^2 + M1*x1 - m1)/(M1*a1*x1^2 - M1*a1*x1 + m1*a1)
    }

    ncx <- seq(from=0, to=1, length.out=ncreso)
    ncI <- data.frame(x=ncx, y=nullcline(x=ncx, alpha=alphas[1], mu=mus[1]), phase="I")
    ncII <- data.frame(x=ncx, y=nullcline(x=ncx, alpha=alphas[2], mu=mus[2]), phase="II")
    ncIII <- data.frame(x=ncx, y=nullcline(x=ncx, alpha=alphas[3], mu=mus[3]), phase="III")
    ncIV <- data.frame(x=ncx, y=nullcline(x=ncx, alpha=alphas[4], mu=mus[4]), phase="IV")

    nc <- rbind(ncI, ncII, ncIII, ncIV)
    nc <- nc[nc$phase == which, ]

    g <- g + geom_path(data=nc, inherit.aes=FALSE, aes(x=x, y=y), alpha=0.5)
    g <- g + geom_path(data=nc, inherit.aes=FALSE, aes(x=y, y=x), alpha=0.5)
  }

  # add equilibria
  g <- g + geom_point(size=2.5, shape=21, bg="white")
  g <- g + geom_point(aes(shape=stability), size=2.5)

  g
}


# Bifurcation diagram for fixed sigma ("tulip") for fully symmetric case
plot.bifu.FS <- function(alpha = 1,
                         mulimits = c(0, 0.5),
                         resolution = 5000) {
  # critical values of mu
  mu_tD <- (alpha + 1)/(3*alpha + 4)
  mu_D <- 1/(alpha + 4)
  mu_D2 <- (sqrt(2*alpha^2 + 4*alpha + 4) - alpha - 2)/alpha^2

  # initialize dataframe to hold all information
  df <- expand.grid(mu=seq(from=mulimits[1], to=mulimits[2], length.out=resolution),
                    equilibrium=c("u", "d1", "d2", "td1", "td2", "od11", "od12", "od21", "od22"),
                    stability=c("sink", "saddle", "source"),
                    value=NA)

  # fill dataframe
  for (i in 1:nrow(df)) {
    # uniform state u
    if (df[i,]$equilibrium == "u") {
      if (df[i,]$stability == "sink") {
        if (df[i,]$mu > mu_tD) {
          df[i,]$value <- 0.5
        }
      }
      if (df[i,]$stability == "saddle") {
        if (df[i,]$mu < mu_tD & df[i,]$mu > mu_D) {
          df[i,]$value <- 0.5
        }
      }
      if (df[i,]$stability == "source") {
        if (df[i,]$mu < mu_D) {
          df[i,]$value <- 0.5
        }
      }
    }

    # tD-equilibria
    if (df[i,]$equilibrium == "td1") {
      if (df[i,]$stability == "sink") {
        if (df[i,]$mu < mu_tD) {
          df[i,]$value <- tDeq(alpha=alpha, mu=df[i,]$mu)[1]
        }
      }
    }
    if (df[i,]$equilibrium == "td2") {
      if (df[i,]$stability == "sink") {
        if (df[i,]$mu < mu_tD) {
          df[i,]$value <- tDeq(alpha=alpha, mu=df[i,]$mu)[2]
        }
      }
    }

    # D-equilibria
    if (df[i,]$equilibrium == "d1") {
      if (df[i,]$stability == "saddle") {
        if (df[i,]$mu < mu_D & df[i,]$mu > mu_D2) {
          df[i,]$value <- Deq(alpha=alpha, mu=df[i,]$mu)[1]
        }
      }
      if (df[i,]$stability == "sink") {
        if (df[i,]$mu < mu_D2 & df[i,]$mu > 0) {
          df[i,]$value <- Deq(alpha=alpha, mu=df[i,]$mu)[1]
        }
      }
    }
    if (df[i,]$equilibrium == "d2") {
      if (df[i,]$stability == "saddle") {
        if (df[i,]$mu < mu_D & df[i,]$mu > mu_D2) {
          df[i,]$value <- Deq(alpha=alpha, mu=df[i,]$mu)[2]
        }
      }
      if (df[i,]$stability == "sink") {
        if (df[i,]$mu < mu_D2 & df[i,]$mu > 0) {
          df[i,]$value <- Deq(alpha=alpha, mu=df[i,]$mu)[2]
        }
      }
    }

    # off-diagonal equilibria
    if (df[i,]$equilibrium == "od11") {
      if (df[i,]$stability == "saddle") {
        if (df[i,]$mu < mu_D2 & df[i,]$mu > 0) {
          df[i,]$value <- offeq(alpha=alpha, mu=df[i,]$mu)[1]
        }
      }
    }
    if (df[i,]$equilibrium == "od12") {
      if (df[i,]$stability == "saddle") {
        if (df[i,]$mu < mu_D2 & df[i,]$mu > 0) {
          df[i,]$value <- offeq(alpha=alpha, mu=df[i,]$mu)[2]
        }
      }
    }
    if (df[i,]$equilibrium == "od21") {
      if (df[i,]$stability == "saddle") {
        if (df[i,]$mu < mu_D2 & df[i,]$mu > 0) {
          df[i,]$value <- offeq(alpha=alpha, mu=df[i,]$mu)[3]
        }
      }
    }
    if (df[i,]$equilibrium == "od22") {
      if (df[i,]$stability == "saddle") {
        if (df[i,]$mu < mu_D2 & df[i,]$mu > 0) {
          df[i,]$value <- offeq(alpha=alpha, mu=df[i,]$mu)[4]
        }
      }
    }
  }

  # construct plot
  g <- ggplot(df, aes(x=mu, y=value, lty=stability)) + geom_path()
  g <- g + scale_x_reverse(expand=c(0,0))
  g <- g + scale_y_continuous(expand=c(0,0), limits=c(0,1))
  g <- g + scale_linetype_manual(values=c("solid", "dashed", "dotted"), name="")
  g <- g + theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.text=element_text(color="black"))
  #g <- g + xlab(expression(mu)) + ylab(expression(italic(x)*", "*italic(y)))
  #g <- g + xlab(expression(mu)) + ylab("x, y")
  g <- g + xlab("") + ylab("")

  # add polygons to indicate phases
  g <- g + geom_polygon(inherit.aes=FALSE, data=data.frame(x=c(0.5, mu_tD, mu_tD, 0.5), y=c(0, 0, 1, 1)), aes(x=x, y=y), fill=mypal[1], alpha=0.3)
  g <- g + geom_polygon(inherit.aes=FALSE, data=data.frame(x=c(mu_tD, mu_D, mu_D, mu_tD), y=c(0, 0, 1, 1)), aes(x=x, y=y), fill=mypal[2], alpha=0.3)
  g <- g + geom_polygon(inherit.aes=FALSE, data=data.frame(x=c(mu_D, mu_D2, mu_D2, mu_D), y=c(0, 0, 1, 1)), aes(x=x, y=y), fill=mypal[3], alpha=0.3)
  g <- g + geom_polygon(inherit.aes=FALSE, data=data.frame(x=c(mu_D2, 0, 0, mu_D2), y=c(0, 0, 1, 1)), aes(x=x, y=y), fill=mypal[4], alpha=0.3)
  g
}


# Quantitative empirical evaluation plot
empirical_eval <- function(empirical_delta = 0.08,
                           mus = c(0.1, 0.2, 0.227, 0.25, 0.3),
                           sigmas = exp(seq(from=log(0.01), to=log(100), length.out=1000))) {
  df <- expand.grid(mu=mus, sigma=sigmas, delta=0)
  df$delta <- sqrt(-(3*df$sigma^2 + 4*df$sigma)*df$mu^2 - 2*(df$sigma^2 + 3*df$sigma + 2)*df$mu + (df$sigma + 1)^2)/(2*(df$sigma*df$mu + df$sigma + 1))

  g <- ggplot(df, aes(x=sigma, y=delta, group=mu, color=mu)) + geom_line()
  g <- g + theme_bw() + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), axis.text=element_text(color="black"), legend.position="none")
  g <- g + scale_x_log10(labels=fancy_scientific) + annotation_logticks(sides="bt", size=0.3)
  g <- g + ylim(0, 0.5)
  g <- g + scale_color_distiller(palette="Set1")
  g <- g + geom_hline(yintercept=empirical_delta, lty=2)
  g <- g + xlab("") + ylab("")

  g
}


# Quantitative empirical evaluation plot; faceted version
empirical_eval_facet <- function(empirical_delta = 0.08,
                           mus = c(0.2, 0.23, 0.25, 0.3),
                           sigmas = exp(seq(from=log(0.01), to=log(100), length.out=1000))) {
  df <- expand.grid(mu=mus, sigma=sigmas, delta=0)
  
  get_delta <- function(sigma, mu) {
    df$delta <- sqrt(-(3*sigma^2 + 4*sigma)*mu^2 - 2*(sigma^2 + 3*sigma + 2)*mu + (sigma + 1)^2)/(2*(sigma*mu + sigma + 1))
  }

  df$delta <- get_delta(df$sigma, df$mu)

  g <- ggplot(df, aes(x=sigma, y=delta, group=mu, color=mu)) + geom_line()
  g <- g + theme_bw() + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), axis.text=element_text(color="black"), legend.position="none")
  g <- g + scale_x_log10(labels=fancy_scientific, expand=c(0,0)) + annotation_logticks(sides="b")
  g <- g + scale_y_continuous(expand=c(0,0), limits=c(0, 0.4))
  g <- g + scale_color_distiller(palette="Dark2")
  g <- g + geom_hline(yintercept=empirical_delta, lty=2)
  g <- g + xlab("") + ylab("")

  # add envelopes
  if (1==0) {
  g <- g + geom_hline(yintercept=get_delta(c(0.00001, 100000), 0.1), color=evalpal[1], lty=2, alpha=0.5)
  g <- g + geom_hline(yintercept=get_delta(c(0.00001, 100000), 0.2), color=evalpal[2], lty=2, alpha=0.5)
  g <- g + geom_hline(yintercept=get_delta(c(0.00001, 100000), 0.23), color=evalpal[3], lty=2, alpha=0.5)
  g <- g + geom_hline(yintercept=get_delta(c(0.00001, 100000), 0.25), color=evalpal[4], lty=2, alpha=0.5)
  g <- g + geom_hline(yintercept=get_delta(c(0.00001, 100000), 0.3), color=evalpal[5], lty=2, alpha=0.5)
  }

  polydf <- NULL
  for (mu in mus) {
    polydf <- rbind(polydf, data.frame(sigma=0, delta=0, mu=mu, x=c(min(sigmas), max(sigmas), max(sigmas), min(sigmas)), y=c(get_delta(0.00001, mu), get_delta(0.00001, mu), get_delta(100000, mu), get_delta(100000, mu))))
  }
  polydf$y <- ifelse(is.na(polydf$y), 0.0, polydf$y)
  g <- g + geom_polygon(data=polydf, aes(x=x, y=y, fill=mu), color=NA, lty=2, alpha=0.1) + facet_wrap(.~mu, nrow=1) + scale_fill_distiller(palette="Dark2")

  if (1==0) {
  g <- g + geom_polygon(inherit.aes=FALSE, data=data.frame(x=c(min(sigmas), max(sigmas), max(sigmas), min(sigmas)), y=c(get_delta(0.00001, 0.1), get_delta(0.00001, 0.1), get_delta(100000, 0.1), get_delta(100000, 0.1))), aes(x=x, y=y), fill=evalpal[5], alpha=0.2)
  g <- g + geom_polygon(inherit.aes=FALSE, data=data.frame(x=c(min(sigmas), max(sigmas), max(sigmas), min(sigmas)), y=c(get_delta(0.00001, 0.2), get_delta(0.00001, 0.2), get_delta(100000, 0.2), get_delta(100000, 0.2))), aes(x=x, y=y), fill=evalpal[4], alpha=0.2)
  g <- g + geom_polygon(inherit.aes=FALSE, data=data.frame(x=c(min(sigmas), max(sigmas), max(sigmas), min(sigmas)), y=c(get_delta(0.00001, 0.23), get_delta(0.00001, 0.23), get_delta(100000, 0.23), get_delta(100000, 0.23))), aes(x=x, y=y), fill=evalpal[3], alpha=0.2)
  g <- g + geom_polygon(inherit.aes=FALSE, data=data.frame(x=c(min(sigmas), max(sigmas), max(sigmas), min(sigmas)), y=c(get_delta(0.00001, 0.25), get_delta(0.00001, 0.25), get_delta(100000, 0.25), get_delta(100000, 0.25))), aes(x=x, y=y), fill=evalpal[2], alpha=0.2)
  g <- g + geom_polygon(inherit.aes=FALSE, data=data.frame(x=c(min(sigmas), max(sigmas), max(sigmas), min(sigmas)), y=c(get_delta(0.00001, 0.3), get_delta(0.00001, 0.3), get_delta(100000, 0.3), get_delta(100000, 0.3))), aes(x=x, y=y), fill=evalpal[1], alpha=0.2)
  }

  g
}
