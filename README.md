# Ch-6

############################################################################

#   R code for in Chapter 6: Bivariate Normal Distribution

############################################################################

setwd("C:/Users/zeldan/Desktop/Projects/Multivariate/Data")

###########################################################################

# Figure 6.1: Swiss agriculture and high exam scores

library(datasets)
library(MVA)

biv <- swiss[,2:3]            #   Extract bivariate data

### pdf(file = "swiss.pdf")
 op <- par(mfrow=c(2, 3), mar = c(4, 4, 4, 0.25), cex.lab = 1.75)  

plot(biv[, 1], biv[, 2], pch = 16,        #  Naive scatterplot
     cex = 1.5, col = "blue", cex.lab = 1.75,
     xlab = "% agriculture", 
     ylab = "% with high exam")

bvbox(biv, xlab = "% agriculture", yaxt = "n", ylab = "", cex.lab = 1.75,
      pch = 16, cex = 1.5, col = "blue", method = "NOT")

library(aplpack)
nlev <- 5
colors <- heat.colors(9)[3 : (nlev + 2)]
plothulls(biv, n.hull = 5, col.hull = colors, cex.lab = 1.75, 
          xlab = "% agriculture", ylab = " ", yaxt = "n", 
          lty.hull = 1 : nlev, density = NA, col = 0, main = " ")
points(biv, pch = 16, cex = 1.5, col = "blue")

 par(op)
### dev.off()


############################################################################


#   Figure 6.2:  Bivariate normal contours

corcon <- function(x, y, correl)
{
  nx <- length(x)
  ny <- length(y)
  z <- matrix(rep(0, nx * ny), nx, ny)
  for (i in 1 : nx)
  {
    for (j in 1 : ny)
    {
      z[i, j] <- dmvnorm(c(x[i], y[j]), c(0, 0),
                        matrix(c(1, correl, correl, 1), 2, 2)) 
    }
  }
  z
}


library(datasets)
library(MVA)
library(mvtnorm)

del <- .05      # how fine the grid
lim <- 3.25     # std normals plotted on +/- lim

### pdf(file = "bivnormal.pdf")
 par(mfrow = c(2, 4), mar = c(5, 0, 5, 0), cex.lab = 1.5)   # Four plots across
contour(corcon(seq(-lim, lim, del), seq(-lim, lim, del), 
        correl = -.5), col = "blue",
        xlab = expression(rho == -.5), cex.lab = 1.5,
        drawlabels = FALSE, axes = FALSE, frame = TRUE)
contour(corcon(seq(-lim, lim, del), seq(-lim, lim, del), 0),
        xlab = expression(rho == 0), cex.lab=1.5, col = "blue",
        drawlabels = FALSE, axes = FALSE, frame = TRUE)
contour(corcon(seq(-lim, lim, del), seq(-lim, lim, del), .5),
        xlab = expression( rho == .5), cex.lab = 1.5, col = "blue",
        drawlabels = FALSE, axes = FALSE, frame = TRUE)
contour(corcon(seq(-lim, lim, del), seq(-lim, lim, del), .9),
        xlab = expression(rho == .9), cex.lab = 1.5, col = "blue",
        drawlabels = FALSE, axes = FALSE, frame = TRUE)
 par(op)
### dev.off()

############################################################################

#  Figure 6.3: Heat map and profile

### pdf(file = "bivnormal2.pdf")
library(MASS)
library(mvtnorm)
library(graphics)

 layout(t(matrix(c(1 : 2, rep(0, 2)), 2, 2)), widths=c(1, 1))

del <- .0125      # how fine the grid
lim <- 1.25      # std normals plotted on +/- lim

image(corcon(seq( -lim, lim, del), seq(-lim, lim, del), 
    correl = 0.8), axes = FALSE) 

del <- .1        # how fine the grid
lim <- 2.7       # std normals plotted on +/- lim

persp(corcon(seq(-lim, lim, del), seq( -lim, lim, del), 
      correl = 0.8),
      axes = FALSE, xlab = " ", ylab = " ", box = FALSE, 
      col = "lightblue", shade = .05)

par(op)

### dev.off()

############################################################################

# Figure 6.4 :  Boxplots for conditional distributions

ag <- swiss[, 2]
ex <- swiss[, 3]
low <- min(ex) - 1                     #  set Y ranges
hi <-  max(ex) + 5
xl <- quantile(ag, probs = c(.25, .5, .75))

### pdf(file = "condbox.pdf")
### op <- par(mfrow=c(2, 2), mar = c(5, 4, 5, 0)) 

plot(ag, ex, pch = 16, cex = 1.25, xaxt = "n", yaxt = "n", col = "red",
     xla b= "Agriculture", ylab = "Examination scores")    
lines(c(xl[1], xl[1]), c(low ,hi), lty = "dashed", col = "blue", lwd = 2)
lines(c(xl[2], xl[2]), c(low, hi), lty = "dashed", col = "blue", lwd = 2)
lines(c(xl[3], xl[3]), c(low, hi), lty = "dashed", col = "blue", lwd = 2)

quart <- rep(1, length(ag))    #  Form quartile categories
for (j in 1 : 3)   quart <- quart + (ag >= xl[j])

text( c(0, 1), c(0, 0), 
      boxplot(ex ~ quart, col = "blue", axes = FALSE, 
              boxwex = .25, names = c("1Q", "2Q", "3Q", "4Q"),
              ylim = range(ex), xlab = "Agriculture quartiles"))

### par(op)
###dev.off()

############################################################################

# Figure 6.5:  Bivariate confidence ellipse

library(datasets)
library(MASS)
library(MVA)

biv <- swiss[, 2 : 3]          #   Extract bivariate data

bivCI <- function(s, xbar, n, alpha, m)
  #  returns m (x,y) coordinates of 1-alpha joint confidence ellipse of mean
  
{
  x <- sin( 2* pi * (0 : (m - 1) )/ (m - 1))    # m points on a unit circle
  y <- cos( 2* pi * (0 : (m - 1)) / (m - 1))
  cv <-  qchisq(1 - alpha, 2)                   # chisquared critical value
  cv <- cv / n                                  # value of quadratic form
  for (i in 1 : m)
  {
    pair <- c(x[i], y[i])        # ith (x,y) pair
    q <- pair %*% solve(s, pair)  # quadratic form
    x[i] <- x[i] * sqrt(cv / q) + xbar[1]
    y[i] <- y[i] * sqrt(cv / q) + xbar[2]
  }
  return(cbind(x, y))
}

### pdf(file = "bivSwiss.pdf")

plot(biv, col = "red", pch = 16, cex.lab = 1.5)
lines(bivCI(var(biv), colMeans(biv), dim(biv)[1], .01, 1000), type = "l",
      col = "blue")
lines(bivCI(var(biv), colMeans(biv), dim(biv)[1], .05, 1000), 
      type = "l", col = "green", lwd = 1)
lines(colMeans(biv)[1], colMeans(biv)[2], pch = 3, cex = .8, type = "p", 
     lwd = 1)

### dev.off()


##########################################################################

# Section 6.5: Maximum Likelihood Estimation II

##########################################################################


# US all-cancer rates for males in 2007 and 2008 by state + DC

cancer <- read.table(file = "Cancer2007_8.dat", header = TRUE, row.names = 1)
cancer[1 : 3, ]       # first few
cancer[49 : 51, ]     # last few

# Simple summary statistics

colMeans(cancer)
sapply(cancer, sd)
var(cancer)
cor(cancer)

####################################################################

#   Figure 6.6:  Plot data and identify extremes

### pdf(file = "cancer2007.pdf") 

plot(cancer[, 1], cancer[, 2], xlab = "2007 Rate", 
     ylab = "2008 Rate", pch = 19, col= "red", cex = 1.5)
ext <- c(3, 9, 29)      # identified extremes
text(cancer[ext, 1], cancer[ext, 2], 
     labels = row.names(cancer)[ext], pos = c(3, 2, 3), col = "blue")

rc <- glm(cancer[, 2] ~ cancer[, 1])$coefficients # regression coefficients
xb <- range(cancer[,1]) * c( .9, 1.1)
yb <- rc[1] + xb * rc[2]
lines(xb, yb, col = "green")
text(538, 496, labels = "regression line", col = "green")

lines(xb, xb, col = "blue")  # diagonal line
text(580, 610, labels = "equal rates line", col = "blue")

### dev.off()

####################################################################
#  Fit some models:

library(mvtnorm)       # library with dmvnorm function
biv5 <- function(par)  # all five parameter, natural parameteriztions
{
  cov <- par[5]* sqrt(par[3] * par[4]) 
  biv5 <- sum(
    -dmvnorm(cancer, mean = c(par[1], par[2]),
             sigma = matrix(c(par[3], cov, cov, par[4]),2, 2), 
             log = TRUE) )
  print(c(par, biv5))
  return(biv5)
}


nlm(biv5, c(45, 45, 1600, 1600, .8), hessian = TRUE)  # fits with warnings

nlm.out <- nlm(biv5, c(45, 45, 1600, 1600, .8), hessian = TRUE)

nlm.out$estimate                    # parameter estimates

nlm.out$hessian                     # estimated Hessian matrix

solve(nlm.out$hessian)              # invert the Hessian matrix

diag(solve(nlm.out$hessian))        # diagonal elements

sqrt(diag(solve(nlm.out$hessian)))  # estimated se"s 

se <- sqrt(diag(solve(nlm.out$hessian))) 
q <- qnorm(.975)                      # 95% interval quantile
upper <- nlm.out$estimate + q * se    # upper 95% CI intervals
lower <- nlm.out$estimate - q * se    # lower 95% CI intervals
summary <- data.frame(              
  cbind(nlm.out$estimate, se, lower, upper),
  row.names = c("mean1", "mean2", "var1", "var2", "rho"))
colnames(summary)[1] <- "estimate"
print(summary, digits = 3)


##########################################################################

# Exercise: Marginal normality is not bivariate normality

checker <- function(n)    # generates n pairs of marginal normals
  # that are nor bivariate normal
{
  for (i in 1 : n)
  {
    x <- rnorm(2)      # pair of independent normals
    if(x[1] > 0)  x[2] <-  abs(x[2])
       else       x[2] <- -abs(x[2])
    if(i == 1) checker <- x
       else checker <- rbind(checker, x)
  }
  return(checker)
}


# solution to uncorrlated pairs:

rad <- 1.83  # found by trial and error 

quad <- function(n)
{   
  for (i in 1 : n)
  {
    x <- rnorm(1)
    y <- rnorm(1)
    dia <- sqrt(x ^ 2 + y ^ 2)
    
    if(dia > rad)    # outside the circle
    {
      if(x>0) y <-  abs(y)
      else y <- -abs(y)
      if(x<0) y <- -abs(y)
      else y <-  abs(y)
    }
    
    if(dia < rad)    # inside the circle
    {
      if(x > 0) y <- -abs(y)
      else y <-  abs(y)
      if(x<0)y <-  abs(y)
      else y <- -abs(y)
    }
    
    if(i == 1) quad <- c(x, y)
    else quad <- rbind(quad, c(x, y))
  }
  
  return(quad)
}

### pdf(file = "quad.pdf")

require(graphics)

xy <- quad(750)                   # Bivariate data sample
boundx <- range(xy[, 1]) * 1.05
boundy <- range(xy[, 2]) * 1.05
plot(xy, xlim = boundx, ylim = boundy, # Plot data with these bounds
     axes = FALSE, cex = .7, pch = 16, # Omit axes, large bold circles
     col = "blue", xlab = "X", yla b= "Y")
rug(xy[, 1], side = 1)                 # Add rug fringes instead of axes
rug(xy[, 2], side = 2)

circ <- (0 : 200) * 2 * pi / 200           # add circle
circ <- cbind( sin(circ), cos(circ)) * rad
lines(circ, col = "red", type = "l")

### dev.off()


###########################################################################


#  Exercise 6.5: Reparameterize to avoid warnings

biv5r <- function(par)  # all five parameter, reparameterized
{
  sig1 <- exp(par[3])
  sig2 <- exp(par[4])
  rho <- par[5] / sqrt(1 + par[5] ^ 2)
  cov <- rho * sig1 * sig2 
  biv5 <- sum(
    -dmvnorm(cancer, mean = c(par[1], par[2]),
             sigma = matrix(c(sig1 ^ 2, cov, cov, sig2 ^ 2), 2, 2), 
             log = TRUE) )
  print(c(par[1 : 2], sig1, sig2, rho, biv5))
  return(biv5)
}

nlm(biv5r, c(45, 45, 7.25, 7.25, 2))  # fits without warnings



############################################################################
#####################  end of this file ####################################
############################################################################

