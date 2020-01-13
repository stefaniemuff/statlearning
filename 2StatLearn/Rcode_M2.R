### Slide 12 ###
# Scatterplot for three different covariance structures of a bivariate model
library(mvtnorm)

sample1 <- rmvnorm(1000,c(0,0),sigma=matrix(c(1,-0.8,-0.8,1),ncol=2))
plot(sample1)

sample2 <- rmvnorm(1000,c(0,0),sigma=matrix(c(1,0,0,1),ncol=2))
plot(sample2)

sample3 <- rmvnorm(1000,c(0,0),sigma=matrix(c(1,0.8,0.8,1),ncol=2))
plot(sample3)

cov(sample3)