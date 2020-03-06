

set.seed(5143313)
n <- 200
beta <- c(-0.5,1)
var_e <- 0.25
var_x <- 2

x <- rnorm(n,0,sqrt(var_x))

eta <- beta[1] + beta[2]*x

y <- eta + rnorm(n,0,sqrt(var_e))

y_bin <- rbinom(n,1,exp(eta)/(1+exp(eta)))

#glm(y_bin ~ x, family="binomial")

y_pois <- rpois(n,exp(eta))

#glm(y_pois ~ x, family="poisson")
