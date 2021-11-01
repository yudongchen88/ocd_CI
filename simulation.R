# this file demonstrates how the ocd_CI algorithm can be used to compute a 
# confidence interval of the changepoint location after declaring a change

library(putils)  # install using install_github('wangtengyao/putils')
# problem parameters
p <- 100; beta <- 1; vartheta <- 1   
z <- 500; s <- 10; theta <- vector.normalise(c(rnorm(s), rep(0, p-s))) * vartheta
# ocd parameters
T_diag <- 13.43; T_sparse <- 61.74   
# # ocd_CI parameters
alpha <- 0.05;                       
d1 <- 1/2*sqrt(log(p/alpha)); d2 <- d1^2 * 4; l <- 0   

set.seed(1)
# initiliase
n <- 0
bunch(A, tail) %=% initialise(p)

# monitoring for changepoint
repeat{
  n <- n + 1
  x <- rnorm(p) + theta * (n > z)
  bunch(stat, A, tail) %=% process_new_obs(x, A, tail, beta)
  if (stat['diag'] > T_diag || stat['sparse'] > T_sparse) break
}
N <- n  # time of declaration
println('Declaration at time: ', N)

# possible post-declaration additional sample
for (n in N + seq_len(l)) {
  x <- rnorm(p) + theta * (n > z)
  bunch(stat, A, tail) %=% process_new_obs(x, A, tail, beta)
}

# compute anchor tail length and anchor coordinate
bunch(anchor_tail_length, anchor_coord) %=% find_anchor(A, tail, beta)

# compute support estimate 
bunch(S_hat, b_tilde) %=% support_estimate(A, tail, anchor_coord, anchor_tail_length, beta, d1)

# compute confidence interval
CI <- conf_int(N, tail, S_hat, b_tilde, d2)
println('Confidence interval: [', paste(CI, collapse=', '), ']')
