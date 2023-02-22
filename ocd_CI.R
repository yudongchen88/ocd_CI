# This R file provides the main function of the ocd_CI algorithm, ocd_CI, which can be used to compute a 
# confidence interval of the changepoint location and provide a support estimate after declaring a change
library('MASS')
source('ingredients.R')
source('auxiliary.R')

ocd_CI <- function(p=100, s=10, z=1000, n_max=50000, vartheta=1, signal_shape=c('random', 'uniform', 'inv_sqrt', 'harmonic'), 
                   spatial=c('identity', 'toeplitz'), rho=0,
                   gamma=30000, beta, ocd_thresholds=c('MC', 'theoretical'),
                   l_choice=c('zero', 'given', 'other'), l=0,
                   d1_choice=c('default', 'variant1', 'variant2', 'other'), d1=0, alpha=0.05, kaul_comp=FALSE,
                   print_type = c('CI only', 'support only', 'both', 'muted')){
  
  # generate change vector
  theta <- vector_of_change(p, s, vartheta, signal_shape)

  # set d1, d2
  d1 <- 0.5*sqrt(log(p/alpha))  # default choice for d1
  if (d1_choice=='variant1') d1 <- sqrt(2*log(p/alpha))
  if (d1_choice=='variant2') d1 <- sqrt(2*log(p))
  if (d1_choice=='other') d1 <- d1
  d2 <- d1^2 * 4
  
  # set extra sample size l
  l <- 0  # defulat choice is zero
  if (l_choice=='given') l <- ceiling(2*s*log2(2*p)*log(p)*beta^{-2})
  if (l_choice=='other') l <- l
  
  # set ocd declaration thresholds
  if (ocd_thresholds=='theoretical'){
    T_diag <- log(16*p*gamma*log2(4*p))
    T_off <- 8*log(16*p*gamma*log2(2*p))
  } else {
    thresh_table <- read.csv('thresh.csv')
    param_str <- paste0(apply(thresh_table[,1:5], 1, toString))
    if (spatial=='identity'){
      spatial_flag <- 0
      rho <- 0
    } else {
      spatial_flag <- 1
    }
    cur_str <- toString(c(p, gamma, beta, spatial_flag, rho))
    T_diag <- thresh_table[cur_str == param_str, 6]
    T_off <- thresh_table[cur_str == param_str, 7]
  }
  
  # initiliase
  n <- 0
  bunch(A, tail) %=% initialise(p)
  
  if (kaul_comp) X <- matrix(,p,0)
  
  # monitoring for changepoint
  if (spatial=='identity'){
    repeat{
      n <- n + 1
      x <- rnorm(p) + theta * (n > z)
      if (kaul_comp) X <- cbind(X, x)
      bunch(stat, A, tail) %=% process_new_obs(x, A, tail, beta)
      if (stat['diag'] > T_diag || stat['off'] > T_off) break
    }
    N <- n  # time of declaration
  }
  if (spatial=='toeplitz'){
    repeat{
      n <- n + 1
      x <- mvrnorm(n=1, rep(0, p), Sigma=toeplitz(rho^(0:(p-1)))) + theta * (n > z)
      if (kaul_comp) X <- cbind(X, x)
      bunch(stat, A, tail) %=% process_new_obs(x, A, tail, beta)
      if (stat['diag'] > T_diag || stat['off'] > T_off) break
    }
    N <- n  # time of declaration
  }
  
  if (kaul_comp){
    kaul_CI <- kaul(X)
    if (print_type != 'muted'){
      println('Confidence interval based on Kaul et al. (2021): [', paste(kaul_CI, collapse=', '), ']')
    }
  }
  
  # possible post-declaration additional sample
  additional <- rep(0, p)
  for (n in N + seq_len(l)) {
    additional <- additional + rnorm(p) + theta * (n > z)
  }
  A <- A + additional
  tail <- tail + l
  colnames(A) <- paste0('tail=', sort(unique(as.vector(tail))))
  
  # compute anchor tail length and anchor coordinate
  bunch(anchor_tail_length, anchor_coord) %=% find_anchor(A, tail, beta)
  
  # compute support estimate
  # S_hat is the support (default, without the anchor coordinate)
  # S_hat_a is the support plus the anchor coordinate
  bunch(S_hat, b_tilde) %=% support_estimate(A, tail, anchor_coord, anchor_tail_length, beta, d1)
  S_hat <- S_hat[S_hat!=anchor_coord]
  S_hat_a <- sort(c(S_hat, anchor_coord))
  
  # compute confidence interval CI
  CI <- conf_int(N, tail, S_hat, b_tilde, beta, d2)
  
  # print statistical findings (support and confidence interval)
  if (print_type=='CI only'){
    println('Confidence interval ocd_CI: [', paste(CI, collapse=', '), ']')
  } else if (print_type=='support only'){
    println('Support coordinates (without anchor): ', paste(S_hat, collapse = ', '), '\n',
            'Support coordinates (with anchor): ', paste(S_hat_a, collapse = ', '))
  } else if (print_type=='both'){
    println('Confidence interval ocd_CI: [', paste(CI, collapse=', '), ']')
    println('Support coordinates (without anchor): ', paste(S_hat, collapse = ', '), '\n',
            'Support coordinates (with anchor): ', paste(S_hat_a, collapse = ', '))
  }
  
  if (kaul_comp){
    return(list(CI=CI, kaul_CI=kaul_CI, support=S_hat, support_aug=S_hat_a))
  } else {
    return(list(CI=CI, support=S_hat, support_aug=S_hat_a))
  }
}
  

# test example for each table
# Table 1
ocd_CI(p=100, s=10, vartheta=1, signal_shape='random', spatial='identity', beta=1,
       ocd_thresholds='MC', l_choice='zero', d1_choice='default', alpha=0.05, kaul_comp = TRUE,
       print_type ='both')
  
# Table 2
# ocd_CI(p=100, s=5, vartheta=1, signal_shape='harmonic', spatial='identity', beta=1,
#       ocd_thresholds='MC', l_choice='given', d1_choice='variant1', alpha=0.05, print_type='support only')

# Table S1
# ocd_CI(p=100, s=10, vartheta=1, signal_shape='random', spatial='toeplitz', rho=0.75, beta=1,
#       ocd_thresholds='MC', l_choice='zero', d1_choice='default', alpha=0.05, print_type='CI only')

# Figure 2 left panel
#ocd_CI(p=100, s=5, vartheta=2, signal_shape='uniform', spatial='identity', beta=2,
#       ocd_thresholds='MC', l_choice='zero', d1_choice='variant1', alpha=0.05, print_type ='support only')
# FIgure 2 right panel
#ocd_CI(p=100, s=20, vartheta=2, signal_shape='inv_sqrt', spatial='identity', beta=2,
#       ocd_thresholds='MC', l_choice='zero', d1_choice='variant2', alpha=0.05, print_type='support only')

