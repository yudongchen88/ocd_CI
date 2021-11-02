# this file demonstrates how the ocd_CI algorithm can be applied to study 
# market movements leading up to the financial crisis of 2007–2008

library(putils)  # install using install_github('wangtengyao/putils')
source('algorithms.R')

# output SNP sectors
output_snp_sector <- function(test_std, supp_est){
  snp_sector <- read.csv('./snp500_sector.csv', stringsAsFactors = FALSE)
  snp_sector$Sector <- as.factor(snp_sector$Sector)
  snp_sector <- snp_sector[snp_sector$Symbol %in% colnames(test_std), c(1,4)]
  rownames(snp_sector) <- match(snp_sector$Symbol, colnames(test_std))
  
  snp_support_stocks <- snp_sector[supp_est, ]
  print(snp_support_stocks$Sector)
  sector_prop <- table(snp_support_stocks$Sector)/table(snp_sector$Sector)
  return(sector_prop)  
}

#######################################
# process S&P500 data
# read data
df <- read.csv('./snp500.csv', stringsAsFactors = FALSE)

# calculate log returns
all_dates <- df$X
df <- apply(log(df[, -1]), 2, diff, 1, 1)
all_dates <- all_dates[-1]

# train test split and standardise data
last_train_date='2006-12-31'; first_test_date='2007-01-01'
train_ind <- all_dates <= last_train_date
test_ind <- all_dates >= first_test_date
train_dates <- all_dates[train_ind]
test_dates <- all_dates[test_ind]
df_std <- t((t(df) - apply(df[train_ind, ], 2, mean))/apply(df[train_ind, ], 2, mad))
train_std <- df_std[train_ind, ]
test_std <- df_std[test_ind, ]

# clip standardised data at +/-3
test_std <- pmax(pmin(test_std, 3), -3)
rownames(test_std) <- test_dates

#---------------------------------------------------------------
# ocd_CI to be applied to the testing dataset 
# problem parameters
p <- ncol(test_std); beta <- 50; gamma <- 1000; alpha = 0.05
# ocd parameters
T_diag <- log(16*p*gamma*log2(4*p)); T_sparse <- 8*log(16*p*gamma*log2(2*p)) 
# # ocd_CI parameters (with l = 0)
alpha <- 0.05;                       
d1 <- 1/2*sqrt(log(p/alpha)); d2 <- d1^2 * 4

# Run 4 changepoints
for (i in 1:4){
  if (i == 1){
    first_test_date = '2007-01-01'
  } else {
    first_test_date = rownames(test_std)[CI[2] + 10] # 10 days of cooldown
  }
  test_std = test_std[rownames(test_std) >= first_test_date, ]

  # initiliase
  n <- 0
  bunch(A, tail) %=% initialise(p)
  
  # monitoring for changepoint
  repeat{
    n <- n + 1
    x <- as.numeric(test_std[n,])
    bunch(stat, A, tail) %=% process_new_obs(x, A, tail, beta)
    if (stat['diag'] > T_diag || stat['sparse'] > T_sparse) break
  }
  N <- n  # time of declaration
  println('Declaration on ', rownames(test_std)[N])
  
  # compute anchor tail length and anchor coordinate
  bunch(anchor_tail_length, anchor_coord) %=% find_anchor(A, tail, beta)
  
  # compute support estimate 
  bunch(S_hat, b_tilde) %=% support_estimate(A, tail, anchor_coord, anchor_tail_length, beta, d1)
  println('Estimated support of the change stocks: ', paste(colnames(test_std)[S_hat], collapse = ', '))
  
  # compute confidence interval
  CI <- conf_int(N, tail, S_hat, b_tilde, beta, d2)
  CI_dates <- c(as.Date(rownames(test_std)[CI[1]]), as.Date(rownames(test_std)[CI[2]]))
  println('Confidence interval: [', paste(CI_dates, collapse=', '), ']')
  
  # print out the sector distribution of the stocks in the support
  print(output_snp_sector(test_std,S_hat))
}

