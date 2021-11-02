# this file demonstrates how the ocd_CI algorithm can be applied to a dataset of
# weekly deaths number in the United States between January 2017 and June 2020

library(putils)  # install using install_github('wangtengyao/putils')
source('algorithms.R')

# process US weekly deaths data
# get the day numbers in a year for any given week
getDays <- function(end_date){
  days <- as.numeric(strftime(end_date, format = "%j"))
  days <- seq(-6, 0) + rep(days, each=7)
  days[days <= 0] <- days[days <= 0] + 365
  return(matrix(days, nrow=7))
}

# 1. estimate the average number of deaths on each day of the yearfrom these derived daily
#    death numbers using a Gaussian kernel with a bandwidth of 20 days
# 2. apply a square-root transformation to stabilise the variance
# 3. weekly excess deaths are computed on a square-root scale
process_state <- function(state, last_training_date='2019-06-30'){
  state_deaths <- US_deaths[US_deaths$region_code==state, 
                            c('end_date', 'total_deaths')]
  state_deaths$end_date <- as.Date(state_deaths$end_date)
  state_deaths_train <- state_deaths[state_deaths$end_date <= last_training_date, ]
  deaths <- rep(state_deaths_train$total_deaths / 7, each=7)
  days <- getDays(state_deaths_train$end_date)
  mean_deaths <- aggregate(deaths ~ as.vector(days), FUN=mean)
  bw <- 20
  dist <- outer(1:365, 1:365, 
                function(i,j){pmin(abs(i - j), i + 365 - j, j + 365 - i)})
  wt <- exp(-dist^2/(2*bw^2))
  wt <- wt / rowSums(wt)
  smooth_deaths <- wt %*% mean_deaths$deaths
  days <- getDays(state_deaths$end_date)
  sqrt_total_deaths <- sqrt(state_deaths$total_deaths) - sqrt(colSums(matrix(smooth_deaths[days], nrow=7)))
  ave <- mean(sqrt_total_deaths[state_deaths$end_date <= last_training_date])
  std <- sd(sqrt_total_deaths[state_deaths$end_date <= last_training_date])
  sqrt_total_deaths <- (sqrt_total_deaths - ave) / std
  return(sqrt_total_deaths)
}

# read data
US_deaths <- read.csv('./US_weekly_deaths.csv')

# use data up to the end of 2020
df <- data.frame(end_date=unique(US_deaths$end_date)[1:207])

# use data up to 2019-06-30 as training, rest as monitoring
last_training_date='2019-06-30'

# process the data in each state
for (state in unique(US_deaths$region_code)){
  df[[state]] <- process_state(state)[1:207]
}
df$end_date = as.Date(as.character(df$end_date))
rownames(df) = df$end_date

# construct testing data
first_test_week <- which(df$end_date>last_training_date)[1]
df_test = df[first_test_week:207, 2:52]

#---------------------------------------------------------------
# ocd_CI to be applied to the testing dataset 
# problem parameters
p <- ncol(df_test); beta <- 50; gamma <- 1000
# ocd parameters
T_diag <- log(16*p*gamma*log2(4*p)); T_sparse <- 8*log(16*p*gamma*log2(2*p)) 
# # ocd_CI parameters (with l = 0)
alpha <- 0.05;                       
d1 <- 1/2*sqrt(log(p/alpha)); d2 <- d1^2 * 4

# initiliase
n <- 0
bunch(A, tail) %=% initialise(p)

# monitoring for changepoint
repeat{
  n <- n + 1
  x <- as.numeric(df_test[n,])
  bunch(stat, A, tail) %=% process_new_obs(x, A, tail, beta)
  if (stat['diag'] > T_diag || stat['sparse'] > T_sparse) break
}
N <- n  # time of declaration
println('Declaration on week ending: ', rownames(df_test)[N])

# compute anchor tail length and anchor coordinate
bunch(anchor_tail_length, anchor_coord) %=% find_anchor(A, tail, beta)

# compute support estimate 
bunch(S_hat, b_tilde) %=% support_estimate(A, tail, anchor_coord, anchor_tail_length, beta, d1)
println('Estimated support of the change states: ', paste(colnames(df_test)[S_hat], collapse = ', '))

# compute confidence interval
CI <- conf_int(N, tail, S_hat, b_tilde, beta, d2)
CI_dates <- c(as.Date(rownames(df_test)[CI[1]-1])+1, as.Date(rownames(df_test)[CI[2]]))
println('Confidence interval: [', paste(CI_dates, collapse=', '), ']')


# Plot selected 12 states: figure 4 of Chen, Wang and Samworth (2021)
palet <- matplotlib_palette(10)
par(mfrow=c(1, 1), mar=c(3, 1.5, 0.2, 0.2), mgp=c(1.5,0.5,0), tcl=-0.03)
plot(c(as.Date(df$end_date[1])-6, as.Date(df$end_date[181])), c(0, 56), pch=' ', xlab='', ylab='', 
     xaxt = 'n', yaxt = 'n', xaxs = 'i', yaxs='i', frame.plot=F)
polygon(c(as.Date('2020-03-15'), as.Date('2020-03-15'), as.Date('2020-03-28'), as.Date('2020-03-28')), 
        c(-3, 59, 59, -3), col = "#00FFFF")
plot_states = c('FL', 'NC', 'TX', 'WY', 'CA', 'WA', 'DC', 'LA', 'MI', 'CT', 'NY', 'NJ')
for (j in 1:12){
  y <- df[1:181, c(plot_states[j])]
  y <- 12*(y - max(df[1:181,2:52])) / diff(range(df[1:181,2:52])) + 10 + 4*j
  if (j %in% 8:12){
    lines(df$end_date[1:181], y, col = 'red')
  }
  else {
    lines(df$end_date[1:181], y)
  }
}
axis.Date(1, at=c(as.Date("2017/01/14"), seq(as.Date('2017/03/01'), as.Date('2020/01/01'), by="2 mon"),
                  as.Date('2020/03/08'), as.Date('2020/03/28'), as.Date('2020/05/01'),
                  as.Date('2020/06/27')), format="%Y-%m-%d", las = 2, cex.axis = 0.5, tck = -0.03)
axis(2, at = seq(2.7, 46.92, by = 4.02), labels = plot_states, las = 2, cex.axis = 0.7, tck = 0, lwd = 0)
abline(v=df$end_date[first_test_week]-6, lty = 2)
abline(v=as.Date('2020-03-15'), col = "#00FFFF")
abline(v=as.Date('2020-03-28'), col = "#00FFFF")
