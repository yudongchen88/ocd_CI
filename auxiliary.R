# This R file includes all auxiliary functions

# Norm of a vector
vector.norm <- function(v, q = 2, na.rm = FALSE){
  if (na.rm) v <- na.omit(v)
  if (q > 0) (sum(abs(v)^q))^(1/q)
  else NaN
}


# Normalise a vector
vector.normalise <- function(v, q = 2, na.rm = FALSE){
  return(v / vector.norm(v, q, na.rm))
}


# generate a vector, with given dimension, sparsity and signal_shape
vector_of_change <- function(p, s=p, vartheta=1, signal_shape='random'){
  signal_template <- list(rep(1,s), rnorm(s), seq_len(s), seq_len(s)^(-1), seq_len(s)^(-0.5),
                          seq_len(s)^(-2))
  names(signal_template) <- c('uniform', 'random', 'linear', 'harmonic', 'inv_sqrt','inv_quad')
  theta <- c(sample(c(-1, 1), s, replace=T)*signal_template[[signal_shape]], rep(0, p-s))
  theta <- vector.normalise(theta) * vartheta
  return(theta)
}


# cusum transformation
cusum.transform <- function(x){
  x <- as.matrix(x)
  if (dim(x)[2] == 1) x <- t(x) # treat univariate time series as a row vector
  p <- dim(x)[1] # dimensionality of the time series
  n <- dim(x)[2] # time length of the observation
  
  leftsums <- t(apply(x, 1, cumsum))
  rightsums <- leftsums[, n] - leftsums
  
  t <- 1:(n - 1)
  
  # constructing CUSUM matrix
  rightmeans <- sweep(rightsums[, t, drop=FALSE], 2, n - t, '/')
  leftmeans <- sweep(leftsums[, t, drop=FALSE], 2, t, '/')
  cusum <- sweep(rightmeans - leftmeans, 2, sqrt(t * (n-t) / n), '*')
  return(cusum)
}


# sparse SVD
sparse.svd <- function(Z, lambda){
  # with Frobenius norm constraint, the sparse vector is obtained by soft
  # thresholding
  Mhat <- sign(Z) * pmax(0, abs(Z) - lambda)
  
  # compute the leading left singular vector of Mhat
  if (sum(Mhat^2)!=0){
    # compute the leading left singular vector
    vector.proj <- svd(Mhat)$u[,1]
    
  } else {
    # if the thresholded matrix is zero, return a random vector
    vector.proj <- rnorm(ncol(Z))
    vector.proj <- vector.proj / vector.norm(vector.proj)
  }
  return(vector.proj)
}


# InspectChangepoint (Wang and Samworth, 2018)
locate.change <- function(x, lambda, schatten=2)
{
  x <- as.matrix(x)
  if (dim(x)[2] == 1) x <- t(x) # treat univariate time series as a row vector
  p <- dim(x)[1] # dimensionality of the time series
  n <- dim(x)[2] # time length of the observation
  if (missing(lambda)) lambda <- sqrt(log(log(n)*p)/2)
  x1 <- x
  x2 <- x
  
  # construct cusum matrix of x
  cusum.matrix1 <- cusum.transform(x1)
  cusum.matrix2 <- cusum.matrix1
  
  # estimate changepoint
  if (lambda >= max(abs(cusum.matrix1))) lambda <- max(abs(cusum.matrix1)) - 1e-10
  
  vector.proj <- sparse.svd(cusum.matrix1, lambda);
  cusum.proj <- t(cusum.matrix2)%*%vector.proj
  
  
  ret <- NULL
  ret$changepoint <- which.max(abs(cusum.proj))
  ret$cusum <- max(abs(cusum.proj))
  ret$vector.proj <- as.numeric(vector.proj)
  
  return(ret)
}


# Choosing thresholding parsmeter using BIC (Kaul et al., 2021)
Q_tune <- function(tau, X, mean_pre_thresh, mean_post_thresh){
  X_pre_tau <- X[ , 1:tau]
  X_post_tau <- X[ , (tau+1):ncol(X)]
  Q_value <- (vector.norm(X_pre_tau - mean_pre_thresh, 2))^2 + (vector.norm(X_post_tau - mean_post_thresh, 2))^2
  return(Q_value)
}


BIC_tune <- function(lambda, X_pre, X_post){
  
  mean_pre <- apply(as.matrix(X_pre), 1, mean)
  mean_pre_thresh <- sign(mean_pre) * pmax(abs(mean_pre) - lambda, 0)
  X_pre_thresh <- X_pre - mean_pre_thresh
  v_pre <- (vector.norm(X_pre_thresh, 2))^2
  
  mean_post <- apply(as.matrix(X_post), 1, mean)
  mean_post_thresh <- sign(mean_post) * pmax(abs(mean_post) - lambda, 0)
  X_post_thresh <- X_post - mean_post_thresh
  v_post <- (vector.norm(X_post_thresh, 2))^2
  
  count <- sum((mean_pre_thresh != 0) | (mean_post_thresh != 0))
  
  BIC_value <- v_pre + v_post + count * log(ncol(X_pre) + ncol(X_post))
  return(BIC_value)
}


# Find first non-zero element in vector
find.first <- function(v){
  match(TRUE, v, nomatch = NA)
}


# Print line
println <- function(...){
  .Internal(cat(c(list(...), '\n'), file=stdout(), sep='', fill=FALSE, labels=NULL, append=FALSE))
}


# Multiple assigment
'%=%' <- function(l, r) UseMethod('%=%')


# Binary operator
'%=%.lbunch' = function(l, r) {
  Envir = as.environment(-1)
  if (length(r) > length(l))
    warning("RHS has more args than LHS. Only first", length(l), "used.")
  if (length(l) > length(r))  {
    warning("LHS has more args than RHS. RHS will be repeated.")
    r <- extendToMatch(r, l)
  }
  for (i in seq_along(l)) {
    do.call('<-', list(l[[i]], r[[i]]), envir=Envir)
  }
}


# Grouping the left hand side in multiple assignment
bunch <- function(...) {
  List <- as.list(substitute(list(...)))[-1L]
  class(List) <- 'lbunch'
  return(List)
}


# Matplotlib palette
matplotlib_palette <- function(n=0, scheme='default', visualise=FALSE){
  default_palette <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd",
                       "#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf")
  bright_palette <- c("#0e4897","#17813f","#1d99b4","#1f9ee8","#25ca7a","#471c7c",
                      "#68c7ed","#6d4e98","#73af38","#7f1273","#9e1653","#ab0077",
                      "#b01426","#b1b2b4","#c1d430","#cc0b24","#e10064","#e12653",
                      "#e34e9d","#e46b07","#fbee29","#fcc125")
  rainbow_palette <- c("#BF4D4D","#BF864D","#BFBF4D","#86BF4D","#4DBF4D",
                       "#4DBF86","#4DBFBF","#4D86BF","#4D4DBF","#864DBF",
                       "#BF4DBF","#BF4D86")
  full_palette <- switch(scheme,
                         'default' = default_palette,
                         'bright' = bright_palette,
                         'rainbow' = rainbow_palette)
  
  if (n == 0) {
    ret <- full_palette
  } else {
    reps <- ceiling(n / length(full_palette))
    ret = character()
    for (rep in 1:reps){
      mod_color <- unname(sapply(full_palette, function(c)mix_color(c, '#ffffff', (rep-1)/reps)))
      if (rep==reps) mod_color <- head(mod_color, n - length(full_palette)*(reps-1))
      ret <- c(ret, mod_color)
    }
  }
  
  if (visualise) barplot(rep(1, length(ret)), col=ret, axes=F, border=F,
                         names.arg=seq_along(ret), cex.names=0.8)
  return(ret)
}


# Color mixing
mix_color <- function(color1, color2, lambda=0.5){
  R1 <- strtoi(paste0('0x', substr(color1, 2, 3)))
  G1 <- strtoi(paste0('0x', substr(color1, 4, 5)))
  B1 <- strtoi(paste0('0x', substr(color1, 6, 7)))
  R2 <- strtoi(paste0('0x', substr(color2, 2, 3)))
  G2 <- strtoi(paste0('0x', substr(color2, 4, 5)))
  B2 <- strtoi(paste0('0x', substr(color2, 6, 7)))
  R <- round(R1 * (1-lambda) + R2 * lambda, 0)
  G <- round(G1 * (1-lambda) + G2 * lambda, 0)
  B <- round(B1 * (1-lambda) + B2 * lambda, 0)
  return(rgb(R, G, B, maxColorValue = 255))
}


