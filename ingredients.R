# This R file includes main ingredient functions of the ocd_CI procedure,
# as well as an implementation of a procedure based on Kaul et al. (2021).

# Initialise the partial sum matrix 'A' and tail length matrix 'tail'
initialise <- function(p){
  A <- matrix(0, p, 1) # single column, corresponding to tail = 0
  colnames(A) <- 'tail=0'
  tail <- matrix(0, p, floor(log2(p))*2+4) # all tail lengths initialised to 0
  return(list(A=A, tail=tail))
}

# Compute all relevant scales for a given \ell_2 norm lower bound of signal
get_scales <- function(beta, p){
  L <- floor(log2(p))  # number of different scales
  tmp <- beta/sqrt(2^(0:(L+1)) * log2(2*p))
  B <- c(tmp, -tmp) # diadic grid of univariate alternative values
  return(B)
}

# Process a new observation using the ocd algorithm
process_new_obs <- function(x_new, A, tail, beta, unique_tail){
  p <- length(x_new)   # dimension of data
  a <- sqrt(2*log(p))  # hard threshold parameter
  B <- get_scales(beta, p)
  
  # find unique elements in the tail matrix, and establish a match
  unique_tail <- sort(unique(as.vector(tail)))
  tail_loc <- matrix(match(as.vector(tail), unique_tail), p)
  
  # update A and tail first
  A <- A + x_new
  tail <- tail + 1
  
  # expand A into a p x |B| matrix, and compute a p x |B| matrix R of tail loglik ratio stats
  A_expand <- matrix(A[cbind(rep(1:p, length(B)), as.vector(tail_loc))], p)
  R <- t(t(A_expand) * B - t(tail) * B^2 / 2)
  
  # reset some tail and CUSUM in A to zero
  tail <- tail * (R > 0)
  unique_tail_new <- sort(unique(as.vector(tail)))
  if (unique_tail_new[1] == 0){
    A <- A[, match(unique_tail_new[-1] - 1, unique_tail)]
    A <- cbind(0, A)
  } 
  colnames(A) <- paste0('tail=', unique_tail_new)
  
  # compute test stats S_diag, S_off
  G <- sweep(A^2, 2, pmax(1, unique_tail_new), '/')  # normalisation
  G0 <- G * (G > a^2) # hard thresholded
  colsum_off <- colSums(G0)
  
  for (i in seq_along(unique_tail_new)) {
    colsum_off[i] <- colsum_off[i] - min(G0[apply(tail == unique_tail_new[i], 1, any), i])
  }
  
  S_diag <- max(R, 0)
  S_off <- max(colsum_off)
  stat <- setNames(c(S_diag, S_off), c('diag', 'off'))
  
  return(list(stat=stat, A=A, tail=tail))
}

# Compute the anchor coordinate and anchor scale
find_anchor <- function(A, tail, beta){
  p <- nrow(tail)   # dimension of data
  a <- sqrt(2*log(p))  # hard threshold parameter
  B <- get_scales(beta, p)
  
  unique_tail <- sort(unique(as.vector(tail)))
  tail_loc <- matrix(match(as.vector(tail), unique_tail), p)
  
  # expand A into a p x |B| matrix, and compute a p x |B| matrix R of tail loglik ratio stats
  A_expand <- matrix(A[cbind(rep(1:p, length(B)), as.vector(tail_loc))], p)
  R <- t(t(A_expand) * B - t(tail) * B^2 / 2)
  
  # compute test stats S_diag, S_dense and S_off
  G <- sweep(A^2, 2, pmax(1, unique_tail), '/')  # normalisation
  G0 <- G * (G > a^2) # hard thresholded
  colsum_off <- colSums(G0)
  
  anchor_co_off <- rep(0, length(colsum_off))
  
  for (i in seq_along(unique_tail)) {
    filter <- which(apply(tail == unique_tail[i], 1, any)) # coordinates using this tail length
    colsum_off[i] <- colsum_off[i] - min(G0[filter, i])
    anchor_co_off[i] <- filter[which.min(G0[filter, i])] # anchor coord for ith unique tail
  }
  
  tailind <- unname(which.max(colsum_off)) # anchor tail index in the unique tails
  anchor_tail_length <- unique_tail[tailind] # anchor tail length
  anchor_co <- anchor_co_off[tailind] # anchor coordinate
  
    
  tail_trigger = unique_tail[tailind]
  return(list(anchor_tail_length=anchor_tail_length, anchor_coord=anchor_co))
}

# Find the largest (in magnitude) scale below (abs(E) - d1) / sqrt(t) in B
biggest_scale <- function(E, B, d1, t){
  filter1 <- abs(B) <= (abs(E) - d1)/sqrt(t) # condition on \tilde{b}_j
  filter2 <- sign(B) == sign(E)
  ind <- match(TRUE, filter1 & filter2, nomatch = NA)
  return(B[ind])
}

# Construct support estimate of the vector of mean change at time of declaration
support_estimate <- function(A, tail, anchor_coord, anchor_tail_length, beta, d1){
  p <- nrow(tail)   # dimension of data
  B <- get_scales(beta, p)
  
  E_jhat <- A[, paste0('tail=', anchor_tail_length)] / sqrt(anchor_tail_length)
  b_tilde <- sapply(E_jhat, biggest_scale, B, d1, anchor_tail_length)
  b_tilde[anchor_coord] <- NA
  S_hat <- which(!is.na(b_tilde))
  
  return(list(S_hat=S_hat, b_tilde=b_tilde))
}

# Construct a confidence interval at the time of declaration for changepoint location
conf_int <- function(N, tail, S_hat, b_tilde, beta, d2){
  p <- nrow(tail)   # dimension of data
  B <- get_scales(beta, p)
  left_end <- 0
  for (j in S_hat){
    left_end <- max(N - tail[j, match(b_tilde[j], B)] - d2 / b_tilde[j]^2, 
                    left_end)
  }
  return(c(ceiling(left_end), N))
}

# Different types of supports
support <- function(p, s, vartheta, beta, signal_shape){
  b_min = beta/sqrt(2^(floor(log2(2*p)))*log2(2*p))
  theta = vector_of_change(p, s, vartheta, signal_shape)
  s_beta = sum(abs(theta) >= b_min)
  s_dyadic = 2^(0: floor(log2(p)))
  S_grid = vartheta/sqrt(s_dyadic*log2(2*p))
  s_eff = s_dyadic[find.first(colSums(outer(abs(theta), S_grid, '>=')) >= s_dyadic)]
  return(c(s_beta, s_eff))
}


#--------------------------------------------------------------

# An implementation of a procedure based on Kaul et al. (2021).
# Kaul CI
kaul <- function(X){
  ins <- locate.change(X)
  cp <- ins$changepoint
  N <- ncol(X)
  if (cp <= 0) cp <- 1
  if (cp >= N) cp <- N-1
  
  # divide data based on the estimated changepoint location
  X_pre <- X[ ,1:cp]
  X_post <- X[ , (cp+1):N]
  
  # tune lambda based on BIC
  Lambda_all <- seq(0.02, 0.5, 0.02)
  BIC_ret <- sapply(Lambda_all, BIC_tune, X_pre, X_post)
  lambda <- Lambda_all[which.min(BIC_ret)]
  
  mean_pre <- apply(as.matrix(X_pre), 1, mean)
  mean_post <- apply(as.matrix(X_post), 1, mean)
  mean_pre_thresh <- sign(mean_pre) * pmax(abs(mean_pre) - lambda, 0)
  mean_post_thresh <- sign(mean_post) * pmax(abs(mean_post) - lambda, 0)
  
  # refine change point estimation
  Q_ret <- sapply(1:(N-1), Q_tune, X, mean_pre_thresh, mean_post_thresh)
  cp_refine <- which.min(Q_ret)
  
  X_pre_refine <- X[ ,1:cp_refine]
  X_post_refine <- X[ , (cp_refine+1):N]
  
  mean_pre_refine <- apply(as.matrix(X_pre_refine), 1, mean) * (mean_pre_thresh != 0)
  mean_post_refine <- apply(as.matrix(X_post_refine), 1, mean) * (mean_post_thresh != 0)
  
  # estimate vartheta
  vartheta_est <- vector.norm(mean_pre_refine - mean_post_refine, 2)
  
  kaul_CI_left <- max(floor(cp_refine - 11.03 / vartheta_est^2), 0)
  kaul_CI_right <- min(ceiling(cp_refine + 11.03 / vartheta_est^2), N)
  return(c(kaul_CI_left, kaul_CI_right))
}

