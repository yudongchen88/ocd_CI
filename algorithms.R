# initialise the partial sum matrix 'A' and tail length matrix 'tail'
initialise <- function(p){
  A <- matrix(0, p, 1) # single column, corresponding to tail = 0
  colnames(A) <- 'tail=0'
  tail <- matrix(0, p, floor(log2(p))*2+4) # all tail lengths initialised to 0
  return(list(A=A, tail=tail))
}

# compute all relevant scales for a given \ell_2 norm lower bound of signal
get_scales <- function(beta, p){
  L <- floor(log2(p))  # number of different scales
  tmp <- beta/sqrt(2^(0:(L+1)) * log2(2*p))
  B <- c(tmp, -tmp) # diadic grid of univariate alternative values
  return(B)
}

# process a new observation using the ocd algorithm
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
  
  # compute test stats S_diag, S_dense and S_sparse
  G <- sweep(A^2, 2, pmax(1, unique_tail_new), '/')  # normalisation
  G0 <- G * (G > a^2) # hard thresholded
  colsum_sparse <- colSums(G0)
  
  for (i in seq_along(unique_tail_new)) {
    colsum_sparse[i] <- colsum_sparse[i] - min(G0[apply(tail == unique_tail_new[i], 1, any), i])
  }
  
  S_diag <- max(R, 0)
  S_sparse <- max(colsum_sparse)
  stat <- setNames(c(S_diag, S_sparse), c('diag', 'sparse'))
  
  return(list(stat=stat, A=A, tail=tail))
}

# compute the anchor coordinate and anchor scale
find_anchor <- function(A, tail, beta){
  p <- nrow(tail)   # dimension of data
  a <- sqrt(2*log(p))  # hard threshold parameter
  B <- get_scales(beta, p)
  
  unique_tail <- sort(unique(as.vector(tail)))
  tail_loc <- matrix(match(as.vector(tail), unique_tail), p)
  
  # expand A into a p x |B| matrix, and compute a p x |B| matrix R of tail loglik ratio stats
  A_expand <- matrix(A[cbind(rep(1:p, length(B)), as.vector(tail_loc))], p)
  R <- t(t(A_expand) * B - t(tail) * B^2 / 2)
  
  # compute test stats S_diag, S_dense and S_sparse
  G <- sweep(A^2, 2, pmax(1, unique_tail), '/')  # normalisation
  G0 <- G * (G > a^2) # hard thresholded
  colsum_sparse <- colSums(G0)
  
  anchor_co_sparse <- rep(0, length(colsum_sparse))
  
  for (i in seq_along(unique_tail)) {
    filter <- which(apply(tail == unique_tail[i], 1, any)) # coordinates using this tail length
    colsum_sparse[i] <- colsum_sparse[i] - min(G0[filter, i])
    anchor_co_sparse[i] <- filter[which.min(G0[filter, i])] # anchor coord for ith unique tail
  }
  
  tailind <- unname(which.max(colsum_sparse)) # anchor tail index in the unique tails
  anchor_tail_length <- unique_tail[tailind] # anchor tail length
  anchor_co <- anchor_co_sparse[tailind] # anchor coordinate
  
    
  tail_trigger = unique_tail[tailind]
  return(list(anchor_tail_length=anchor_tail_length, anchor_coord=anchor_co))
}

# find the largest (in magnitude) scale below (abs(E) - d1) / sqrt(t) in B
biggest_scale <- function(E, B, d1, t){
  filter1 <- abs(B) <= (abs(E) - d1)/sqrt(t) # condition on \tilde{b}_j
  filter2 <- sign(B) == sign(E)
  ind <- match(TRUE, filter1 & filter2, nomatch = NA)
  return(B[ind])
}

# construct support estimate of the vector of mean change at time of declaration
support_estimate <- function(A, tail, anchor_coord, anchor_tail_length, beta, d1){
  p <- nrow(tail)   # dimension of data
  B <- get_scales(beta, p)
  
  E_jhat <- A[, paste0('tail=', anchor_tail_length)] / sqrt(anchor_tail_length)
  b_tilde <- sapply(E_jhat, biggest_scale, B, d1, anchor_tail_length)
  b_tilde[anchor_coord] <- NA
  S_hat <- which(!is.na(b_tilde))
  
  return(list(S_hat=S_hat, b_tilde=b_tilde))
}

# construct a confidence interval at the time of declaration for changepoint location
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
