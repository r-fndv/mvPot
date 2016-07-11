#' Generating vectors for lattice rules
#'
#' Compute an efficient generating vector for quasi-Monte Carlo estimation.
#'
#' The function computes a generating vector for efficient multivariate integral estimation
#' based on D. Nuyens and R. Cools (2004).
#'
#' @param p Prime number corresponding to the number of samples to use in the quasi-Monte Carlo procedure.
#' @param d Dimension of the multivariate integral to estimate.
#' @param bt Tuning parameter for finding the vector. See D. Nuyens and R. Cools (2004) for more details.
#' @param gm Tuning parameter for finding the vector. See D. Nuyens and R. Cools (2004) for more details.
#' @return A \code{d}-dimensional generating vector for lattice rules with \code{p} samples.
#' @examples
#' #Define the number of sample. Muste be a prime number.
#' p <- 499
#'
#' #Choose a dimension
#' d <- 300
#'
#' #Compute the generating vector
#' genVecQMC(p,d)
#'
#' @export


genVecQMC <- function(p, d, bt = rep(1,d), gm = c(1, (4/5)^(0:(d-2)))){

  if(numbers::isPrime(p) == 0){
    p <- numbers::previousPrime(p)
  }

  q <- 1
  w <- 1
  z <- (1:d)
  m <- (p - 1) / 2

  #Find primitive root
  g <- prmrot(p)
  perm <- 1:m

  for(j in (1:(m-1))){
    perm[j+1] = (g*perm[j]) %% p
  }

  perm <- pmin(p - perm, perm)
  c <- (perm /p)^2 - perm / p + 1/6
  fc <- stats::fft(c)
  bump <- rep(1,d)

  for(s in 2:d){
    q = q * (bt[s - 1] + gm[s - 1]*c[c(w:1, m:(w+1))])
    w <- which.min(Re(stats::fft( fc * stats::fft(q), inverse = TRUE) / length(fc)))
    bump[s] <- w
    z[s] <- perm[w]
  }
  return(z)
}


prmrot <- function(p){

  pm <- p - 1
  fp <- unique(gmp::factorize(pm))
  n <- length(fp)
  r <- 2
  k <- 1

  while(k <= n){
    denom <- as.integer(pm / fp[k])
    rd <- r
    for(i in 2:denom){
      rd <- (rd * r) %% p
    }# computes r^denom mod p
    if(rd == 1) {
      r <- r + 1
      k <- 0;
    }
    k <- k + 1;
  }

  return(r)
}
