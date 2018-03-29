#' Simulation of extremal Student generalized Pareto vectors
#' 
#' The algorithm is described in section 4 of Thibaud and Opitz. It uses the Cholesky decomposition of the matrix \code{Sigma}
#' to generate samples on the unit sphere and an accept-reject algorithm to decide which samples to retain. If \code{normalize = TRUE} (the default), 
#' the vector is scaled by the exponent measure \eqn{\kappa} so that the maximum of the sample is greater than \eqn{\kappa}.
#' 
#' @param n sample size
#' @param Sigma correlation matrix
#' @param nu degrees of freedom parameter
#' @param normalize logical; should unit Pareto samples be returned?
#' @param matchol Cholesky matrix \eqn{\mathbf{A}}{A} such that \eqn{\mathbf{A}\mathbf{A}^\top = \boldsymbol{\Sigma}}{AA^t = \Sigma}. Corresponds to \code{t(chol(Sigma))}. Default to \code{NULL}, in which case the Cholesky root is computed within the function.
#' @references Thibaud, E. and T. Opitz (2015). Efficient inference and simulation for elliptical Pareto processes. Biometrika, 102(4), 855-870.
#' @author Emeric Thibaud, Leo Belzile
#' @export
#' @importFrom stats cov2cor rnorm runif
#' @examples 
#' loc <- expand.grid(1:4, 1:4)
#' Sigma <- exp(-as.matrix(dist(loc))^1.5)
#' rExtremalStudentParetoProcess(10000, Sigma, nu = 2)
rExtremalStudentParetoProcess <- function(n, Sigma, nu, normalize = TRUE, matchol = NULL){
  d <- nrow(Sigma)
  # Check the input is indeed a correlation matrix
  if(!isTRUE(all.equal(as.vector(diag(Sigma)), rep(1, d)))){
    if(all(c(isSymmetric(Sigma), eigen(Sigma, only.values = TRUE)$eigen > 1e-10))){
      warning("`Sigma` is not a corrrelation matrix.")
      Sigma <- cov2cor(Sigma)
    } else{
      stop("The correlation function provided by the user is not valid.")
    }
  }
  if (!is.null(matchol)){
    A <- matchol
  } else {
    A <- t(chol(Sigma))
  }
  
  # Container for simulations
  P <- matrix(0, nrow = n, ncol = d)
  #Compute the exponent measure
  if (normalize){
    genVec <- mvPot::genVecQMC(p = 499, d - 1)
    kap <- sum(sapply(1:d, function(j){mvPot::mvTProbQuasiMonteCarlo(p = genVec$primeP,
                                                  upperBound = 1 - Sigma[-j,j],
                                                  cov = (Sigma[-j, -j] - Sigma[-j, j, drop = FALSE] %*% Sigma[j, -j, drop = FALSE]) / (nu + 1),
                                                  nu = nu + 1, genVec = genVec$genVec)[1]}))
  } else {
    kap <- 1
  }
  #Simulation algorithm (Theorem 5.3 in Thomas Opitz's thesis, p. 80)
  N <- n # Number of simulations
  ni <- 0L # Number of remaining simulations
  while(ni < n){
    R <- 1 / runif(N) # Unit Pareto variates
    V <- matrix(rnorm(N * d), ncol = d) # Standard Normal variates
    normV <- sqrt(rowSums(V^2)) # Compute L2 norm
    U <- V / normV # Create uniform on sphere
    tAU <- t(A %*% t(U)) 
    tAU[tAU < 0] <- 0 # Truncated negative components
    W <- R * (tAU^nu)
    In <- apply(W, 1, max) > 1 # Check accept-reject condition
    acc <- sum(In) # Number of simulations accepted
      if(acc > 0){
        P[(ni + 1):min(n, ni + acc), ] <- W[which(In)[1:(min(n - ni, acc))], ]
        ni <- ni + acc;
        N <- max(5, ceiling((n - ni) * N / acc)) #Number of simulations for next rounds
      }
  }
    return(kap * P)
  
}

