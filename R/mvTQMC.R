#' Multivariate t distribution function
#'
#' Estimate the multivariate t distribution function with quasi-Monte Carlo method.
#'
#' The function uses a quasi-Monte Carlo procedure based on randomly shifted
#' lattice rules to estimate the distribution function a multivariate normal distribution
#' as described in Genz, A. and Bretz, F.(2009).
#'
#' @author Raphael de Fondeville
#' @param p Number of samples used for quasi-Monte Carlo estimation. Must be a prime number.
#' @param upperBound Vector of probabilities, i.e., the upper bound of the integral.
#' @param cov Covariance matrix of the multivariate normal distribution. Must be semi-positive definite.
#' WARNING: for performance in high-dimensions, no check is performed on the matrix. It is the user responsibility to ensure
#' that this property is verified.
#' @param nu Degrees of freedom of the t distribution.
#' @param genVec Generating vector for the quasi-Monte Carlo procedure. Can be computed using \code{genVecQMC}.
#' @return An estimate of the distribution function along with empirical Monte Carlo error.
#' @examples
#'
#' #Define locations
#' loc <- expand.grid(1:4, 1:4)
#' ref <- sample.int(16, 1)
#'
#' #Defined degrees of freedom
#' nu <- 3
#'
#' #Compute variogram matrix
#' variogramMatrix <- ((sqrt((outer(loc[,1],loc[,1],"-"))^2 +
#' (outer(loc[,2],loc[,2],"-"))^2)) / 2)^(1.5)
#'
#' #Define an upper boud
#' upperBound <- variogramMatrix[-ref,ref]
#'
#' #Compute covariance matrix
#' cov <-  (variogramMatrix[-ref,ref]%*%t(matrix(1, (nrow(loc) - 1), 1)) +
#' t(variogramMatrix[ref,-ref]%*%t(matrix(1, (nrow(loc) - 1), 1))) -
#' variogramMatrix[-ref,-ref])
#'
#' #Compute generating vector
#' p <- 499
#' latticeRule <- genVecQMC(p, (nrow(loc) - 1))
#'
#' #Estimate the multivariate distribution function
#' mvTProbQuasiMonteCarlo(latticeRule$primeP, upperBound, cov, nu, latticeRule$genVec)
#' @export
#' @useDynLib mvPot mvTProbCpp
#' @references Genz, A. and Bretz, F. (2009). Computations of Multivariate Normal and t Probabilities, volume 105. Springer, Dordrecht.
#'
#'             Genz, A. (2013). QSILATMVTV \url{http://www.math.wsu.edu/faculty/genz/software/software.html}

mvTProbQuasiMonteCarlo = function(p, upperBound, cov, nu, genVec){

  if(!is.numeric(p) | length(p) > 1) {
    stop('p must be a prime number.')
  }

  if(!is.numeric(upperBound) | !is.vector(upperBound)) {
    stop('upperBound must be a numeric vector.')
  }

  if(!is.numeric(cov) | !is.matrix(cov)) {
    stop('upperBound must be a numeric matrix')
  }

  if(nrow(cov) != length(upperBound) | ncol(cov) != length(upperBound)){
    stop('upperBound and cov must have the same dimension')
  }

  if(!is.numeric(genVec) | !is.vector(genVec) | length(upperBound) != length(genVec)){
    stop('genVec must be a numeric vector of the same size as upperBound')
  }

  if(!is.numeric(nu) | nu <= 0 | length(nu) > 1 ){
    stop('nu must be a positive number')
  }

  tmp <-.C(mvTProbCpp,
           as.integer(p),
           as.integer(length(upperBound)),
           as.double(cov),
           as.double(upperBound),
           as.double(nu),
           as.double(genVec),
           est = double(length=1),
           err = double(length=1)
  )

  c(estimate = tmp$est, error = tmp$err)
}
