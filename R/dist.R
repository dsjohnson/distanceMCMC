#' Distance sampling distributions that can be used in \code{nimble} models
#'
#' \code{ddist} distance sampling distributions for half-normal and hazard rate models.
#'
#' @aliases ddist rdist
#'
#' @name ddist
#' @param x distance observations. either a single value (dHN) or a vector of values (dHN_V)
#' @param sigma scale parameter for detection models
#' @param b Power parameter for hazard rate model.
#' @param w right truncation distance for integration of the likelihood function
#' @param log if TRUE, return the log-likelihood
#' @param n   number of random values to generate
#' @param prob Logical. Should the function return the marginal detection
#' probability instead of effective strip width (ESW)
#' @param pars Parameters for hazard rate function. Used in `hr_integrand()`. Not intended
#' for practitioner use.
#' @author Devin Johnson
#' @rdname ddist
#' @import nimble
#' @importFrom stats pnorm
#' @export
ddist_hn <- nimbleFunction(
  run = function(x = double(0), sigma = double(0), w = double(0), log = integer(0, default = 0)) {
    returnType(double(0))
    if (x < 0 | x > w | sigma <= 1e-10) {
      if (log) return(-Inf) else return(0)
    }
    log_g <- -x^2 / (2 * sigma^2)
    
    # Avoid catastrophic cancellation in pnorm(w, 0, sigma) - 0.5 when sigma is very large
    if (w / sigma < 1e-5) {
      # When sigma >> w, half-normal approaches Uniform(0, w), so esw ≈ w
      log_esw <- log(w)
    } else {
      log_esw <- 0.5 * log(2 * pi) + log(sigma) + log(pnorm(w, 0, sigma) - 0.5)
    }
    
    logL <- log_g - log_esw
    if (log) return(logL) else return(exp(logL))
  }
)


#' @rdname ddist
#' @export
ddist_hr <- nimbleFunction(
	run = function(x = double(0),
								 sigma = double(0),
								 b = double(0),
								 w = double(0),
								 log = double(0, default = 0)) {
		if(x > w){
			L <- 0
		} else{
			pars = c(sigma, b)
			g <- 1 - exp(-(x/sigma)^(-b))
			esw <- nimIntegrate(hr_integrand, lower = 0, upper = w, param = pars)[1]
			L <- g/esw
		}
		if(log) return(log(L))
		else return(L)
		returnType(double(0))
	}
)


### ----------------------------------------------------------------------------
### Distance detection functions that need numerical integration
### ----------------------------------------------------------------------------

#' @rdname ddist
#' @export
hr_integrand <- nimbleFunction(
	run = function(x = double(1), pars = double(1)){
		res <- 1 - exp(-(x/pars[1])^(-pars[2]))
		return(res)
		returnType(double(1))
	}
)

### ----------------------------------------------------------------------------


#' @rdname ddist
#' @importFrom stats pnorm rbinom rnorm
#' @export
rdist_hn <- nimbleFunction(
	run = function(n = integer(),
								 sigma = double(0),
								 w = double(0)
	) {
		x <- w + 1
		it <- 1
		while(x <= w & it<=100){
			x <- rnorm(1,0,sigma)
			x <- abs(x)
			it <- it+1
		}
		returnType(double(0))
		return(x)
	}
)

#' @rdname ddist
#' @importFrom stats pnorm rbinom runif
#' @export
rdist_hr <- nimbleFunction(
	run = function(n = integer(0),
								 sigma = double(0),
								 b = double(0),
								 w = double(0)
	) {
		k <- 0
		it <- 1
		while(k == 0 & it<=100){
			x <- runif(1, 0, w)
			g <- 1 - exp(-(x/sigma)^(-b))
			k <- rbinom(1,1,g)
			it <- it+1
		}
		returnType(double(0))
		return(x)
	}
)


#' @rdname ddist
#' @importFrom stats pnorm
#' @export
esw_dist_hn <- nimbleFunction(
	run = function(sigma = double(1),
								 w = double(0),
								 prob=double(0, default=0)) {
		esw <- 2*(pnorm(w, 0, sigma)-0.5) * sqrt(pi/2) * sigma
		returnType(double(1))
		if(prob==1) return(esw/w)
		else return(esw)
	}
)

#' @rdname ddist
#' @importFrom stats pnorm
#' @export
esw_dist_hr <- nimbleFunction(
	run = function(sigma = double(1),
								 b = double(1),
								 w = double(0),
								 prob=double(0, default=0)) {
		n <- length(sigma)
		esw <- numeric(n)
		for(i in 1:n){
			pars <- c(sigma[i], b[i])
			esw[i] <- nimIntegrate(hr_integrand, lower = 0, upper = w, param = pars)[1]
		}
		returnType(double(1))
		if(prob==1) return(esw/w)
		else return(esw)
	}
)


