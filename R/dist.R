#' Distance sampling distributions that can be used in \code{nimble} models
#'
#' \code{ddist} distance sampling distributions for half-normal and hazard rate models.
#'
#' @aliases ddist rdist
#'
#' @name ddist
#' @param x distance observations. either a single value (dHN) or a vector of values (dHN_V)
#' @param sigma scale parameter for detection models
#' @param b Power parameter for hazard rate model. Ignored for Half-normal (`form==1`)
#' @param w right truncation distance for integration of the likelihood function
#' @param form Integer. Half-normal = 1, Hazard rate = 2
#' @param log if TRUE, return the log-likelihood
#' @param n   number of random values to generate
#' @param prob Logical. Should the function return the marginal detection
#' probability instead of effective strip width (ESW)
#' @author Devin Johnson
#' @rdname ddist
#' @import nimble
#' @importFrom stats pnorm
#' @export
ddist <- nimbleFunction(
	run = function(x = double(0),
								 sigma = double(0),
								 b = double(0, default=1),
								 w = double(0),
								 form = integer(0, default=1),
								 log = double(0, default = 0)) {
		if(form==1){
			g <- exp(-x^2/(2*sigma^2))
			esw <- 2*(pnorm(w, 0, sigma)-0.5) * sqrt(pi/2) * sigma
		} else{
			g <- 1 - exp(-(x/sigma)^(-b))
			pars <- c(sigma, b)
			esw <- nimIntegrate(hr_integrand, lower = 0, upper = w, param = pars)[1]
		}
		L <- g/esw
		if(log) return(log(L))
		else return(L)
		returnType(double(0))
	}
)

### ----------------------------------------------------------------------------
### Distance detection functions that need numerical integration
### ----------------------------------------------------------------------------

hr_integrand <- nimbleFunction(
	run = function(x = double(1), pars = double(1)){
		res <- 1 - exp(-(x/pars[1])^(-pars[2]))
		return(res)
		returnType(double(1))
	}
)

### ----------------------------------------------------------------------------


#' @rdname ddist
#' @importFrom stats pnorm rbinom runif
#' @export
rdist<- nimbleFunction(
	run = function(n = integer(),
								 sigma = double(0),
								 b = double(0, default=1),
								 w = double(0),
								 form = integer(0, default=1)
	) {
		k <- 0
		it <- 1
		while(k == 0 & it<=100){
			x <- runif(1, 0, w)
			if(form==1) g <- exp(-x^2/(2*sigma^2))
			else g <- 1 - exp(-(x/sigma)^(-b))
			d <- rbinom(1,1,g)
			it <- it+1
			if(d==1){
				out <- x
			}
		}
		returnType(double(0))
		return(out)
	}
)

#' @rdname ddist
#' @importFrom stats pnorm
#' @export
esw_dist <- nimbleFunction(
	run = function(sigma = double(0),
								 b = double(0, default=1),
								 w = double(0),
								 form = integer(0, default=1),
								 prob=double(0, default=0)) {
		if(form==1){
			esw <- 2*(pnorm(w, 0, sigma)-0.5) * sqrt(pi/2) * sigma
		} else{
			pars <- c(sigma, b)
			esw <- nimIntegrate(hr_integrand, lower = 0, upper = w, param = pars)[1]
		}
		returnType(double(0))
		if(prob==1) return(esw/w)
		else return(esw)
	}
)

