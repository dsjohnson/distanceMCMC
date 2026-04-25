# # This file registers all distributions when the package is loaded.
.onAttach <- function(libname, pkgname) {

	registerDistributions(list(
		ddist = list(
			BUGSdist = "ddist(sigma, b, w, form)",
			Rdist = "ddist(sigma, b, w, form)",
			discrete = FALSE,
			types = c('value = double(0)',
								'sigma = double(0)',
								'b = double(0)',
								'w = double(0)',
								'form = integer(0)'),
			mixedSizes = TRUE,
			pqAvail = FALSE)), verbose = F
	)

}
