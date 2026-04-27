# # This file registers all distributions when the package is loaded.
.onAttach <- function(libname, pkgname) {

	registerDistributions(list(
		ddist_hn = list(
			BUGSdist = "ddist_hn(sigma, w)",
			Rdist = "ddist_hn(sigma, w)",
			discrete = FALSE,
			types = c('value = double(0)',
								'sigma = double(0)',
								'w = double(0)'),
			mixedSizes = TRUE,
			pqAvail = FALSE)), verbose = F
	)

	registerDistributions(list(
		ddist_hr = list(
			BUGSdist = "ddist_hr(sigma, b, w)",
			Rdist = "ddist_hr(sigma, b, w)",
			discrete = FALSE,
			types = c('value = double(0)',
								'sigma = double(0)',
								'b = double(0)',
								'w = double(0)'),
			mixedSizes = TRUE,
			pqAvail = FALSE)), verbose = F
	)

}
