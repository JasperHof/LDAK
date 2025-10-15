#' Compute transformed phenotype based on LDAK output.
#'
#' The transformation are used in combination with the original
#' covariate to compute a smoothed phenotype for association analysis.
#'
#' @param pheno Phenotype matrix or dataframe in PLINK format
#' @param cov Covariate vector
#' @param size Number of subgroups for stratifying the phenotype
#' @param trans Transformation matrix obtained from LDAK
#' @param spar Smoothing parameter (default: 0.2)
#' @return A transformed phenotypes which can be used for GWAS analysis
#' @importFrom stats lm predict quantile smooth.spline
#' @export
make_transform <- function(pheno = NULL, cov = NULL, size = numeric(), trans = NULL, spar = 0.2) {

  # Compute subset - and quants used in creating phenotype
  cov_cases = cov[which(pheno[,3] == 1)]
  quants = quantile(unique(as.numeric(cov_cases)), (1:(size-1))/size, na.rm = T)

  # Derive the mean quantile per subgroup - this might be non-trivial for discrete data
  quants_split <- cut(as.numeric(cov_cases),
                      breaks = c(-Inf, quants, Inf),  # use -Inf and Inf to cover all values
                      labels = FALSE,                 # numeric codes
                      right = TRUE)
  medians = sapply(1:size, function(i) median(cov_cases[quants_split == i], na.rm = T))
  quant_medians = sapply(1:size, function(i) (2 * rank(cov_cases)[which.min(abs(cov_cases - medians[i]))] - 1) / (2 * max(rank(cov_cases))))

  # Compute new phenotype by smoothing transformed phenotype
  transformed <- matrix(NA, dim(pheno)[1], 5)

  for (i in 1:5) {

    # Fit a quadratic function
    y <- as.numeric(trans[i,])
    fit <- smooth.spline(quants_medians, y, spar = 0.2)

    # Map covariate values to quantiles, and compute predicted transformation
    quantiles <- (2 * rank(cov_cases) - 1) / (2 * max(rank(cov_cases)))
    y_pred <- predict(fit, quantiles)$y

    # Add these values to new transformed phenotype matrix
    transformed[which(pheno[,3] == 1),i] = y_pred
  }

  # Fill in missing values and scale
  transformed[is.na(transformed)] = 0
  transformed = scale(transformed)

  # Perform Gram-Schmidt to ensure orthogonality of transformed phenotypes
  for (j in 2:5) {
    for (k in 1:(j-1)) {
      proj = sum(transformed[,j] * transformed[,k]) / sum(transformed[,k]^2)
      transformed[,j] = transformed[,j] - proj * transformed[,k]
    }
    # Normalize
    transformed[,j] = transformed[,j] / sqrt(sum(transformed[,j]^2))
  }

  # Write as new phenotype
  return(cbind(pheno[,1], pheno[,2], transformed))
}
