#' Compute transformed phenotype based on LDAK output.
#'
#' The transformation are used in combination with the original
#' covariate to compute a smoothed phenotype for association analysis.
#'
#' @param pheno Phenotype matrix or dataframe in PLINK format
#' @param cov Covariate vector
#' @param size Number of subgroups for stratifying the phenotype
#' @param trans Transformation matrix obtained from LDAK
#' @return The output of the system call
#' @export
transform <- function(pheno = NULL, cov = NULL, size = numeric(), trans = NULL) {

  # Compute subset
  cov_cases = cov[which(pheno[,3] == 1)]

  # Compute new phenotype by smoothing transformed phenotype
  transformed <- matrix(NA, dim(pheno)[1], 5)
  quants_medians <- (2 * (1:size)-1)/(2*size)

  for (i in 1:5) {

    # Fit a quadratic function
    y <- as.numeric(transformation[i,])
    fit <- lm(y ~ quants_medians + I(quants_medians^2))

    # Map covariate values to quantiles, and compute predicted transformation
    quantiles <- (2 * rank(cov_cases) - 1) / (2 * max(rank(cov_cases)))
    y_pred <- predict(fit, newdata = data.frame(quants_medians = quantiles))

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
