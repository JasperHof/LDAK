#' Prepare a stratified phenotype for multivariate LDAK.
#'
#' The binary outcome is stratified based on quantiles of the
#' associated covariate.
#'
#' @param pheno Phenotype matrix or dataframe in PLINK format
#' @param cov Covariate vector in PLINK format
#' @param size Number of subgroups for stratifying the phenotype
#' @return The output of the system call
#' @export
create_subtype <- function(pheno = NULL, cov = NULL, size = 5) {
  if (dim(pheno)[2] != length(cov)) {
    stop("Covariate vector should be of same length as phenotype")
  }


}
