#' Prepare a stratified phenotype for multivariate LDAK.
#'
#' The binary outcome is stratified based on quantiles of the
#' associated covariate.
#'
#' @param pheno Phenotype matrix or dataframe in PLINK format
#' @param cov Covariate vector
#' @param size Number of subgroups for stratifying the phenotype
#' @return The output of the system call
#' @export
make_subtype <- function(pheno = NULL, cov = NULL, size = numeric()) {

  if ( dim(pheno)[1] != length(cov) ) {
    stop("Covariate should be a vector, whose length matches the number of individuals in phenotype")
  }
  if ( !all(pheno[,3] %in% c(0,1,NA)) ) {
    stop("Phenotypes should be 0, 1, or NA (missing)")
  }

  # Compute quantiles of cases - restrict to unique covariate values
  quants = quantile(unique(as.numeric(cov[which(pheno[,3] == 1)])), (1:(size-1))/size, na.rm = T)

  quants_split <- cut(as.numeric(cov),
                    breaks = c(-Inf, quants, Inf),  # use -Inf and Inf to cover all values
                    labels = FALSE,                 # numeric codes
                    right = TRUE)
  quants_group <- sapply(1:size, function(i) as.numeric(quants_split == i))

  # Set controls to zero and cases of other subgroup to missing
  quants_group[which(pheno[,3] == 0),] = 0
  quants_group[quants_group == 0] = NA
  quants_group[rowSums(is.na(quants_group)) == size,] = 0

  # Write as new phenotype
  return(cbind(pheno[,1], pheno[,2], quants_group))
}
