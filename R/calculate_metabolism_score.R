#' Calculate Metabolism Signature Score
#'
#' This function computes the metabolism-related gene signature score using different methods.
#'
#' @importFrom IOBR calculate_sig_score
#'
#' @param eset_tpm Data frame. TPM-transformed RNA-seq expression matrix, output from `plot_TME_barplot()`.
#' @param method Character. Method for score calculation ("pca", "ssgsea", "zscore", "integration"). Default: "pca".
#' @param mini_gene_count Integer. Minimum gene count required per signature. Default: 2.
#'
#' @return A data frame containing metabolism scores per sample.
#'
#' @examples
#' \dontrun{
#' sig_meta <- calculate_metabolism_score(eset_tpm, method = "pca", mini_gene_count = 2)
#' }
#'
#' @export
calculate_metabolism_score <- function(eset_tpm, method = "pca", mini_gene_count = 2) {
  # Input validation
  if (!is.data.frame(eset_tpm)) stop("Error: 'eset_tpm' must be a data frame.")
  if (!method %in% c("pca", "ssgsea", "zscore", "integration")) {
    stop("Error: 'method' must be one of 'pca', 'ssgsea', 'zscore', or 'integration'.")
  }

  # Calculate signature score using IOBR::calculate_sig_score

  sig_meta <- calculate_sig_score(
    pdata = NULL,
    eset = eset_tpm,
    signature = signature_metabolism,
    method = method, # It can also be ssgsea or Z-score method
    mini_gene_count = mini_gene_count
  )

  return(sig_meta)
}
