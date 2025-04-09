#' Tumor Microenvironment (TME) Barplot Visualization
#'
#' This function performs immune infiltration analysis using various deconvolution methods
#' and visualizes the immune cell composition in a stacked barplot.
#'
#' @description
#' - Converts count matrix to TPM using `IOBR::count2tpm()`, or allows direct TPM input.
#' - Supports multiple deconvolution methods: `"cibersort"`, `"mcpcounter"`, `"xcell"`, `"quantiseq"`.
#' - Generates a stacked barplot to visualize immune infiltration across sample groups.
#'
#' @details
#' - **Input Data:** Users can provide either a raw count matrix (`input_type = "count"`) or TPM data (`input_type = "TPM"`).
#' - **Group Definition:** Users must provide a valid `group_as` vector matching the sample order.
#' - **Visualization:** Saves a PDF file showing the fraction of immune cell types across groups.
#'
#' @importFrom IOBR count2tpm deconvo_tme cell_bar_plot
#' @importFrom dplyr mutate group_by
#'
#' @param data_matrix Data frame. `Raw` count matrix or `TPM` data (genes × samples).
#' @param input_type Character. `"count"` (convert count to TPM) or `"TPM"` (use TPM directly). Default is `"count"`.
#' @param method Character. The immune infiltration method (`"cibersort"`, `"mcpcounter"`, `"xcell"`, `"quantiseq"`).
#' @param group_as Character vector. Clarify the specific names of different treatment groups and the respective numbers of normal and cancer groups in each of them.
#' @param output_path Character. Directory to save the plot (default: `"./TME"`).
#' @param height Numeric. Height of the output PDF (default: `5.2`).
#' @param width Numeric. Width of the output PDF (default: `10`).
#'
#' @return A list containing:
#' - **`eset_tpm`**: The TPM-transformed expression matrix.
#' - **`method_model`**: Immune infiltration results.
#' - **A stacked barplot PDF file**: Saved at `output_path`, showing immune cell fractions per group.
#'
#' @examples
#' \dontrun{
#' # Example usage (using raw count matrix)
#' group_as <- c(rep("TCGA_normal", 50), rep("TCGA_tumor", 374), rep("HBV Pten KO_normal", 3), rep("HBV Pten KO_tumor", 3))
#'
#' # TME barplot
#' TME_barplot_result <- plot_TME_barplot(data_matrix = merge.data, input_type = "count", method = "cibersort", group_as = group_as)
#'
#' # Extract TPM matrix for further analysis
#' eset_tpm <- TME_barplot_result$eset_tpm
#' }
#'
#' @export
plot_TME_barplot <- function(data_matrix,
                             input_type = "count",
                             method,
                             group_as,
                             output_path = "./TME",
                             height = 5.2,
                             width = 10) {
  # **Check input data format**
  if (!is.matrix(data_matrix) && !is.data.frame(data_matrix)) {
    stop("Error: `data_matrix` must be a data frame or matrix.")
  }

  if (!is.character(group_as) || length(group_as) != ncol(data_matrix)) {
    stop("Error: `group_as` must be a character vector with the same length as the number of samples.")
  }

  if (!method %in% c("cibersort", "mcpcounter", "xcell", "quantiseq")) {
    stop("Error: `method` must be one of 'cibersort', 'mcpcounter', 'xcell', or 'quantiseq'.")
  }

  if (!input_type %in% c("count", "TPM")) {
    stop("Error: `input_type` must be either 'count' or 'TPM'.")
  }

  # **Convert count matrix to TPM (if required)**
  if (input_type == "count") {
    message("Converting count matrix to TPM...")
    eset_tpm <- count2tpm(countMat = data_matrix, idType = "SYMBOL", org = "hsa")
  } else {
    message("Using provided TPM matrix...")
    eset_tpm <- as.matrix(data_matrix)
  }

  # **Perform immune infiltration analysis**
  method_model <- switch(method,
    cibersort = deconvo_tme(eset = eset_tpm, method = "cibersort", arrays = FALSE, perm = 200),
    mcpcounter = deconvo_tme(eset = eset_tpm, method = "mcpcounter"),
    xcell = deconvo_tme(eset = eset_tpm, method = "xcell", arrays = FALSE),
    quantiseq = deconvo_tme(eset = eset_tpm, method = "quantiseq", arrays = FALSE, scale_mrna = TRUE)
  )

  # **Compute group averages**
  method_model_merge <- aggregate(method_model[, -1], by = list(group_as), FUN = mean)
  colnames(method_model_merge)[1] <- "ID"

  # **Ensure correct factor ordering**
  method_model_merge$ID <- factor(method_model_merge$ID, levels = unique(group_as))

  # **Create output directory**
  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
  }

  # **Generate PDF output**
  pdf_file <- file.path(output_path, paste0(method, "_cell_bar_plot.pdf"))
  grDevices::pdf(pdf_file, height = height, width = width)

  cell_bar_plot(
    input = method_model_merge,
    title = paste0(method, "Fraction"),
    pattern = "_cell"
  )

  grDevices::dev.off()

  # **Save immune infiltration data**
  csv_file <- file.path(output_path, paste0(method, "_result.csv"))
  utils::write.csv(method_model, csv_file)

  message("TME barplot saved to: ", pdf_file)
  message("TME results saved to: ", csv_file)

  # **Return TPM matrix and immune infiltration data**
  return(list(eset_tpm = eset_tpm, method_model = method_model))
}
