#' Prepare PCA Data for TCGA and Mouse Models
#'
#' This function processes gene expression data from TCGA and multiple mouse models
#' to generate a merged expression matrix suitable for PCA analysis.
#'
#' @description
#' - Reads expression data from TCGA and multiple mouse cancer models.
#' - Ensures gene expression data are standardized and merged based on homologous genes.
#' - Allows optional batch correction using 'ComBat'.
#'
#' @details
#' - Uses gene names as row identifiers to merge datasets.
#' - Sample grouping ('group_all') and batch effects ('batch') are dynamically assigned.
#' - If 'batch_correction = TRUE', batch correction is applied using 'ComBat'.
#'
#' @importFrom sva ComBat
#' @importFrom tools file_ext
#'
#' @param tcga_file Character. Path to TCGA expression data file.
#' @param mouse_files Named list of file paths for mouse model expression data.
#' @param sample_counts Named list specifying the number of normal and tumor samples for each dataset.
#'                      Each dataset must have a 'normal' and 'tumor' entry, e.g.,
#'                      'list("TCGA" = c(normal = 50, tumor = 374), "Our_Model" = c(normal = 3, tumor = 2))'
#' @param batch_correction Logical. Whether to apply batch correction (default: 'TRUE').
#'
#' @return A list containing:
#' \describe{
#'   \item{merged_data}{Merged expression matrix (genes × samples).}
#'   \item{group_all}{Sample group labels (normal/tumor for each dataset).}
#'   \item{batch}{Batch assignment for each sample (dataset source).}
#'   \item{corrected_data}{Batch-corrected expression matrix (if 'batch_correction = TRUE').}
#' }
#'
#' @examples
#' \dontrun{
#' mouse_files <- list(
#'   "Our_Model" = "./data/Our_Model.homo.csv",
#'   "GSE172629" = "./data/GSE172629.homo.tsv"
#' )
#' sample_counts <- list(
#'   "TCGA" = c(normal = 50, tumor = 374),
#'   "Our_Model" = c(normal = 3, tumor = 2),
#'   "GSE172629" = c(normal = 3, tumor = 3)
#' )
#' pca_data <- prepare_PCA_data(
#'   tcga_file = "./data/TCGA.tsv",
#'   mouse_files = mouse_files,
#'   sample_counts = sample_counts,
#'   batch_correction = TRUE
#' )
#' }
#'
#' @export
prepare_PCA_data <- function(tcga_file, mouse_files, sample_counts, batch_correction = TRUE) {
  # Helper function to auto-detect file type
  read_expression_file <- function(file) {
    ext <- tolower(tools::file_ext(file))
    if (ext == "csv") {
      read.csv(file, stringsAsFactors = FALSE)
    } else if (ext == "tsv") {
      read.table(file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    } else {
      stop("Unsupported file format: must be .csv or .tsv")
    }
  }

  # Read TCGA data
  tcga_data <- read_expression_file(tcga_file)

  # Read mouse model data
  mouse_list <- lapply(mouse_files, read_expression_file)

  # Merge all datasets using gene names
  merged_data <- tcga_data
  for (i in mouse_list) {
    merged_data <- merge(merged_data, i, by = "X")
  }
  rownames(merged_data) <- merged_data$X
  merged_data <- merged_data[, -1] # Remove redundant index column

  # Create sample group labels ('group_all')
  group_all <- unlist(lapply(names(sample_counts), function(dataset) {
    normal_count <- sample_counts[[dataset]]["normal"]
    tumor_count <- sample_counts[[dataset]]["tumor"]
    c(rep(paste0(dataset, "_normal"), normal_count), rep(paste0(dataset, "_tumor"), tumor_count))
  }))

  # Generate sample names ('colnames(merged_data)')
  tcga_sample_count <- sample_counts[["TCGA"]]["normal"] + sample_counts[["TCGA"]]["tumor"]
  tcga_colnames <- colnames(merged_data)[1:tcga_sample_count] # Preserve TCGA column names

  other_sample_names <- unlist(lapply(names(sample_counts)[-1], function(dataset) {
    normal_count <- sample_counts[[dataset]]["normal"]
    tumor_count <- sample_counts[[dataset]]["tumor"]
    c(
      paste0(rep(paste0(dataset, "_normal"), normal_count), c(1:normal_count)),
      paste0(rep(paste0(dataset, "_tumor"), tumor_count), c(1:tumor_count))
    )
  }))

  if ((length(tcga_colnames) + length(other_sample_names)) != ncol(merged_data)) {
    stop("Error: Mismatch between calculated sample names and merged data columns.")
  }
  colnames(merged_data) <- c(tcga_colnames, other_sample_names)

  # Assign batch effects
  batch <- unlist(lapply(names(sample_counts), function(dataset) {
    rep(dataset, sum(sample_counts[[dataset]]))
  }))

  # Batch correction using 'ComBat'
  corrected_data <- if (batch_correction) {
    ComBat(dat = as.matrix(merged_data), batch = batch)
  } else {
    NULL
  }

  return(list(
    merged_data = merged_data,
    group_all = group_all,
    batch = batch,
    corrected_data = corrected_data
  ))
}
