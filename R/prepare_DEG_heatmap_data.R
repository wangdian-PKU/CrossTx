#' Prepare DEG Data for Heatmap Visualization
#'
#' This function preprocesses TCGA and mouse model DEG data by filtering significant genes,
#' performing homologene conversion, and formatting logFC values.
#'
#' @description
#' This function selects significant DEGs from TCGA data, applies homologene conversion,
#' and standardizes logFC values for visualization in a heatmap.
#'
#' @details
#' - Filters genes with `logFC > 1` and `FDR < 0.05` in TCGA.
#' - Performs homologene conversion to retain only mouse-human homologous genes.
#' - Formats logFC values:
#'   - `1` for TCGA high-expression genes.
#'   - `-1` for TCGA low-expression genes.
#'   - `NA` for non-significant genes.
#'
#' @importFrom dplyr filter arrange mutate
#' @importFrom homologene homologene
#'
#' @param tcga_file Character. Path to the TCGA DEG file.
#' @param mouse_files Named list of file paths for mouse model DEG files.
#' @param inTax Numeric. Input taxonomy ID (default: 9606 for human).
#' @param outTax Numeric. Output taxonomy ID (default: 10090 for mouse).
#'
#' @return A list containing:
#'   - `processed_tcga`: Processed TCGA DEG matrix.
#'   - `processed_mouse`: Named list of processed mouse DEG matrices.
#'   - `tcga_gene_list`: List of homologous genes used in processing.
#'
#' @examples
#' \dontrun{
#' # Example usage
#' tcga_path <- "./data/TCGA_DEG.tsv"
#' mouse_files <- list(
#'   "Our_Model" = "./data/Our_Model_DEG.tsv",
#'   "GSE172629" = "./data/GSE172629_DEG.tsv"
#' )
#'
#' # Run data processing
#' processed_data <- prepare_DEG_heatmap_data(tcga_path, mouse_files)
#'
#' # View processed TCGA data
#' head(processed_data$processed_tcga)
#'
#' # View processed mouse data
#' head(processed_data$processed_mouse$Our_Model)
#' }
#'
#' @export
prepare_DEG_heatmap_data <- function(tcga_file, mouse_files, inTax = 9606, outTax = 10090) {
  # **Read TCGA DEG data**
  tcga_data <- read.table(tcga_file, header = TRUE, sep = "\t", row.names = 1)[, c(1, 4)]
  # **Filter genes with logFC > 1 & FDR < 0.05**
  tcga_filtered <- tcga_data %>%
    filter(abs(logFC) >= 1 & FDR < 0.05) %>%
    arrange(desc(logFC))
  # **Get homologous genes**
  homolog_genes <- homologene(rownames(tcga_filtered), inTax = inTax, outTax = outTax)[, c(1:2)]
  # **Keep only homologous genes**
  tcga_filtered <- tcga_filtered[rownames(tcga_filtered) %in% homolog_genes[, 1], ]
  # **Count up-regulated and down-regulated genes**
  pos_count <- sum(tcga_filtered$logFC > 0)
  neg_count <- sum(tcga_filtered$logFC < 0)
  # **Define change_colnames()**
  change_colnames <- function(data, name) {
    colnames(data) <- paste(name, colnames(data), sep = "_")
    # **Filter non-significant genes**
    data[which(abs(data[, 1]) < 1 & data[, 2] < 0.05), 1] <- NA
    # **Up-regulated genes (logFC > 1)**
    data[1:pos_count, 1][which(data[1:pos_count, 1] > 1)] <- 1
    data[1:pos_count, 1][which(data[1:pos_count, 1] < 1)] <- NA
    # **Down-regulated genes (logFC < -0.5)**
    data[(pos_count + 1):(pos_count + neg_count), 1][which(data[(pos_count + 1):(pos_count + neg_count), 1] < -0.5)] <- -1
    data[(pos_count + 1):(pos_count + neg_count), 1][which(data[(pos_count + 1):(pos_count + neg_count), 1] > -0.5)] <- NA
    return(data)
  }
  # **Process mouse data**
  processed_mouse <- lapply(names(mouse_files), function(name) {
    mouse_data <- read.table(mouse_files[[name]], header = TRUE, sep = "\t", row.names = 1)
    # **Filter strictly according to TCGA row names**
    mouse_data <- mouse_data[rownames(tcga_filtered), c(1, 4)] # Only keep logFC (column 1) and FDR (column 4)
    # **Call change_colnames() for processing**
    mouse_data <- change_colnames(mouse_data, name)
    return(mouse_data)
  })
  # **Return processed data**
  return(list(
    processed_tcga = tcga_filtered,
    processed_mouse = processed_mouse,
    tcga_gene_list = rownames(tcga_filtered)
  ))
}
