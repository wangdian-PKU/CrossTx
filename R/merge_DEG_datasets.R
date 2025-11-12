#' Merge Differential Expression Gene (DEG) Datasets
#'
#' This function loads multiple DEG datasets (both human and mouse models),
#' merges them into a single standardized format, and annotates each dataset
#' with a group label for comparison.
#'
#' @description
#' The function takes in:
#' 1. A DEG result file from TCGA (human cancer dataset)
#' 2. A list of DEG result files from various mouse cancer models (GEO data or user-provided data)
#'
#' It merges all datasets into a single `merge.diff` dataframe with standardized
#' column names and labels, making it suitable for visualization and further analysis.
#'
#' @details
#' - Supported file formats: .tsv (tab-separated) and .csv (comma-separated).
#' - Required columns: 'gene', 'logFC', 'FDR'.
#' - If a file format is not supported, an error message will be displayed.
#' - DEG threshold: You can customize 'logFC' and 'FDR' thresholds (default: logFC > 1, FDR < 0.05).
#'
#' @importFrom dplyr mutate bind_rows case_when
#' @importFrom tools file_ext
#'
#' @param human_file_path Character. Path to the TCGA human cancer DEG .tsv or .csv file.
#' @param mouse_files Named list of file paths to mouse cancer model DEG .tsv or .csv files.
#'                    Each element should be named according to the dataset source.
#' @param logFC_value Numeric. Threshold for log2 fold change (default: 1).
#' @param FDR_value Numeric. Threshold for adjusted p-value (default: 0.05).
#'
#' @return A merged data frame (`merge.diff`) with the following columns:
#' \describe{
#'   \item{gene}{Gene identifier}
#'   \item{logFC}{Log fold change of gene expression}
#'   \item{FDR}{False discovery rate (adjusted p-value)}
#'   \item{contrast}{Dataset source (e.g., "TCGA", "Our", "GSE172629", etc.)}
#'   \item{change}{Gene expression change status: "Up", "Down", or "Stable"}
#' }
#'
#' @examples
#' \dontrun{
#' # Example file paths (replace with actual file paths)
#' human_tcga_path <- "./data/TCGA_DEG.tsv"
#' mouse_files <- list(
#'   "Our_Model" = "./data/Our_Model_DEG.csv", # CSV format
#'   "GSE172629" = "./data/GSE172629_DEG.tsv", # TSV format
#'   "GSE208279" = "./data/GSE208279_DEG.tsv"
#' )
#'
#' # Load and merge DEG data
#' merged_data <- merge_DEG_datasets(human_tcga_path, mouse_files)
#'
#' # Display the first few rows
#' head(merged_data)
#' }
#'
#' @export
merge_DEG_datasets <- function(human_file_path, mouse_files, logFC_value = 1, FDR_value = 0.05) {
  # read files
  read_deg_file <- function(file_path) {
    file_ext <- tolower(file_ext(file_path))

    if (file_ext == "tsv") {
      data <- read.table(file_path, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
    } else if (file_ext == "csv") {
      data <- read.csv(file_path, header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
    } else {
      stop(paste("Error: Unsupported file format for", file_path, "Only .tsv and .csv files are allowed."))
    }

    return(data)
  }

  # read TCGA DEG data
  tcga_data <- read_deg_file(human_file_path) %>%
    mutate(contrast = "TCGA", gene = rownames(.))

  # read mouse_model DEG data
  mouse_list <- lapply(names(mouse_files), function(name) {
    read_deg_file(mouse_files[[name]]) %>%
      mutate(contrast = name, gene = rownames(.))
  })

  # merge data
  merge.diff <- bind_rows(tcga_data, bind_rows(mouse_list)) %>%
    mutate(
      change = case_when(
        logFC > logFC_value & FDR < FDR_value ~ "Up",
        logFC < -logFC_value & FDR < FDR_value ~ "Down",
        TRUE ~ "Stable"
      ),
      contrast = factor(contrast, levels = c("TCGA", names(mouse_files)))
    )

  return(merge.diff)
}
