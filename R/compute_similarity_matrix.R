#' Compute Similarity Matrix Between TCGA and Mouse Model Data
#'
#' This function calculates a **transcriptional similarity matrix** between **TCGA samples** and **mouse models**.
#'
#' @description
#' - **Formats** the input expression matrix by setting gene IDs as the first column.
#' - Extracts **TCGA and mouse model expression data** from a preprocessed matrix.
#' - Computes a **similarity matrix** using `ws.trom()`, comparing TCGA with mouse models.
#' - Filters similarity data **to retain only relevant comparisons** (tumor vs. tumor, normal vs. normal).
#' - **Incorporates TCGA clinical data** by merging with `merge_clinic_CNV`.
#'
#' @importFrom dplyr mutate across all_of arrange
#' @importFrom TROM ws.trom
#'
#' @param expression_matrix The Batch-corrected expression data with genes as row names from `prepare_PCA_data()`.
#' @param group_all The group information of TCGA and all mouse models from `prepare_PCA_data()`.
#' @param tcga_sample_size Integer. **Number of TCGA samples** in the dataset.
#' @param tcga_normal_count Integer. **Number of TCGA normal samples**.
#' @param merge_clinic_CNV Data frame. Processed TCGA **clinical + CNV** data from `prepare_similarity_data()`.
#' @param cnv_genes Character vector. Gene names for CNV annotations (e.g., `c("TP53", "PTEN")`).
#' @param clinical_vars Character vector. **User-defined clinical variables** to extract (e.g., `c("Hepatitis B", "Non-Alcoholic")`).
#'
#' @return A list containing:
#' \describe{
#'   \item{similarity_matrix}{The full transcriptional similarity matrix before background removal.}
#'   \item{dm_trom_del_backgroud_clinic}{The filtered similarity matrix merged with clinical data.}
#' }
#'
#' @examples
#' \dontrun{
#' similarity_results <- compute_similarity_matrix(
#'   expression_matrix = pca_data$corrected_data,
#'   group_all = pca_data$group_all,
#'   tcga_sample_size = 424,
#'   tcga_normal_count = 50,
#'   merge_clinic_CNV = merge_clinic_CNV,
#'   cnv_genes = c("TP53", "PTEN"),
#'   clinical_vars = c("Hepatitis B", "Non-Alcoholic")
#' )
#' }
#'
#' @export
compute_similarity_matrix <- function(expression_matrix,
                                      group_all,
                                      tcga_sample_size,
                                      tcga_normal_count,
                                      merge_clinic_CNV,
                                      cnv_genes,
                                      clinical_vars) {
  # **Step 1: Prepare Expression Data**
  newData <- as.data.frame(expression_matrix)
  newData$gene_id <- rownames(newData) # Store row names as a separate column
  newData <- newData[, c(ncol(newData), 1:(ncol(newData) - 1))] # Move gene_id column to the first position
  colnames(newData)[1] <- "gene_id"

  # **Step 2: Extract TCGA and Mouse Model Data**
  tcga_data <- newData[, 1:(tcga_sample_size + 1)] # Extract TCGA data
  model_data <- newData[, c(1, (tcga_sample_size + 2):ncol(newData))] # Extract mouse model data

  # **Step 3: Compute Similarity Matrix**
  similarity_matrix <- ws.trom(
    sp_gene_expr = tcga_data, single = FALSE,
    sp_gene_expr2 = model_data, provide = FALSE,
    save_overlap_genes = FALSE
  )

  # **Step 4: Define Filtering Function**
  del_background <- function(similarity_matrix, tissue1, tissue2) {
    tissue1_index <- grepl(tissue1, colnames(similarity_matrix))
    tissue2_index <- grepl(tissue2, colnames(similarity_matrix))

    similarity_matrix[grepl(tissue2, group_all[1:tcga_sample_size]), tissue1_index] <- NA
    similarity_matrix[grepl(tissue1, group_all[1:tcga_sample_size]), tissue2_index] <- NA

    return(similarity_matrix) # Retain only tumor-tumor and normal-normal comparisons
  }

  # **Step 5: Apply Filtering**
  filtered_similarity_matrix <- del_background(similarity_matrix, "tumor", "normal") %>%
    as.data.frame()

  # **Step 6: Incorporate Clinical Data**
  sort_columns <- c(clinical_vars, cnv_genes)
  dm_trom_del_backgroud_clinic <- filtered_similarity_matrix %>%
    mutate(
      PATIENT_ID = substr(rownames(.), 1, 12),
      tissue = c(
        rep("normal", tcga_normal_count),
        rep("tumor", tcga_sample_size - tcga_normal_count)
      )
    ) %>%
    merge(., merge_clinic_CNV, by = "PATIENT_ID") %>%
    arrange(tissue, across(all_of(sort_columns)))

  return(list(
    similarity_matrix = similarity_matrix, # Full similarity matrix before filtering
    dm_trom_del_backgroud_clinic = dm_trom_del_backgroud_clinic # Processed similarity matrix with clinical data
  ))
}
