#' Prepare TCGA Clinical and CNV Data for Similarity Analysis
#'
#' This function processes **TCGA clinical data** and **CNV (Copy Number Variation) data** to generate a
#' combined dataset for similarity analysis with mouse models.
#'
#' @description
#' - Extracts **user-specified clinical variables** (e.g., `Hepatitis`, `Non-Alcoholic`).
#' - Extracts **CNV information** for user-defined genes (e.g., `TP53`, `PTEN`).
#' - Returns a **merged dataset** containing **patient ID, clinical data, and CNV mutations**.
#'
#' @details
#' - **Input Data:** Requires clinical data (`.tsv`) and CNV data (`.tsv`).
#' - **User-defined Clinical Variables:** Extracted from risk factor history.
#' - **CNV Data:** Filters CNV values for selected genes.
#' - **Returns:** Merged data containing **clinical variables** and **CNV alterations**.
#'
#' @importFrom dplyr select mutate across
#'
#' @param clinic_file Character. Path to the TCGA clinical data file.
#' @param cnv_file Character. Path to the TCGA CNV data file.
#' @param clinical_vars Character vector. **User-defined clinical variables** to extract (e.g., `c("Hepatitis B", "Non-Alcoholic")`).
#' @param cnv_genes Character vector. **Genes to extract from CNV data** (e.g., `c("TP53", "PTEN")`).
#' @param risk_factor_column Character. **Column name in clinical data that contains risk factors** (e.g., `"HISTORY_HEPATO_CARCINOMA_RISK_FACTORS"`).
#'
#' @return A **data frame** containing:
#' - `PATIENT_ID`: Unique patient identifier.
#' - **Selected clinical variables** (e.g., HBV, Non-Alcoholic).
#' - **Selected CNV gene alterations** (e.g., TP53, PTEN).
#'
#' @examples
#' \dontrun{
#' merge_clinic_CNV <- prepare_similarity_data(
#'   clinic_file = "./lihc_tcga/data_bcr_clinical_data_patient.tsv",
#'   cnv_file = "./lihc_tcga/data_CNA.tsv",
#'   clinical_vars = c("Hepatitis B", "Non-Alcoholic"),
#'   cnv_genes = c("TP53", "PTEN"),
#'   risk_factor_column = "HISTORY_HEPATO_CARCINOMA_RISK_FACTORS"
#' )
#' }
#'
#' @export
prepare_similarity_data <- function(clinic_file,
                                    cnv_file,
                                    clinical_vars,
                                    cnv_genes,
                                    risk_factor_column) {
  # **Step 1: Read TCGA Clinical Data**
  clinic <- read.table(clinic_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  clinic$PATIENT_ID <- gsub("-", ".", clinic$PATIENT_ID) # Standardize PATIENT_ID format

  # Extract essential clinical variables
  clinic_select <- clinic %>%
    select(PATIENT_ID)

  # **Dynamically Extract User-Defined Clinical Variables**
  for (var in clinical_vars) {
    clinic_select[[var]] <- ifelse(grepl(var, clinic[[risk_factor_column]], ignore.case = TRUE), var, NA)
  }

  clinic_select[clinic_select == "[Not Available]"] <- NA # Replace unavailable values with NA

  # **Step 2: Read TCGA CNV Data**
  cnv_data <- read.table(cnv_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)


  # Extract relevant CNV genes
  cnv_selected <- cnv_data[cnv_genes, , drop = FALSE]

  # **Step 3: Merge Clinical and CNV Data**
  merge_clinic_CNV <- rbind(substring(colnames(cnv_data), 1, 12), cnv_selected) %>%
    t() %>%
    merge(., clinic_select, by.x = "1", by.y = "PATIENT_ID", all.y = TRUE) %>%
    mutate(across(all_of(cnv_genes), ~ factor(., levels = c("-2", "-1", "0", "1"))))

  colnames(merge_clinic_CNV)[1] <- "PATIENT_ID"

  return(merge_clinic_CNV)
}
