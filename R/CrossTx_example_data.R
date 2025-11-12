#' Example input dataset for the CrossTx transcriptome similarity workflow
#'
#' A deeply nested list containing all necessary input data to demonstrate the full
#' transcriptome similarity analysis pipeline implemented in the CrossTx package.
#' This includes expression, DEG, homolog, and clinical datasets across three cancer types:
#' COAD, LUAD, and BRCA, each with TCGA data and 2–3 matched mouse models.
#'
#' @format A nested list of 3 cancer types:
#' \describe{
#'   \item{COAD}{Contains \code{TCGA} and three GSE mouse models: \code{GSE245724}, \code{GSE263502}, \code{GSE285835}}
#'   \item{LUAD}{Same structure as COAD}
#'   \item{BRCA}{Same structure as COAD}
#'   \item{Sub-items (per model)}{
#'     Each cancer type includes the following named elements:
#'     \describe{
#'       \item{cleancounts}{Normalized expression matrix (genes × samples)}
#'       \item{edgerout}{Differentially expressed genes computed via edgeR}
#'       \item{homo}{Homologous gene mappings between human and mouse}
#'       \item{clinical_unique}{TCGA-only: Clinical annotation table}
#'       \item{CNV}{TCGA-only: Copy number variation matrix}
#'     }
#'   }
#' }
#'
#' @details This dataset is intended for demonstration and reproducibility
#' purposes. All example files are preprocessed and provided in .tsv format
#' within the nested list structure. Users can directly load and use
#' this object to test all functions in the CrossTx package.
#'
#' @usage data(CrossTx_example_data)
#'
#' @examples
#' data(CrossTx_example_data)
#'
#' # Show available cancer types
#' names(CrossTx_example_data)
#'
#' # Show TCGA and mouse models under COAD
#' names(CrossTx_example_data$COAD)
#'
#' # Inspect edgeR DEG results from TCGA-COAD
#' head(CrossTx_example_data$COAD$TCGA$edgerout)
#'
#' # Inspect cleancounts from LUAD-GSE263502
#' head(CrossTx_example_data$LUAD$GSE263502$cleancounts)
#'
#' @source Curated and formatted by Dian Wang for this package
"CrossTx_example_data"
