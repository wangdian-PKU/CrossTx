#' Gene Ontology (GO) Enrichment Analysis
#'
#' This function performs GO enrichment analysis (BP, MF, CC)
#' for differentially expressed genes (DEGs) obtained from edgeR analysis.
#'
#' @description
#' - Extracts significant DEGs (logFC > 1 & FDR < 0.05).
#' - Reads input data automatically from .tsv or .csv files.
#' - Converts gene IDs based on the selected species (hsa, mmu, rno).
#' - Supports user-defined ontology (BP, MF, CC).
#' - Performs GO enrichment analysis and saves the plot to a PDF file.
#'
#' @details
#' - Input Data: The function requires an edgeR DEG results file.
#' - Supported File Formats: .tsv (tab-separated) and .csv (comma-separated).
#' - Supported Species: hsa (human), mmu (mouse), rno (rat).
#' - Ontology Options (ont):
#'   - "BP": Biological Process (default)
#'   - "MF": Molecular Function
#'   - "CC": Cellular Component
#' - Visualization: GO term similarity is visualized using simplifyGO().
#'
#' @importFrom clusterProfiler enrichGO
#' @importFrom simplifyEnrichment simplifyGO GO_similarity
#' @importFrom dplyr filter mutate
#' @importFrom AnnotationDbi toTable
#' @importFrom tools file_ext
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom org.Mm.eg.db org.Mm.eg.db
#' @importFrom org.Rn.eg.db org.Rn.eg.db
#'
#' @param deg_file Character. Path to the edgeR DEG results file (.tsv or .csv).
#' @param species Character. One of hsa, mmu, or rno.
#' @param org_db Character. The corresponding `OrgDb` annotation database (e.g., "org.Hs.eg.db").
#' @param ont Character. Ontology type for GO enrichment ("BP", "MF", "CC"). Default is "BP".
#' @param column_title Character. The title used for GO similarity visualization.
#' @param width Integer. The width of PDF (default: 10.36).
#' @param height Integer. The height of PDF (default: 6.35).
#' @param output_path Character. Directory to save the PDF file (default: "./enrichment").
#'
#' @return Saves the GO enrichment analysis plot as a .pdf file.
#'
#' @examples
#' \dontrun{
#' # Example usage for human (TCGA data)
#' GO_enrichment_analysis(
#'   deg_file = "./data/TCGA_DEG.tsv",
#'   species = "hsa",
#'   org_db = "org.Hs.eg.db",
#'   ont = "BP",
#'   column_title = "TCGA GO BP terms"
#' )
#' }
#'
#' @export
GO_enrichment_analysis <- function(deg_file,
                                   species,
                                   org_db,
                                   ont = "BP",
                                   column_title,
                                   width = 10.36,
                                   height = 6.35,
                                   output_path = "./enrichment") {
  # Check if the "species" parameter is valid
  species_list <- c("hsa", "mmu", "rno")
  if (!species %in% species_list) {
    stop("Error: "species" must be one of 'hsa', 'mmu', or 'rno'.")
  }

  # Check if the "ont" parameter is valid
  if (!ont %in% c("BP", "MF", "CC")) {
    stop("Error: "ont" must be one of 'BP', 'MF', or 'CC'.")
  }

  # Automatically detect file format and read data
  file_ext <- tolower(file_ext(deg_file)) # Extract file extension
  if (file_ext == "tsv") {
    deg_data <- read.table(deg_file, header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE)
  } else if (file_ext == "csv") {
    deg_data <- read.csv(deg_file, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
  } else {
    stop("Error: Unsupported file format. Please provide a .tsv or .csv file.")
  }

  # Ensure required columns are present
  required_cols <- c("logFC", "FDR")
  if (!all(required_cols %in% colnames(deg_data))) {
    stop("Error: Input file must contain 'logFC' and 'FDR' columns.")
  }

  # Ensure row names are gene identifiers
  if (is.null(rownames(deg_data)) || any(rownames(deg_data) == "")) {
    stop("Error: The input file must have gene identifiers as row names.")
  }

  # Load species-specific gene ID mapping
  if (species == "hsa") {
    EG2Symbol <- toTable(org.Hs.egSYMBOL)
  } else if (species == "mmu") {
    EG2Symbol <- toTable(org.Mm.egSYMBOL)
  } else {
    EG2Symbol <- toTable(org.Rn.egSYMBOL)
  }

  # Filter significant DEGs (logFC > 1 & FDR < 0.05)
  geneLists <- deg_data %>%
    filter(abs(logFC) > 1 & FDR < 0.05) %>%
    mutate(symbol = rownames(.))

  # Convert gene IDs
  results <- merge(geneLists, EG2Symbol, by = "symbol", all.x = TRUE)
  results <- na.omit(results)
  id <- results$gene_id # Extract gene ID column

  # Perform GO enrichment analysis
  GO_results <- enrichGO(
    OrgDb = get(org_db), gene = id,
    ont = ont, pvalueCutoff = 0.05, readable = TRUE
  )

  GO_result_df <- as.data.frame(GO_results)

  go_terms <- na.omit(GO_result_df$ID)
  go_terms <- go_terms[grepl("^GO:\\d{7}$", go_terms)]

  if (length(go_terms) < 2) {
    stop("Error: The number of effective GO terms is too small (< 2) to conduct go similarity analysis. 
         Please check whether the input differential expression matrix is too small or the annotation fails.")
  }

  mat <- GO_similarity(go_terms, ont = ont)


  # Create output directory
  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
  }

  # Define output file name
  pdf_file <- file.path(output_path, paste0(column_title, ".pdf"))

  # Save GO enrichment analysis result to PDF
  grDevices::pdf(pdf_file, width = width, height = height)
  simplifyGO(mat, column_title = column_title)
  grDevices::dev.off()

  message("GO enrichment analysis saved to: ", pdf_file)
}
