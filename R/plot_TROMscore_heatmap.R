#' Plot TROM Score Heatmap for TCGA and Mouse Models
#'
#' This function visualizes transcriptional similarity scores between TCGA samples and mouse models.
#'
#' @details
#' - Returns 'heatmap_obj' (heatmap plot) and 'merge_long' (processed similarity scores).
#' - 'merge_long' is used as the input for 'plot_TROMscore_boxplot()'.
#'
#' @importFrom ComplexHeatmap Heatmap draw HeatmapAnnotation rowAnnotation
#' @importFrom circlize colorRamp2
#' @importFrom dplyr group_by summarise ungroup mutate select
#'
#' @param dm_trom_del_backgroud_clinic Data frame. Output from 'compute_similarity_matrix()', containing 'TROM scores'.
#' @param trom_matrix Data frame. 'Raw ws.trom similarity matrix' before background removal.Output from 'compute_similarity_matrix()'.
#' @param tcga_sample_size Integer. Total number of TCGA samples.
#' @param tcga_normal_count Integer. Number of TCGA normal samples.
#' @param model_group_count Integer. Total number of Mouse model samples.
#' @param condition Named list. User-defined sample counts per mouse model condition (e.g., 'list("HBV pten p53 ko" = 5, "DEN Cl4" = 21)').
#' @param clinical_vars Character vector. Clinical variables used in 'prepare_dm_trom_del_backgroud_clinic()'.
#' @param cnv_genes Character vector. Copy number variation (CNV) genes used in 'prepare_dm_trom_del_backgroud_clinic()'.
#' @param group_all The group information of TCGA and all mouse models from 'prepare_PCA_data()'.
#' @param width Integer. The width of PDF (default: 11.84).
#' @param height Integer. The height of PDF (default: 7.11).
#' @param output_path Character. File path for saving the heatmap (default: "./similarity/TROM_heatmap.pdf").
#' @param cluster_rows Logical. Whether to cluster rows (default: 'FALSE').
#' @param cluster_cols Logical. Whether to cluster columns (default: 'FALSE').
#'
#' @return A list containing:
#' \describe{
#'   \item{heatmap_obj}{The generated heatmap.}
#'   \item{merge_long}{Processed similarity scores for downstream analysis.}
#' }
#'
#' @examples
#' \dontrun{
#' heatmap_results <- plot_TROMscore_heatmap(
#'   dm_trom_del_backgroud_clinic = similarity_results$dm_trom_del_backgroud_clinic,
#'   trom_matrix = similarity_results$similarity_matrix,
#'   tcga_sample_size = 424,
#'   tcga_normal_count = 50,
#'   model_group_count = 52,
#'   condition = list("HBV pten p53 ko" = 5, "DEN Cl4" = 16),
#'   clinical_vars = c("Hepatitis B", "Non-Alcoholic"),
#'   cnv_genes = c("TP53", "PTEN"),
#'   group_all = pca_data$group_all,
#'   output_path = "./similarity/TROM_heatmap.pdf"
#' )
#' }
#'
#' @export
plot_TROMscore_heatmap <- function(dm_trom_del_backgroud_clinic,
                                   trom_matrix,
                                   tcga_sample_size,
                                   tcga_normal_count,
                                   model_group_count,
                                   condition,
                                   clinical_vars,
                                   cnv_genes,
                                   group_all,
                                   width = 11.84,
                                   height = 7.11,
                                   output_path = "./similarity/TROM_heatmap.pdf",
                                   cluster_rows = FALSE,
                                   cluster_cols = FALSE) {
  # Check input validity
  if (!is.data.frame(dm_trom_del_backgroud_clinic)) stop("Error: 'dm_trom_del_backgroud_clinic' must be a data frame.")
  if (!is.numeric(tcga_sample_size) || tcga_sample_size <= 0) stop("Error: 'tcga_sample_size' must be a positive integer.")
  if (!is.numeric(tcga_normal_count) || tcga_normal_count <= 0) stop("Error: 'tcga_normal_count' must be a positive integer.")
  if (nrow(dm_trom_del_backgroud_clinic) == 0) warning("Warning: 'dm_trom_del_backgroud_clinic' is empty. No data to plot.")

  # Step 1: Extract Similarity Matrix
  df <- dm_trom_del_backgroud_clinic[, 2:(model_group_count + 1)] %>% as.matrix()

  # Step 2: Define Color Mapping
  col_fun <- colorRamp2(c(0, 1, 5, 10, 20, 30), c("white", "#3f84a2", "#cae2a7", "#f8e4a1", "#d85b39", "#cb3025"))

  # Step 3: Process merge_long Data
  group_model <- group_all[(tcga_sample_size + 1):length(group_all)]

  merge_long <- data.frame()

  for (i in unique(group_model)) {
    if (grepl("normal", i)) {
      y <- trom_matrix[1:tcga_normal_count, group_model == i] %>% apply(2, mean)
      g1 <- rep("normal", length(y))
    } else {
      y <- trom_matrix[(tcga_normal_count + 1):(tcga_sample_size), group_model == i] %>% apply(2, mean)
      g1 <- rep("tumor", length(y))
    }
    x <- rep(i, length(y))

    data <- data.frame(x = x, y = y, group = g1)

    merge_long <- rbind(merge_long, data)
  }

  merge_long <- merge_long %>%
    mutate(condition = factor(rep(names(condition), unlist(condition)), levels = unique(names(condition)))) %>%
    mutate(x = factor(x, levels = unique(x)), group = factor(group, levels = c("normal", "tumor")))

  # Step 4: Compute Mean Scores
  merge_long_2 <- aggregate(merge_long$y, by = list(merge_long$x), FUN = mean)
  merge_long <- merge(merge_long, merge_long_2, by.x = "x", by.y = "Group.1") %>%
    mutate(x.y = round(x.y, 1))
  colnames(merge_long) <- c("x", "y", "group", "condition", "score")

  # Step 5: Define Col Annotations
  col_anno <- HeatmapAnnotation(
    model_tissue = merge_long$group,
    col = list(model_tissue = c("normal" = "#79c2f2", "tumor" = "#e5616e"))
  )

  # Step 6: Define Dynamic Row Annotations
  annotation_data <- list()
  color_mapping <- list()

  # TCGA tumor/normal annotation
  annotation_data[["TCGA_tissue"]] <- c(rep("tumor", tcga_sample_size - tcga_normal_count), rep("normal", tcga_normal_count))
  color_mapping[["TCGA_tissue"]] <- c("normal" = "#79c2f2", "tumor" = "#e5616e")

  # Iterate over clinical variables and CNV genes
  for (var in c(clinical_vars, cnv_genes)) {
    if (var %in% colnames(dm_trom_del_backgroud_clinic)) {
      annotation_data[[var]] <- rev(dm_trom_del_backgroud_clinic[[var]]) # Extract variable column

      # Define color mapping based on variable type
      if (var %in% clinical_vars) {
        color_mapping[[var]] <- c("Positive" = "#f6d593", "Negative" = "#ddeaac") # Default for clinical vars
      } else if (var %in% cnv_genes) {
        color_mapping[[var]] <- c("-2" = "#1010ff", "-1" = "#109eff", "0" = "#bababa", "1" = "#e274dd") # CNV colors
      }
    }
  }

  # Generate Row Annotations
  row_anno <- rowAnnotation(
    df = as.data.frame(annotation_data),
    col = color_mapping
  )

  # Step 7: Generate and Save Heatmap
  dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)

  grDevices::pdf(output_path, height = height, width = width)
  heatmap_obj <- Heatmap(df[rev(1:nrow(df)), ],
    col = col_fun,
    na_col = "white",
    cluster_rows = cluster_rows, cluster_columns = cluster_cols,
    row_names_side = "left",
    heatmap_legend_param = list(
      title = "Socre",
      legend_direction = "vertical"
    ),
    row_names_gp = gpar(fontsize = 10, font = 3),
    column_names_gp = gpar(fontsize = 10, font = 3),
    row_split = factor(rep(c("Tumor", "Normal"), c(tcga_sample_size - tcga_normal_count, tcga_normal_count))),
    column_split = factor(rep(names(condition), unlist(condition)), levels = unique(names(condition))),
    border = TRUE,
    row_gap = unit(0, "mm"), column_gap = unit(0, "mm"),
    show_column_names = FALSE,
    top_annotation = col_anno,
    left_annotation = row_anno
  )
  draw(heatmap_obj)

  message("TROMscore_heatmap saved to: ", output_path)


  grDevices::dev.off()

  return(list(heatmap_obj = heatmap_obj, merge_long = merge_long))
}
