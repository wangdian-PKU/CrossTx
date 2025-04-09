#' Heatmap for DEG Analysis Based on TCGA High and Low Expression Groups
#'
#' This function generates a **heatmap** of log fold changes (logFC)
#' for differentially expressed genes (DEGs) across multiple datasets.
#'
#' @importFrom ComplexHeatmap Heatmap draw
#'
#' @param deg_data List. Output `processed data` from `prepare_DEG_heatmap_data()`, containing processed TCGA and mouse DEG data.
#' @param mouse_files Named list. The original input `mouse_files` from `merge_DEG_datasets()`, needed for column names.
#' @param col_names Character vector. Custom column names for the heatmap. Default uses `names(mouse_files)`.
#' @param cluster_rows Logical. Whether to cluster rows (default: FALSE).
#' @param cluster_cols Logical. Whether to cluster columns (default: FALSE).
#' @param color_palette Named vector. Colors for the heatmap (default: c("1" = "#ff7676", "-1" = "#66d4ff", "NA" = "white")).
#' @param width Integer. The width of pdf (default: 7.9).
#' @param height Integer. The height of pdf (default: 5.95).
#' @param output_path Character. File path for saving the plot (default: "./DEG/").
#'
#' @return A matrix (`merge_edger`) containing merged DEG expression data as the input object of `plot_DEG_barplot()`
#'
#' @examples
#' \dontrun{
#' # example usage
#' merge_edger <- plot_DEG_heatmap(deg_data = processed_data, mouse_files = mouse_files)
#' }
#'
#' @export
plot_DEG_heatmap <- function(deg_data,
                             mouse_files,
                             col_names = NULL, # Allow user-defined column names
                             cluster_rows = FALSE,
                             cluster_cols = FALSE,
                             color_palette = c("1" = "#ff7676", "-1" = "#66d4ff", "NA" = "white"),
                             width = 7.9,
                             height = 5.95,
                             output_path = "./DEG/") {
  # Ensure input is correctly formatted
  if (!is.list(deg_data) || !all(c("processed_tcga", "processed_mouse", "tcga_gene_list") %in% names(deg_data))) {
    stop("Error: Input data must be the output from `prepare_DEG_heatmap_data()`.")
  }

  # Extract processed mouse DEG data
  processed_mouse <- deg_data$processed_mouse

  # Merge mouse model data
  merge_edger <- do.call(cbind, lapply(processed_mouse, function(x) x[, 1]))

  # Ensure matrix formatting
  merge_edger <- as.matrix(merge_edger)

  # Restore row names using `tcga_gene_list`
  rownames(merge_edger) <- deg_data$tcga_gene_list

  # Assign column names
  if (is.null(col_names)) {
    colnames(merge_edger) <- names(mouse_files) # Default: use `names(mouse_files)`
  } else {
    if (length(col_names) != ncol(merge_edger)) {
      stop("Error: `col_names` length does not match number of columns in heatmap.")
    }
    colnames(merge_edger) <- col_names # User-defined column names
  }

  # Ensure heatmap is not empty
  if (nrow(merge_edger) == 0 || ncol(merge_edger) == 0) {
    stop("Error: Heatmap matrix is empty. Check input data formatting.")
  }

  # Generate output file name
  output_file <- paste0(output_path, "heatmap.pdf")

  # Save heatmap to file
  grDevices::pdf(output_file, width = width, height = height)

  # **Update legend labels**
  legend_labels <- c("Highly expression in tumor", "Low expression in tumor")
  names(legend_labels) <- c("1", "-1")

  # **Ensure Heatmap is drawn**
  heatmap_obj <- Heatmap(merge_edger,
    col = color_palette,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_cols,
    na_col = color_palette["NA"],
    show_row_names = FALSE,
    column_names_rot = 0,
    column_names_centered = TRUE,
    show_heatmap_legend = TRUE,
    heatmap_legend_param = list(title = "Expression", at = c(1, -1), labels = legend_labels)
  )

  draw(heatmap_obj) # Ensure the heatmap is actually rendered
  grDevices::dev.off()

  message("Heatmap saved to: ", output_file)

  # Return merged matrix for downstream analysis
  return(merge_edger)
}
