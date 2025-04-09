#' Barplot for DEG Count Comparison
#'
#' This function creates a **barplot** comparing the number of upregulated and
#' downregulated genes across multiple mouse HCC models.
#'
#' @description
#' This function visualizes the **distribution of DEGs** in different mouse
#' HCC models relative to TCGA-based gene expression categories.
#'
#' @details
#' - **Input Data:** The function requires `merge_edger` from `plot_DEG_heatmap()`.
#' - **Barplot Groups:**
#'   - **X-axis:** Mouse model datasets (contrast).
#'   - **Y-axis:** Number of genes.
#'   - **Colors:**
#'     - **Red (`up`)**: Genes following TCGA high-expression pattern.
#'     - **Blue (`down`)**: Genes following TCGA low-expression pattern.
#'
#' @importFrom ggplot2 ggplot geom_bar aes scale_fill_manual xlab ylab theme_minimal theme element_text
#' @importFrom tidyr pivot_longer
#'
#' @param heatmap_data Matrix. The `merge_edger` output from `plot_DEG_heatmap()`.
#' @param colors Named list. Colors for upregulated (`"up"`) and downregulated (`"down"`) genes.
#'               Default: `list(up = "#ff7676", down = "#66d4ff")`.
#' @param width Integer. The width of ggplot_draw (default: 8).
#' @param height Integer. The height of ggplot_draw (default: 4).
#' @param output_path Character. File path for saving the plot (default: "./DEG/").
#'
#' @return A `ggplot` object representing the barplot visualization.
#'
#' @examples
#' \dontrun{
#' # Generate barplot
#' barplot <- plot_DEG_barplot(heatmap_data = merge_edger)
#' print(barplot)
#' }
#'
#' @export
plot_DEG_barplot <- function(heatmap_data,
                             colors = list(up = "#ff7676", down = "#66d4ff"),
                             width = 8,
                             height = 4,
                             output_path = "./DEG/") {
  # Ensure input data is a matrix from `plot_DEG_heatmap()`
  if (!is.matrix(heatmap_data)) {
    stop("Error: `heatmap_data` must be a matrix from `plot_DEG_heatmap()`.")
  }

  # **Transform data format**
  merge_edger_long <- heatmap_data %>%
    as.data.frame() %>%
    pivot_longer(cols = everything(), names_to = "contrast", values_to = "foldchange") %>%
    na.omit() %>%
    mutate(
      foldchange = factor(ifelse(foldchange == -1, "down", "up"), levels = c("up", "down")),
      contrast = factor(contrast, levels = colnames(heatmap_data)) # Ensure contrast order
    )

  # **Generate barplot**
  p <- ggplot(merge_edger_long) +
    geom_bar(aes(x = contrast, fill = foldchange)) +
    scale_fill_manual(name = "", values = colors) +
    xlab("") +
    ylab("Number of genes") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # Save plot in the specified format
  file_ext <- paste0(output_path, "barplot.pdf")
  ggsave(file_ext, plot = p, width = width, height = height)
  message("Barplot saved to: ", file_ext)


  return(p)
}
