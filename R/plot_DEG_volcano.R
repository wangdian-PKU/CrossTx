#' Volcano Plot for DEG Analysis
#'
#' Creates a volcano plot for visualizing differential gene expression (DEG) data.
#'
#' @description
#' This function generates a **volcano plot** that displays the log2 fold change (logFC)
#' against the -log10 false discovery rate (FDR) of differentially expressed genes.
#'
#' @details
#' - Highlights significantly **upregulated**, **downregulated**, and **stable** genes.
#' - Users can **customize colors, axis labels, title, and save format**.
#'
#' @param deg_data Data frame containing DEG results. Output result from `merge_DEG_datasets()`
#'                 Must include columns: `logFC`, `FDR`, `contrast`, `gene`, `change`.
#' @param title Character. Title of the volcano plot (default: "Volcano Plot of DEG Analysis").
#' @param xlab Character. Label for the x-axis (default: "Group").
#' @param ylab Character. Label for the y-axis (default: "Log2 Fold Change").
#' @param colors Named list with colors for upregulated (`"Up"`), downregulated (`"Down"`), and stable (`"Stable"`) genes.
#'               Default: `list(Up = "#e6550d", Down = "#3182bd", Stable = "#636363")`.
#' @param point_size Numeric. Size of points in the plot (default: 2).
#' @param alpha_range Numeric vector. Range for transparency based on -log10(FDR) (default: c(0.3, 1)).
#' @param x_angle Numeric. Angle of x-axis labels (default: 45).
#' @param legend_position Character. Position of legend ("right", "top", "bottom", "left", "none"). Default: "right".
#' @param width Integer. The width of ggplot_draw (default: 10).
#' @param height Integer. The height of ggplot_draw (default: 6).
#' @param output_path Character. File path for saving the plot (default: "./DEG/").
#'
#' @importFrom ggplot2 ggplot geom_point aes geom_hline scale_color_manual scale_alpha labs theme_bw theme element_text ggsave
#'
#' @return A `ggplot` object representing the volcano plot.
#' \describe{
#'   \item{Visual Output}{A scatter plot with logFC on the y-axis and contrast groups on the x-axis.}
#'   \item{Customization}{Users can further modify the plot using `ggplot2` functions.}
#'   \item{File Output}{The plot is also saved to the specified `output_path` in the chosen format.}
#' }
#'
#' @examples
#' \dontrun{
#' # Generate default volcano plot
#' p <- plot_DEG_volcano(deg_data = merge.data)
#' p # Display the plot
#'
#' # Modify the returned plot
#' p + theme_minimal()
#' }
#'
#' @export
plot_DEG_volcano <- function(deg_data,
                             title = "Volcano Plot of DEG Analysis",
                             xlab = "Group",
                             ylab = "Log2 Fold Change",
                             colors = list(Up = "#e6550d", Down = "#3182bd", Stable = "#636363"),
                             point_size = 2,
                             alpha_range = c(0.3, 1),
                             x_angle = 45,
                             legend_position = "right",
                             width = 10,
                             height = 6,
                             output_path = "./DEG/") {
  # Ensure directory exists
  dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)

  # Create volcano plot
  p <- ggplot(deg_data) +
    geom_point(aes(x = contrast, y = logFC, color = change, size = -log10(FDR), alpha = -log10(FDR)),
      position = "jitter"
    ) +
    geom_hline(yintercept = c(-1, 1), linetype = "dotdash", color = "grey30") +
    scale_color_manual(values = colors) +
    scale_alpha(range = alpha_range) +
    labs(title = title, x = xlab, y = ylab) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x = element_text(angle = x_angle, vjust = 1, hjust = 1),
      legend.position = legend_position
    )

  # Save plot in the specified format
  file_ext <- paste0(output_path, "volcano_plot.pdf")
  ggsave(file_ext, plot = p, width = width, height = height)
  message("Volcano plot saved to: ", file_ext)

  return(p)
}
