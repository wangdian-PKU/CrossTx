#' Plot TROM Score Boxplot and Tileplot
#'
#' This function generates **two subplots**:
#' 1. **Boxplot** of transcriptional similarity scores (`TROM score`).
#' 2. **Tileplot** showing the **mean similarity scores per condition**.
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot xlab ylab scale_color_manual scale_fill_manual theme_classic geom_tile scale_fill_gradient2 guides guide_colorbar ggsave
#' @importFrom dplyr filter mutate
#' @importFrom patchwork plot_layout
#'
#' @param merge_long Data frame. Output from `plot_TROMscore_heatmap()`, containing **TROM scores**.
#' @param output_path Character. File path to save the figure (default: `"./similarity/TROM_boxplot"`).
#' @param combine_plots Logical. Whether to **combine tileplot & boxplot** into one figure (default: `TRUE`).
#'
#' @return Saves the plot(s) and returns the ggplot object.
#'
#' @examples
#' \dontrun{
#' plot_TROMscore_boxplot(
#'   merge_long = heatmap_results$merge_long,
#'   output_path = "./similarity/TROM_boxplot",
#'   combine_plots = TRUE
#' )
#' }
#'
#' @export
plot_TROMscore_boxplot <- function(merge_long,
                                   output_path = "./similarity/TROM_boxplot",
                                   combine_plots = TRUE) {
  # **Check if `merge_long` is valid**
  if (!is.data.frame(merge_long)) {
    stop("Error: `merge_long` must be a data frame.")
  }
  if (!all(c("condition", "score") %in% colnames(merge_long))) {
    stop("Error: `merge_long` must contain `condition` and `score` columns.")
  }
  if (nrow(merge_long) == 0) {
    warning("Warning: `merge_long` is empty. No data to plot.")
  }

  # **Step 1: Generate Boxplot**
  p_boxplot <- ggplot(merge_long) +
    geom_boxplot(aes(x, y, fill = group, color = group), alpha = 0.5, outlier.shape = NA) +
    facet_grid(. ~ condition, scales = "free_x", space = "free_x") +
    scale_color_manual(name = "", values = c("#5aadd0", "#af7aa1")) +
    scale_fill_manual(name = "", values = c("#5aadd0", "#af7aa1")) +
    xlab("") +
    ylab("Similarity score (TROM)") +
    geom_vline(xintercept = 2.6, alpha = 0.5, color = "grey") +
    scale_y_continuous(expand = c(0, 0)) +
    theme_classic() +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "right",
      legend.key.size = unit(15, "mm")
    )

  # tileplot_data preparation
  tile_data <- aggregate(merge_long$y, by = list(merge_long$condition), FUN = mean)

  colnames(tile_data) <- c("y", "x")

  tile_data <- data.frame(x = tile_data$x, y = tile_data$y)

  tile_data$y <- factor(tile_data$y, levels = unique(tile_data$y))

  # **Step 2: Generate Tileplot**
  p_tileplot <- ggplot(tile_data, aes(x = c(1:nrow(tile_data)), y = 1, fill = x)) +
    geom_tile(color = "white", linewidth = 1) +
    scale_x_discrete("", expand = c(0, 0)) +
    scale_y_discrete("", expand = c(0, 0)) +
    scale_fill_gradient2(low = "#4178a7", mid = "#ffffff", high = "#e2201c", midpoint = 1) +
    guides(fill = guide_colorbar("Score", title.position = "top")) +
    theme(
      panel.background = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y.left = element_text(size = 12),
      axis.ticks = element_blank(),
      legend.direction = "vertical",
      legend.position = "right"
    ) +
    geom_text(aes(label = unique(y), color = I(c("white", rep("black", (nrow(tile_data) - 1))))), angle = 90, size = 5)



  # **Step 3: Combine Plots (Optional)**
  final_plot <- if (combine_plots) p_tileplot / p_boxplot + plot_layout(heights = c(1.2, 5)) else list(p_tileplot, p_boxplot)

  # **Step 4: Save Plot**
  file_ext <- paste0(output_path, ".pdf")
  ggsave(file_ext, plot = final_plot, height = ifelse(combine_plots, 7, 6.83), width = 9.51, dpi = 300)

  message("Plot saved to: ", file_ext)
  return(final_plot)
}
