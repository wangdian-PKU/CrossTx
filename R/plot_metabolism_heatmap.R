#' Plot Metabolism Signature Heatmap
#'
#' This function generates a heatmap of metabolism-related signature scores.
#'
#' @param sig_meta Data frame. Output from `calculate_metabolism_score()`.
#' @param output_path Character. File path to save the heatmap (default: "./metabolism/").
#' @param clustering_method Character. Clustering method ("manhattan", "canberra"). Default: "manhattan".
#' @param width Numeric. Width of the output PDF (default: 10).
#' @param height Numeric. Height of the output PDF (default: 22).
#' @param row_name_width Numeric. Maximum width for row names in centimeters (default: 7).
#' @param right_padding Numeric. Right margin in centimeters to avoid label cutoff (default: 6).
#' @param row_height_pt Numeric. Height of each heatmap row in points (default: 5).
#'
#' @importFrom ComplexHeatmap Heatmap draw
#' @importFrom grid gpar unit
#'
#' @return Saves the heatmap and returns the heatmap object.
#'
#' @export
plot_metabolism_heatmap <- function(sig_meta,
                                    output_path = "./metabolism",
                                    clustering_method = "manhattan",
                                    width = 10,
                                    height = 22,
                                    row_name_width = 7,
                                    right_padding = 6,
                                    row_height_pt = 5) {
  if (!is.data.frame(sig_meta)) stop("Error: `sig_meta` must be a data frame.")
  if (!clustering_method %in% c("manhattan", "canberra")) {
    stop("Error: `clustering_method` must be 'manhattan' or 'canberra'.")
  }

  sig_meta_merge <- aggregate(sig_meta[, -1], by = list(group_as), FUN = mean)
  rownames(sig_meta_merge) <- sig_meta_merge$Group.1
  df <- sig_meta_merge[, -1] %>% t()

  heatmap_obj <- Heatmap(
    df,
    clustering_distance_columns = clustering_method,
    clustering_distance_rows = clustering_method,
    row_names_gp = gpar(fontsize = 8),
    row_names_max_width = unit(row_name_width, "cm"),
    height = unit(nrow(df) * row_height_pt, "pt")
  )

  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
  }

  output_file <- file.path(output_path, "heatmap.pdf")

  grDevices::pdf(output_file, width = width, height = height, useDingbats = FALSE)
  draw(heatmap_obj,
    heatmap_legend_side = "right",
    padding = unit(c(2, 2, 2, right_padding), "cm")
  )
  grDevices::dev.off()

  message("Heatmap saved to: ", output_file)

  return(heatmap_obj)
}
