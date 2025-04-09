#' PCA Plot Using `tinyarray::draw_pca()`
#'
#' This function performs PCA analysis and visualizes the results using `tinyarray::draw_pca()`.
#'
#' @importFrom tinyarray draw_pca
#'
#' @param pca_data List. Output from `prepare_PCA_data()`, containing merged expression data.
#' @param batch_correction Logical. Whether to use batch-corrected data (default: `TRUE`).
#' @param colors Character vector. Colors for different groups.
#'
#' @examples
#' \dontrun{
#' # draw PCA
#' plot_PCA_tinyarray(pca_data = pca_data, batch_correction = T)
#' }
#'
#' @export
plot_PCA_tinyarray <- function(pca_data,
                               batch_correction = TRUE,
                               colors = NULL) {
  # **Ensure input data is correct**
  if (!is.list(pca_data) || !all(c("merged_data", "group_all", "batch", "corrected_data") %in% names(pca_data))) {
    stop("Error: Input must be the output from `prepare_PCA_data()`.")
  }

  # **Select PCA data**
  data_for_pca <- if (batch_correction && !is.null(pca_data$corrected_data)) {
    pca_data$corrected_data
  } else {
    pca_data$merged_data
  }

  # **Ensure data is a numeric matrix**
  if (!is.matrix(data_for_pca)) {
    data_for_pca <- as.matrix(data_for_pca)
  }

  # **Check for NA values**
  if (any(is.na(data_for_pca))) {
    warning("Warning: PCA data contains NA values, which may cause issues in PCA plotting.")
  }

  # **Convert sample groups to factor**
  group_all_factor <- as.factor(pca_data$group_all)

  # **Check and adjust colors**
  num_groups <- length(levels(group_all_factor))
  if (is.null(colors)) {
    colors <- c(
      "#2874C5", "#F87669", "#E6B707", "#868686", "#92C5DE", "#F4A582",
      "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494",
      "#B3B3B3", "#2874C5", "#F87669"
    )[1:num_groups] # Select first `num_groups` colors
  }

  if (length(colors) != num_groups) {
    warning("Warning: The number of colors does not match the number of groups. Adjusting automatically.")
    colors <- rep(colors, length.out = num_groups)
  }

  # plot PCA
  draw_pca(
    exp = data_for_pca,
    group_list = group_all_factor,
    color = colors
  )
}
