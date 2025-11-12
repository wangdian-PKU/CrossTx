#' Plot Tumor Microenvironment (TME) Boxplots for Multiple Signatures
#'
#' This function computes immune infiltration scores using CIBERSORT and generates
#' boxplots for multiple immune signatures.
#'
#' @param eset_tpm Data frame. TPM-transformed RNA-seq expression matrix, output from 'plot_TME_barplot()'.
#' @param group_as Character vector. 'Identical' with the 'group_as' in 'plot_TME_barplot()'.
#' @param condition_list Named list. Sample count for each mouse model (e.g., 'list("DEN" = 6, "GSExxx" = 8)').
#' @param width Integer. Width of saved plots (default: 10).
#' @param height Integer. Height of saved plots (default: 6).
#' @param output_path Character. Directory path to save boxplots (default: "./TME_boxplots").
#'
#' @importFrom ggplot2 geom_boxplot labs aes theme_bw theme element_text element_line element_blank ggsave
#' @importFrom dplyr mutate select all_of
#' @importFrom tidyr separate
#' @importFrom IOBR deconvo_tme
#'
#' @return Saves multiple boxplots and returns a list of ggplot objects.
#'
#' @examples
#' \dontrun{
#' group_as <- c(rep("HBV_normal", 3), rep("HBV_tumor", 3), rep("DEN_normal", 4), rep("DEN_tumor", 4))
#' condition_list <- list("HBV" = 6, "DEN" = 8)
#' plot_TME_boxplot(eset_tpm = TME_barplot_result$eset_tpm, group_as = group_as, condition_list = condition_list)
#' }
#'
#' @export
plot_TME_boxplot <- function(eset_tpm,
                             group_as,
                             condition_list,
                             width = 10,
                             height = 6,
                             output_path = "./TME_boxplots") {
  # Step 1: Validate Inputs
  if (!is.data.frame(eset_tpm)) stop("Error: 'eset_tpm' must be a data frame.")
  if (length(group_as) != ncol(eset_tpm)) stop("Error: 'group_as' length must match number of samples (columns in 'eset_tpm').")

  # Step 2: Compute Immune Cell Fractions
  message("Computing immune infiltration scores using CIBERSORT...")
  method_model <- deconvo_tme(
    eset = eset_tpm,
    method = "cibersort",
    arrays = FALSE,
    perm = 200
  )

  # Step 3: Create Output Directory
  if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)

  # Step 4: Initialize plot list
  plot_list <- list()

  # Step 5: Generate boxplot for each signature
  for (sig_style in colnames(method_model)[-1]) {
    message("Generating boxplot for: ", sig_style)

    sig_score <- method_model %>%
      select(ID, all_of(sig_style)) %>%
      mutate(group = group_as) %>%
      separate(col = group, into = c("condition", "group"), sep = "_") %>%
      mutate(condition = factor(condition, levels = unique(condition)))

    colnames(sig_score)[2] <- "sig"

    p <- ggplot(data = sig_score) +
      geom_boxplot(aes(x = group, y = sig, fill = condition)) +
      labs(y = paste(sig_style, "score")) +
      theme_bw() +
      theme(
        axis.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 90),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.justification = c(1, 1)
      )

    # Step 6: Save each plot
    file_ext <- file.path(output_path, paste0(sig_style, ".pdf"))
    ggsave(file_ext, plot = p, height = height, width = width)

    message("Saved boxplot to: ", file_ext)
    plot_list[[sig_style]] <- p
  }

  return(plot_list)
}
