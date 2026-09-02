#' Plot Finlay–Wilkinson regression lines highlighting Top-N and checks
#'
#' Draws genotype-specific regression profiles across the environment gradient,
#' highlighting the Top-N genotypes (by mean predicted value across the gradient)
#' and the check genotypes.
#'
#' Input compatibility:
#' - Uses `post_prob$across$reg` created from the Legendre basis (`gradient_basis`)
#'   in `prob_sup_bglr`, where predictions across the adaptive gradient are computed as
#'   #$# \hat{y}_{g}(l) = \sum_{d=0}^{D} \beta_{g,d}\,\mathrm{deg}_d(l) #$#,
#'   with `deg0 = 1` and `deg1..degD` being orthonormalized Legendre basis columns.
#'
#' @param post_prob An object produced by `prob_sup_bglr()`. It must contain
#'   `post_prob$across$reg`, a data frame with columns:
#'   - `gradient`: environmental index values,
#'   - `value`: predicted response along the gradient,
#'   - `ID`: genotype ID,
#'   - `check`: logical flag indicating checks.
#' @param top_n Integer. Number of top genotypes to highlight (default = 10).
#'
#' @details
#' - The Top-N genotypes are defined by the highest mean of `value` across the gradient.
#' - Lines are colored and thickened for `TopN` and `Check` categories; others are thin and grey.
#' - Labels are added at the last gradient point only for highlighted lines.
#' - The legend labels adapt dynamically to `top_n` (e.g., "Top10", "Top30").
#'
#' @return A ggplot object.
#'
#' @examples
#' \donttest{
#' # Assuming you already ran prob_sup_bglr():
#' p1 <- plot_fw(post_prob = probs, top_n = 10)
#' print(p1)
#' }
#'
#' @seealso prob_sup_bglr
#' @export
#' @importFrom dplyr group_by summarize arrange slice_head pull mutate case_when filter ungroup
#' @importFrom ggplot2 ggplot aes geom_line scale_color_manual scale_linewidth_manual guides theme_bw theme element_text
#' @importFrom ggrepel geom_text_repel
plot_fw <- function(post_prob,
                    top_n = 10){
  
  library(ggplot2)
  library(ggrepel)
  
  df <- post_prob$across$reg 
   # dplyr::mutate(value = value + gradient) |> 
  top_label <- paste0("Top", top_n)
  
  top_ids <- df |> 
    dplyr::group_by(ID) |> 
    dplyr::summarize(metric = mean(value, na.rm = TRUE), .groups = "drop") |> 
    dplyr::arrange(desc(metric)) |> 
    dplyr::slice_head(n = top_n) |> 
    dplyr::pull(ID)
  
  df_plot <- df |> 
    dplyr::mutate(category = dplyr::case_when(
      check ~ "Check",
      ID %in% top_ids ~ top_label,
      TRUE ~ "Other"
    ))
  
  df_lab <- df_plot |> 
    dplyr::group_by(ID, category) |> 
    dplyr::filter(gradient == max(gradient)) |> 
    dplyr::ungroup() |> 
    dplyr::filter(category != "Other")
  
  # Dynamic scales for color and linewidth
  color_map <- c(Check = "#D62728", Other = "grey70")
  color_map[top_label] <- "#1F77B4"
  lw_map <- c(Check = 0.9, Other = 0.4)
  lw_map[top_label] <- 0.9
  
  ggplot(df_plot, aes(x = gradient, y = value, group = ID,
                      color = category, linewidth = category)) +
    geom_line() +
    geom_text_repel(data = df_lab, aes(label = ID),
                    min.segment.length = 0, size = 3,
                    box.padding = 0.25, point.padding = 0.2, show.legend = FALSE) +
    scale_color_manual(values = color_map) +
    scale_linewidth_manual(values = lw_map) +
    guides(color = guide_legend(title = "Category"),
           linewidth = "none") +
    theme_bw() +
    theme(legend.position = "right",
          axis.title = element_text(size = 11),
          axis.text = element_text(size = 9))
}



#' Lollipop plot of across-environment top-tail probability (intercept-based)
#'
#' Plots the probability of being in the top tail (as defined in `prob_sup_bglr`, e.g., top 20%)
#' for genotype intercepts, highlighting Top-N genotypes and checks.
#'
#' Legendre compatibility:
#' - Intercepts are genotype-level #$# \alpha_g #$# means; slope terms (Legendre coefficients)
#'   do not enter here. Top-N is ranked by the mean intercept.
#'
#' @param post_prob Output from `prob_sup_bglr()` containing `across$bayesian_stats`
#'   with columns `ID`, `check`, `intercept`, and `prob_top`.
#' @param top_n Integer. Number of genotypes to highlight by highest average intercept (default = 30).
#'
#' @details
#' - The Top-N set is formed by ranking genotypes on the mean `intercept`.
#' - The y-axis shows `prob_top`: probability of being in the global top tail defined by `int`
#'   (from `prob_sup_bglr(control$int)`).
#'
#' @return A ggplot object.
#'
#' @examples
#' \donttest{
#' p2 <- plot_prob_intercept(post_prob = probs, top_n = 30)
#' print(p2)
#' }
#'
#' @seealso prob_sup_bglr
#' @export
#' @importFrom dplyr group_by summarize arrange slice_head pull mutate case_when filter
#' @importFrom ggplot2 ggplot aes geom_segment geom_point labs ggtitle theme_bw theme element_text reorder
plot_prob_intercept <- function(post_prob,
                                top_n = 30){
  
  library(ggplot2)
  
  df <- post_prob$across$bayesian_stats
  top_label <- paste0("Top", top_n)
  
  top_ids <- df |> 
    dplyr::group_by(ID) |> 
    dplyr::summarize(metric = mean(intercept, na.rm = TRUE), .groups = "drop") |> 
    dplyr::arrange(desc(metric)) |> 
    dplyr::slice_head(n = top_n) |> 
    dplyr::pull(ID)
  
  df_plot <- df |> 
    dplyr::mutate(category = dplyr::case_when(
      check ~ "Check",
      ID %in% top_ids ~ top_label,
      TRUE ~ "Other"
    )) |> 
    dplyr::filter(category %in% c(top_label, "Check"))
  
  ggplot(df_plot, aes(x = reorder(.data$ID, -.data$prob_top), y = .data$prob_top)) +
    geom_segment(aes(xend = reorder(.data$ID, -.data$prob_top), 
                     y = 0, yend = .data$prob_top), 
                 linewidth = 1, color = "gray60") +
    geom_point(aes(fill = category), size = 3, shape = 21, color = "black") +
    labs(x = "Genotype", y = "Probability") +
    ggtitle("Probability of Being in the Top Tail (Intercept)")+
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}



#' Lollipop plot of slope sensitivity probability (Pr(beta > 1))
#'
#' Plots the probability of degree-1 slope being greater than 1 for Top-N genotypes and checks,
#' highlighting sensitivity to the environment index relative to 1.
#'
#' Legendre note:
#' - Under orthonormal Legendre basis, the degree-1 term multiplies #$# \tilde P_1(x) #$#,
#'   so a threshold at 0 is more natural for directionality; `> 1` is kept for backward compatibility.
#'
#' @param post_prob Output from `prob_sup_bglr()` containing `across$bayesian_stats`
#'   with columns `ID`, `check`, and `slope_gt_1`.
#' @param top_n Integer. Number of genotypes to highlight by highest average slope (default = 30).
#'
#' @details
#' - Top-N set is formed by ranking genotypes on the mean `slope` (degree-1 coefficient).
#'
#' @return A ggplot object.
#'
#' @examples
#' \donttest{
#' p3 <- plot_prob_slope(post_prob = probs, top_n = 30)
#' print(p3)
#' }
#'
#' @seealso prob_sup_bglr
#' @export
#' @importFrom dplyr group_by summarize arrange slice_head pull mutate case_when filter
#' @importFrom ggplot2 ggplot aes geom_segment geom_point labs ggtitle theme_bw theme element_text reorder
plot_prob_slope <- function(post_prob,
                            top_n = 30){
  
  library(ggplot2)
  
  df <- post_prob$across$bayesian_stats
  top_label <- paste0("Top", top_n)
  
  top_ids <- df |> 
    dplyr::group_by(ID) |> 
    dplyr::summarize(metric = mean(slope, na.rm = TRUE), .groups = "drop") |> 
    dplyr::arrange(desc(metric)) |> 
    dplyr::slice_head(n = top_n) |> 
    dplyr::pull(ID)
  
  df_plot <- df |> 
    dplyr::mutate(category = dplyr::case_when(
      check ~ "Check",
      ID %in% top_ids ~ top_label,
      TRUE ~ "Other"
    )) |> 
    dplyr::filter(category %in% c(top_label, "Check"))
  
  ggplot(df_plot, aes(x = reorder(.data$ID, -.data$slope_gt_1), y = .data$slope_gt_1)) +
    geom_segment(aes(xend = reorder(.data$ID, -.data$slope_gt_1), 
                     y = 0, yend = .data$slope_gt_1), 
                 linewidth = 1, color = "gray60") +
    geom_point(aes(fill = category), size = 3, shape = 21, color = "black") +
    labs(x = "Genotype", y = "Probability") +
    ggtitle("Probability of Slope > 1 (Degree 1)")+
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}



#' Scatter of FW parameters: slope vs intercept colored by R2
#'
#' Visualizes genotype-level FW parameters with slope on the x-axis and intercept on the y-axis,
#' coloring by an R2 statistic (e.g., goodness-of-fit of the FW regression per genotype).
#'
#' Legendre-compatible:
#' - The slope here corresponds to the degree-1 coefficient from the Legendre random regression
#'   (as reported in `bayesian_stats$slope`). Intercept is #$# \alpha_g #$#.
#'
#' @param post_prob Output from `prob_sup_bglr()` containing `across$bayesian_stats`
#'   with columns `slope`, `intercept`, and an `R2` column.
#'
#' @details
#' - This function assumes that `post_prob$across$bayesian_stats` contains an `R2` column.
#'   If `R2` was not computed earlier, you should add it before calling this function
#'   (e.g., based on observed vs fitted FW lines per genotype).
#'
#' @return A ggplot object.
#'
#' @examples
#' \donttest{
#' # If bayesian_stats has an R2 column:
#' p4 <- plot_param(post_prob = probs)
#' print(p4)
#' }
#'
#' @seealso prob_sup_bglr
#' @export
#' @importFrom ggplot2 ggplot aes geom_point theme_bw
#' @importFrom viridis scale_color_viridis
plot_param <- function(post_prob){
  
  library(ggplot2)
  
  df <- post_prob$across$bayesian_stats
  
  ggplot(df, aes(x = slope, y = intercept, color = R2))+
    geom_point()+
    viridis::scale_color_viridis()+
    theme_bw()
}



#' Heatmap of within-environment top-tail probabilities
#'
#' Draws a heatmap of environment-specific probabilities for Top-N genotypes,
#' where the probabilities represent within-environment tail probabilities computed in `prob_sup_bglr()`.
#'
#' Legendre-compatible:
#' - Probabilities are based on environment-specific posterior predictions from the Legendre
#'   random regression (`within$perfo`).
#'
#' @param post_prob Output from `prob_sup_bglr()` containing:
#'   - `across$bayesian_stats` with `ID` and `prob_top` for ranking,
#'   - `within$perfo` with columns `env`, `ID`, and `prob` for environment-specific probabilities.
#' @param top_n Integer. Number of genotypes to display (ranked by `prob_top`). Default = 30.
#'
#' @return A ggplot object.
#'
#' @examples
#' \donttest{
#' p5 <- plot_prob_within(post_prob = probs, top_n = 30)
#' print(p5)
#' }
#'
#' @seealso prob_sup_bglr
#' @export
#' @importFrom dplyr arrange slice_head pull filter
#' @importFrom ggplot2 ggplot aes geom_tile labs theme_bw theme element_text scale_fill_viridis_c
plot_prob_within <- function(post_prob,
                             top_n = 30){
  
  library(ggplot2)
  
  mean_perf <- post_prob$across$bayesian_stats
  
  top_ids <- mean_perf |> 
    dplyr::arrange(desc(prob_top)) |> 
    dplyr::slice_head(n = top_n) |> 
    dplyr::pull(ID)
  
  df <- post_prob$within$perfo |> 
    dplyr::filter(ID %in% top_ids)
  
  ggplot(df, aes(x = .data$env, y = reorder(.data$ID, .data$prob), fill = .data$prob)) +
    geom_tile(color = "white") +
    labs(x = "Environment", y = "Genotype", title = "Probability Within Environment") +
    scale_fill_viridis_c(direction = 1, na.value = 'white', limits = c(0,1), option = 'turbo') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}



#' Heatmap of pairwise probabilities of superior performance (Top-N)
#'
#' Visualizes pairwise probabilities among the Top-N genotypes using the
#' matrix of pairwise performance probabilities computed in `prob_sup_bglr()`.
#'
#' @param post_prob Output from `prob_sup_bglr()` containing:
#'   - `across$bayesian_stats` (for ranking by `prob_top`),
#'   - `across$pair_perfo` (square matrix of pairwise probabilities),
#'   - `control$increase` (flag indicating whether higher is better), typically set in `prob_sup_bglr`.
#' @param top_n Integer. Number of top genotypes to display. Default = 30.
#'
#' @details
#' - The function subsets the pairwise matrix to the Top-N genotypes and masks the
#'   lower triangle and diagonal for clearer visualization.
#' - The legend label is determined by the `increase` flag:
#'   `"Pr(X > Y)"` if higher is better, otherwise `"Pr(X < Y)"`.
#'
#' @return A ggplot object.
#'
#' @examples
#' \donttest{
#' p6 <- plot_pair_prob(post_prob = probs, top_n = 30)
#' print(p6)
#' }
#'
#' @seealso prob_sup_bglr
#' @export
#' @importFrom dplyr arrange slice_head pull
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_viridis_c labs theme_bw theme element_text element_blank coord_fixed
plot_pair_prob <- function(post_prob,
                           top_n = 30){
  
  mean_perf <- post_prob$across$bayesian_stats
  control <- attr(post_prob, "control")
  
  top_ids <- mean_perf |> 
    dplyr::arrange(desc(prob_top)) |> 
    dplyr::slice_head(n = top_n) |> 
    dplyr::pull(ID)
  
  pair_mat <- post_prob$across$pair_perfo
  pair_mat <- pair_mat[top_ids, top_ids, drop = FALSE]
  pair_mat[lower.tri(pair_mat, diag = TRUE)] <- NA
  
  melted <- reshape2::melt(pair_mat, na.rm = TRUE)
  colnames(melted) <- c("x", "y", "prob")
  
  ggplot2::ggplot(melted, ggplot2::aes(x = factor(.data$x, levels = top_ids), 
                                       y = factor(.data$y, levels = rev(top_ids)), 
                                       fill = .data$prob)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_viridis_c(limits = c(0, 1), option = "turbo") +
    ggplot2::labs(
      x = "Genotype", 
      y = "Genotype", 
      fill = if(control$increase) "Pr(X > Y)" else "Pr(X < Y)",
      title = "Pairwise Probability of Superior Performance"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
      axis.text.y = ggplot2::element_text(),
      panel.grid = ggplot2::element_blank(),
      legend.position = "right"
    ) +
    ggplot2::coord_fixed()
}
