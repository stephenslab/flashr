#' @title Plot flash object
#'
#' @description \code{plot} method for class \code{'flash'}.
#'
#' @param x Flash object to plot.
#'
<<<<<<< HEAD
#' @param plot_pve Whether to include a scree plot of the proportion of
#'   variance explained by each factor/loading pair.
=======
#' @param plot_scree Whether to include a scree plot of the proportion
#'  of variance explained per factor/loading pair.
>>>>>>> master
#'
#' @param plot_factors Whether to plot the factors indexed by
#'   \code{factor_kset}.
#'
#' @param factor_kset If \code{plot_factors} is \code{TRUE}, then
#'   \code{factor_kset} specifies which factors will be plotted.
#'   Defaults to all factors.
#'
#' @param factor_colors If \code{plot_factors} is \code{TRUE}, then
#'   \code{factor_colors} specifies the colors to be used for the factor
#'   variables. If colors are used, then a legend will be shown
#'   alongside the factor plots.
#'
#' @param factor_legend_size If \code{factor_colors} is not \code{NULL},
#'   then \code{factor_legend_size} specifies the size of the legend
#'   show alongside the factor plots.
#'
#' @param plot_loadings Whether to plot the loadings indexed by
#'   \code{loading_kset}.
#'
#' @param loading_kset If \code{plot_loadings} is \code{TRUE}, then
#'   \code{loading_kset} specifies which loadings will be plotted.
#'   Defaults to all loadings.
#'
#' @param loading_colors If \code{plot_loadings} is \code{TRUE}, then
#'   \code{loading_colors} specifies the colors to be used for the
#'   loading variables. If colors are used, then a legend will be shown
#'   alongside the loading plots.
#'
#' @param loading_legend_size If \code{loading_colors} is not \code{NULL},
#'   then \code{loading_legend_size} specifies the size of the legend
#'   show alongside the loading plots.
#'
#' @param plot_grid_nrow The number of rows to use when displaying
#'   multiple factor/loading plots on a single screen.
#'
#' @param plot_grid_ncol The number of columns to use when displaying
#'   multiple factor/loading plots on a single screen.
#'
#' @param ask Should the user be prompted before displaying each
#'   successive plot?
#'
#' @param ... Additional arguments (not used by this method).
#'
#' @importFrom grDevices dev.interactive devAskNewPage
#' @importFrom graphics plot
#'
#' @export
#'
plot.flash = function(x,
<<<<<<< HEAD
                      plot_pve = TRUE,
=======
                      plot_scree = TRUE,
>>>>>>> master
                      plot_factors = FALSE,
                      factor_kset = 1:x$nfactors,
                      factor_colors = NULL,
                      factor_legend_size = 5,
                      plot_loadings = FALSE,
                      loading_kset = 1:x$nfactors,
                      loading_colors = NULL,
                      loading_legend_size = 5,
                      plot_grid_nrow = 2,
                      plot_grid_ncol = 2,
                      ask = (plot_factors || plot_loadings)
                              && dev.interactive(),
                      ...) {
  if (!plot_pve && !plot_factors && !plot_loadings) {
    stop(paste("Nothing to plot. Set plot_pve, plot_factors, or",
               "plot_loadings to TRUE."))
  }

  if (ask) {
    old_ask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(old_ask))
  }

<<<<<<< HEAD
  if (plot_pve) {
=======
  if(!plot_scree && !plot_factors && !plot_loadings) {
    stop(paste("Nothing to do if plot_scree, plot_factors, and",
               "plot_loadings are all FALSE."))
  }

  if (plot_scree && x$nfactors < 2) {
    warning("Not enough factors to create a scree plot.")
  } else if (plot_scree) {
>>>>>>> master
    plot(plot_pve(x))
  }

  plots_per_screen = plot_grid_nrow * plot_grid_ncol

  if (plot_factors) {
    idx = 1
    while (idx <= length(factor_kset)) {
      next_kset = factor_kset[idx:min(idx + plots_per_screen - 1,
                                      length(factor_kset))]
      idx = idx + plots_per_screen
      plot(plot_kset(x, next_kset, factor_colors, factor_legend_size,
                     plot_grid_ncol, factors = TRUE))
    }
  }

  if (plot_loadings) {
    idx = 1
    while (idx <= length(loading_kset)) {
      next_kset = loading_kset[idx:min(idx + plots_per_screen - 1,
                                       length(loading_kset))]
      idx = idx + plots_per_screen
      plot(plot_kset(x, next_kset, loading_colors, loading_legend_size,
                     plot_grid_ncol, factors = FALSE))
    }
  }
}


# @title Plot factors or loadings
#
# @return A ggplot object for the factors/loadings.
#
# @inheritParams plot.flash
#
# @param factors If TRUE, factors will be plotted. If FALSE, loadings
#   will be plotted.
#
#' @importFrom ggplot2 ggplot aes_string geom_bar scale_fill_manual
#' @importFrom ggplot2 scale_x_discrete ylim theme_grey theme labs
#' @importFrom ggplot2 element_text element_blank facet_wrap guides
#' @importFrom ggplot2 guide_legend
#' @importFrom reshape2 melt
#'
plot_kset = function(f,
                     kset,
                     bar_colors = NULL,
                     legend_size = 5,
                     plot_grid_ncol = 2,
                     factors = TRUE) {
  if (factors) {
    vals = f$ldf$f
  } else {
    vals = f$ldf$l
  }
  vals = vals[, kset, drop = FALSE]
  min_val = min(0, min(vals))
  max_val = max(0, max(vals))

  pve = pmax(round(f$pve, 3), 0.001)[kset]
  data = melt(vals)
  colnames(data) = c("variable", "k", "value")

  var_labels = rownames(vals)
  if (is.null(var_labels)) {
    var_labels = as.character(1:nrow(vals))
  }
  data$variable = factor(data$variable,
                         levels = var_labels,
                         labels = var_labels)

  if (factors) {
    title = "Factor"
  } else {
    title = "Loading"
  }
  plot_titles = paste0(title, " ", kset, "; pve: ", pve)
  data$k = factor(data$k, levels = 1:length(kset), labels = plot_titles)

  if (is.null(bar_colors)) {
    plot_object = ggplot(data, aes_string(x = "variable", y = "value")) +
      geom_bar(stat = "identity", width = 0.6) +
      scale_x_discrete(labels = NULL) +
      ylim(min_val, max_val) +
      theme_grey() +
      theme(legend.text = element_text(size = legend_size)) +
      labs(y = "", x = "") +
      facet_wrap(~k, ncol = plot_grid_ncol)
  } else {
    plot_object = ggplot(data, aes_string(x = "variable", y = "value",
                                          fill = "variable")) +
      geom_bar(stat = "identity", width = 0.6) +
      scale_fill_manual(values = bar_colors) +
      scale_x_discrete(labels = NULL) +
      ylim(min_val, max_val) +
      theme_grey() +
      theme(legend.position="right",
            legend.text = element_text(size = legend_size),
            legend.title = element_blank()) +
      labs(y = "", x = "") +
      facet_wrap(~k, ncol = plot_grid_ncol) +
      guides(fill = guide_legend(ncol = 1,
                                 keyheight = legend_size / 6,
                                 keywidth = legend_size / 15))
  }

  return(plot_object)
}


# @title Plot PVE
#
# @description Create a scree plot giving the proportion of variance
#   explained by each factor.
#
# @param f A flash object.
#
# @return A \pkg{ggplot} plot object.
#
#' @importFrom ggplot2 ggplot aes_string geom_point geom_line ylim labs
#'
plot_pve = function(f) {
  pve_dat = data.frame(factor_index = seq(1, length(f$pve)), PVE = f$pve)

  plot_object = ggplot(pve_dat,
         aes_string(x = "factor_index", y = "PVE"),
         environment = environment()) +
    geom_point(size = 2) +
    geom_line(linetype = "dotdash") +
    ylim(0, NA) +
    labs(title = "Proportion of variance explained per factor/loading",
         x = "factor/loading index",
         y = "")

  return(plot_object)
}
