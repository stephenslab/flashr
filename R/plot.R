#' @title Plot flash object
#'
#' @description \code{plot} method for class \code{'flash'}.
#'
#' @param x Flash object to plot.
#'
#' @param plot_factors Whether to plot the factors indexed by
#'   \code{kset}.
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
#'   \code{kset}.
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
#' @param plots_per_screen The maximum number of factor/loading plots to
#'   be shown per screen.
#'
#' @param facet_wrap_ncol The number of columns to use when displaying
#'   multiple factor/loading plots on a single screen.
#'
#' @param ask Should the user be prompted before displaying each
#'   successive plot?
#'
#' @importFrom grDevices dev.interactive devAskNewPage
#' @importFrom graphics plot
#'
#' @export
#'
plot.flash = function(x,
                      plot_factors = FALSE,
                      factor_kset = 1:x$nfactors,
                      factor_colors = NULL,
                      factor_legend_size = 5,
                      plot_loadings = FALSE,
                      loading_kset = 1:x$nfactors,
                      loading_colors = NULL,
                      loading_legend_size = 5,
                      plots_per_screen = 4,
                      facet_wrap_ncol = 2,
                      ask = (plot_factors || plot_loadings)
                              && dev.interactive()) {
  if (ask) {
    old_ask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(old_ask))
  }

  plot(plot_pve(x))

  if (plot_factors) {
    idx = 1
    while (idx <= length(factor_kset)) {
      next_kset = factor_kset[idx:min(idx + plots_per_screen - 1,
                                      length(factor_kset))]
      idx = idx + plots_per_screen
      plot(plot_kset(x, next_kset, factor_colors, factor_legend_size,
                     facet_wrap_ncol, factors = TRUE))
    }
  }

  if (plot_loadings) {
    idx = 1
    while (idx <= length(loading_kset)) {
      next_kset = loading_kset[idx:min(idx + plots_per_screen - 1,
                                       length(loading_kset))]
      idx = idx + plots_per_screen
      plot(plot_kset(x, next_kset, loading_colors, loading_legend_size,
                     facet_wrap_ncol, factors = FALSE))
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
#' @importFrom ggplot2 ggplot aes geom_bar scale_fill_manual
#' @importFrom ggplot2 scale_x_discrete ylim theme_grey theme labs
#' @importFrom ggplot2 element_text element_blank facet_wrap guides
#' @importFrom reshape2 melt
#'
plot_kset = function(f,
                     kset,
                     bar_colors = NULL,
                     legend_size = 5,
                     facet_wrap_ncol = 2,
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
    plot_object = ggplot(data, aes(x = variable, y = value)) +
      geom_bar(stat = "identity", width = 0.6) +
      scale_x_discrete(labels = NULL) +
      ylim(min_val, max_val) +
      theme_grey() +
      theme(legend.text = element_text(size = legend_size)) +
      labs(y = "", x = "") +
      facet_wrap(~k, ncol = facet_wrap_ncol)
  } else {
    plot_object = ggplot(data, aes(x = variable, y = value, fill = factor(variable))) +
      geom_bar(stat = "identity", width = 0.6) +
      scale_fill_manual(values = bar_colors) +
      scale_x_discrete(labels = NULL) +
      ylim(min_val, max_val) +
      theme_grey() +
      theme(legend.position="right",
            legend.text = element_text(size = legend_size),
            legend.title = element_blank()) +
      labs(y = "", x = "") +
      facet_wrap(~k, ncol = facet_wrap_ncol) +
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
#' @importFrom ggplot2 ggplot aes geom_point geom_line ylim labs
#'
plot_pve = function(f) {
  pve_dat = data.frame(factor_index = seq(1, length(f$pve)), PVE = f$pve)

  plot_object = ggplot(pve_dat,
         aes(x = factor_index, y = PVE),
         environment = environment()) +
    geom_point(size = 2) +
    geom_line(linetype = "dotdash") +
    ylim(0, NA) +
    labs(title = "Proportion of variance explained per factor/loading",
         x = "factor/loading index",
         y = "")

  return(plot_object)
}
