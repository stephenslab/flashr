#' @title Plot the factor/loading pairs from flash results.
#'
#' @description Plots each factor and loading as a barplot.
#'
#' @return List containing:
#'
#'   \item{\code{plot_f}}{A ggplot object for the factors.}
#'
#'   \item{\code{plot_l}}{A ggplot object for the loadings.}
#'
#' @param data An n by p matrix or a flash data object created using
#'   \code{flash_set_data}.
#'
#' @param f A flash fit object.
#'
#' @param kset The indices of the factor/loading pairs to be plotted.
#'
#' @param loading_label If \code{TRUE}, then the row names of the data
#'   matrix will be used to label the loading plots.
#'
#' @param factor_label If \code{TRUE}, then the column names of the data
#'   matrix will be used to label the factor plots.
#'
#' @export
#'
flash_plot_factors = function(data,
                              f,
                              kset = NULL,
                              loading_label = FALSE,
                              factor_label = FALSE) {
  data = handle_data(data, output = "matrix")
  # f is handled by flash_get_pve
  # think about how to handle possible indexing problems caused by
  # dropping zero factors
  kset = handle_kset(kset, f)

  # Plot the expectation of PVE.
  pve = flash_get_pve(f, drop_zero_factors = FALSE)

  plot_f = list()
  plot_l = list()

  for (k in kset) {
    # Plot the factors.
    if (factor_label == TRUE) {
      f_labels = colnames(data)
    } else {
      f_labels = NA
    }
    plot_f[[k]] = plot_one_factor(flash_get_f(f)[, k],
                                  pve[k],
                                  k,
                                  f_labels = f_labels,
                                  y_lab = "factor values")

    # Plot the loadings.
    if (loading_label == TRUE) {
      f_labels = row.names(data)
    } else {
      f_labels = NA
    }
    plot_l[[k]] = plot_one_factor(flash_get_l(f)[, k],
                                  pve[k],
                                  k,
                                  f_labels = f_labels,
                                  y_lab = "loading values")
  }
  return(list(plot_f = plot_f, plot_l = plot_l))
}


# @title Factor plot.
#
# @return A ggplot object for the factors.
#
# @param f Factor to plot.
#
# @param pve PVE for this factor.
#
# @param k The order of the factor.
#
# @param f_labels The labels for the factor.
#
# @param y_lab The name of the Y axis.
#
# @details Plots the factors in a barplot.
#
#' @importFrom ggplot2 ggplot aes_string scale_fill_manual labs geom_bar
#' @importFrom ggplot2 geom_text theme_minimal ylim
#'
plot_one_factor = function(f,
                           pve,
                           k,
                           f_labels = NA,
                           x_lab = "variable",
                           y_lab = "factor values",
                           show_sign.f_legend = FALSE) {
  P = length(f)
  if (any(is.na(f_labels))) {
    f_dat <- data.frame(variable = 1:P, Factor = f,
                        sign.f = factor(sign(f)),
                        hjust = factor(sign(f)))

    plot_f = ggplot(f_dat, aes_string(x = "variable", y = "Factor",
                                      fill = "sign.f"),
                    environment = environment()) +
      geom_bar(stat = "identity", width = 0.5,
               show.legend = show_sign.f_legend) +
      scale_fill_manual(values = c("blue", "red")) +
      theme_minimal() +
      labs(title = paste("Factor", k, "with PVE =", round(pve, 3)),
           x = x_lab, y = y_lab)
  } else {
    f_dat <- data.frame(variable = 1:P, Factor = f,
                        sign.f = factor(sign(f)),
                        variablenames = f_labels,
                        hjust = factor(sign(f)))

    # 120% limit.
    range_f = max(f) - min(f)
    upper_f = max(f, 0) + 0.15 * range_f
    lower_f = min(f, 0) - 0.15 * range_f

    plot_f = ggplot(f_dat,
                    aes_string(x = "variable", y = "Factor",
                               label = "variablenames", fill = "sign.f"),
                    environment = environment()) +
      geom_bar(stat = "identity", width = 0.5,
               show.legend = show_sign.f_legend) +
      geom_text(size = 2.75, angle = 90,
                hjust = as.character(f_dat$hjust),
                nudge_y = sign(f_dat$Factor)*0.1*mean(abs(f_dat$Factor))) +
      scale_fill_manual(values = c("blue", "red")) +
      ylim(lower_f, upper_f) +
      theme_minimal() +
      labs(title = paste("Factor", k, "with PVE =", round(pve, 3)),
           x = x_lab, y = y_lab)
  }

  return(plot_f)
}


#' @title Scree plot
#'
#' @description Create a scree plot giving the proportion of variance
#' explained by each factor.
#'
#' @param f A flash fit object.
#'
#' @param main The main caption to use for the plot.
#'
#' @param drop_zero_factors If \code{TRUE}, then any factor/loadings
#'   that are zero will be removed.
#'
#' @return A \pkg{ggplot} plot object.
#'
#' @importFrom ggplot2 ggplot geom_point geom_line labs aes_string
#'
#' @export
#'
flash_plot_pve = function(f,
                          main = "Scree plot of PVE for each factor",
                          drop_zero_factors = TRUE) {
  # handling is done by flash_get_pve

  pve = flash_get_pve(f, drop_zero_factors)
  pve_dat = data.frame(factor_index = seq(1, length(pve)), PVE = pve)
  p <- ggplot(pve_dat, aes_string("factor_index", "PVE"),
              environment = environment()) +
       geom_point(size = 4) + geom_line(linetype = "dotdash") +
       labs(title = main, x = "factor index")

  return(p)
}
