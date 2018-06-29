#' @title Plot the factor/loading pairs from flash results.
#'
#' @return List containing:
#'
#'   \item{\code{plot_f}}{A ggplot object for the factors.}
#'
#'   \item{\code{plot_l}}{A ggplot object for the loadings.}
#'
#' @param data The flash data object.
#'
#' @param f The flash fit object.
#'
#' @param k The indices of the factor/loading pairs to be plotted.
#'
#' @param loading_label The label to use for the loading plots.
#'
#' @param factor_label The label to use for the factor plots.
#'
#' @details Plots each factor and loading as a barplot.
#'
#' @export
#'
flash_plot_factors =
  function(data, f, k = 1, loading_label = FALSE, factor_label = FALSE) {
    Y = get_Yorig(data)
    sample_name = rownames(Y)
    variable_name = colnames(Y)

    # Plot the expectation of PVE.
    K = flash_get_k(f)
    pve = flash_get_pve(f)

    # Plot the factors.
    if (factor_label == TRUE) {
        plot_f = plot_one_factor(flash_get_f(f), pve[k], k,
          f_labels = colnames(Y), y_lab = "factor values")
    } else {
        plot_f = plot_one_factor(flash_get_f(f), pve[k], k,
          f_labels = NA, y_lab = "factor values")
    }

    # Plot the loadings.
    if (loading_label == TRUE) {
        plot_l = plot_one_factor(flash_get_l(f), pve[k], k,
            f_labels = row.names(Y), y_lab = "loading values")
    } else {
        plot_l = plot_one_factor(flash_get_l(f), pve[k], k,
            f_labels = NA, y_lab = "loading values")
    }
    return(list(plot_f = plot_f, plot_l = plot_l))
}

#' @title Factor plot.
#'
#' @return list of factor, loading and variance of noise matrix
#'
#'    \item{\code{plot_f}} {A ggplot object for the factors.}
#'
#' @param f Factor to plot.
#'
#' @param pve PVE for this factor.
#'
#' @param k The order of the factor.
#'
#' @param f_labels The labels for the factor.
#'
#' @param y_lab The name of the Y axis.
#'
#' @details Plots the factors in a barplot.
#'
#' @importFrom ggplot2 ggplot aes_string scale_fill_manual labs geom_bar
#' @importFrom ggplot2 geom_text theme_minimal ylim
#'
#' @export
#'
plot_one_factor = function(f, pve, k, f_labels = NA, y_lab = "factor values") {
    P = length(f)
    if (any(is.na(f_labels))) {
        f_dat <- data.frame(variable = 1:P, Factor = f,
                            sign.f = factor(sign(f)),
                            hjust = factor(sign(f)))

        plot_f = ggplot(f_dat,aes_string(x = "variable",y = "Factor",
                                         fill = "sign.f"),
                        environment = environment()) +
            geom_bar(stat = "identity", width = 0.5) +
            scale_fill_manual(values = c("blue", "red")) +
            theme_minimal() +
            labs(title = paste("factor",k, "with PVE=", round(pve, 3)),
                 y = y_lab)
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
                                   label = "variablenames",fill = "sign.f"),
                        environment = environment()) +
          geom_bar(stat = "identity", width = 0.5) +
          geom_text(size = 2.75, angle = 90,
                    hjust = as.character(f_dat$hjust),
                 nudge_y = sign(f_dat$Factor)*0.1*mean(abs(f_dat$Factor))) +
          scale_fill_manual(values = c("blue", "red")) +
          ylim(lower_f, upper_f) +
          theme_minimal() +
          labs(title = paste("factor", k, "with PVE=", round(pve, 3)),
               y = y_lab)
    }
    return(plot_f)
}

#' @title Scree plot.
#'
#' @description Create a scree plot giving proportion of variance
#' explained by each factor.
#'
#' @param f The flash fit object.
#'
#' @param main The main caption to use for the plot.
#'
#' @return A ggplot plot object.
#'
#' @importFrom ggplot2 ggplot geom_point geom_line labs aes_string
#'
#' @export
#'
flash_plot_pve = function(f, main = "Scree plot of PVE for each factor") {
  pve = flash_get_pve(f)
  pve_dat = data.frame(factor_index = seq(1, flash_get_k(f)), PVE = pve)
  p <- ggplot(pve_dat, aes_string("factor_index","PVE"),
              environment = environment()) +
    geom_point(size = 4) + geom_line(linetype = "dotdash") +
    labs(title = main)
  return(p)
}
