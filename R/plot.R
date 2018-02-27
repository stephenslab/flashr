#' @title Plot the factors from flash results.
#' 
#' @return List containing:
#' 
#'  \itemize{
#'   \item{\code{plot_f}} {A ggplot object for the factors.}
#'   \item{\code{plot_l}} {A ggplot object for the loadings}
#'  }
#' 
#' @param data The flash data object.
#' 
#' @param f A flash fit object.
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
    
    # plot the expectation of PVE
    K = flash_get_k(f)
    pve = flash_get_pve(f)
    
    # plot the factors
    if (factor_label == TRUE) {
        plot_f = plot_one_factor(flash_get_f(f, k), pve[k], k, f_labels = colnames(Y), y_lab = "factor values")
    } else {
        plot_f = plot_one_factor(flash_get_f(f, k), pve[k], k, f_labels = NA, y_lab = "factor values")
    }
    
    # plot the loadings
    if (loading_label == TRUE) {
        plot_l = plot_one_factor(flash_get_l(f, k), pve[k], k, f_labels = row.names(Y), y_lab = "loading values")
    } else {
        plot_l = plot_one_factor(flash_get_l(f, k), pve[k], k, f_labels = NA, y_lab = "loading values")
    }
    return(list(plot_f = plot_f, plot_l = plot_l))
}

#' @title Factor plot.
#' 
#' @return list of factor, loading and variance of noise matrix
#'  \itemize{
#'   \item{\code{plot_f}} {is a ggplot object for the factors}
#'  }
#' 
#' @param f Factor to plot.
#' 
#' @param pve PVE for this factor.
#' 
#' @param P The length of this factor.
#' 
#' @param k The order of the factor.
#' 
#' @param f_labels The labels for the factor.
#' 
#' @param y_lab The name of the Y axis.
#' 
#' @details Plots the factors in a barplot.
#'
plot_one_factor = function(f, pve, k, f_labels = NA, y_lab = "factor values") {
    P = length(f)
    if (any(is.na(f_labels))) {
        f_dat <- data.frame(variable = 1:P, Factor = f, sign.f = factor(sign(f)), hjust = factor(sign(f)))
        
        plot_f = ggplot2::ggplot(f_dat, aes(x = variable, y = Factor, fill = sign.f)) + geom_bar(stat = "identity", 
            width = 0.5) + scale_fill_manual(values = c("blue", "red")) + theme_minimal() + labs(title = paste("factor", 
            k, "with PVE=", round(pve, 3)), y = y_lab)
    } else {
        f_dat <- data.frame(variable = 1:P, Factor = f, sign.f = factor(sign(f)), variablenames = f_labels, hjust = factor(sign(f)))
        
        # 120% lim
        range_f = max(f) - min(f)
        upper_f = max(f, 0) + 0.15 * range_f
        lower_f = min(f, 0) - 0.15 * range_f
        
        plot_f = ggplot2::ggplot(f_dat, aes(x = variable, y = Factor, label = variablenames, fill = sign.f)) + geom_bar(stat = "identity", 
            width = 0.5) + geom_text(size = 2.75, angle = 90, hjust = as.character(f_dat$hjust), nudge_y = sign(f_dat$Factor) * 
            0.1 * mean(abs(f_dat$Factor))) + scale_fill_manual(values = c("blue", "red")) + ylim(lower_f, upper_f) + 
            theme_minimal() + labs(title = paste("factor", k, "with PVE=", round(pve, 3)), y = y_lab)
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
#' @return A ggplot plot object.
#' 
#' @export
#' 
flash_plot_pve = function(f, main = "Scree plot of PVE for each factor") {
  pve = flash_get_pve(f)
  pve_dat = data.frame(factor_index = seq(1, flash_get_k(f)), PVE = pve)
  p <- ggplot2::ggplot(pve_dat, aes(factor_index, PVE)) +
    geom_point(size = 4) + geom_line(linetype = "dotdash") + 
    labs(title = main)
  return(p)
}
