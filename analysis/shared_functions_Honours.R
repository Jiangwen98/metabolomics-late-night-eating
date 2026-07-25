# ---------------------------------------------------------------------------
# ANALYSIS FUNCTIONS
#
# Utility functions used by the analysis scripts in this repository: quality
# control filtering, transformation, scaling, outlier handling, and plotting.
#
# These functions originate in earlier analysis code from the Drum Laboratory,
# inherited from a former lab member together with the B&B project. They are
# not my own work. Tags in the headers below preserve the distinction made in
# the original source: "adapted from earlier lab code" marks the functions I
# modified for the requirements of this analysis; the others are as inherited.
#
# Included so the analysis scripts can be read and run. See ATTRIBUTION.md.
# ---------------------------------------------------------------------------


#' Plots Single-Variable Line Graphs Across Time (of Different Categories) with error bars -- adapted from earlier lab code
#'
#' @param qtc Dataframe. Dataframe containing the values of the Variables (should be single/median values), their categorical classifications and time (metadata).
#' @param colName_toSortBy String. Column name containing the different categories of data to apply contrasts on.
#' @param xAxisName String. Column name for the X-axis. Generally should be related to Time analysis - indicating the different time points.
#' @param mtbl_name String. Column name of variable to be plotted in that graph.
#' @param y_label String. Y-label for graph. Default is 'Log10(Counts)'.
#' @param x_label String. X-label for graph. Default is 'Time'.
#'
#' @return A ggplot2 plot of the line graph across time.
#' @export
#'
#' @examples plot_graph_line2(qtc, 'Classification', 'Time in Minutes', 'Biomarker 33')
plot_graph_line0 <- function(qtc, colName_toSortBy, xAxisName, mtbl_name, ymin, ymax, y_label = 'Counts', x_label = 'Time after meal (min)', labels,
                             palette = c("lightgreen", "skyblue2", "lightpink", "plum")){
  require(tidyverse)
  
  mtbl_name <- paste0('`', mtbl_name, '`')  ## We cannot eval-parse spaces or slashes!
  
  pd <- position_dodge(width = 0.4)
  
  p1 <- ggplot(data = qtc, aes(x = get(xAxisName), y = eval(parse(text = mtbl_name)),
                               group = get(colName_toSortBy))) +
    theme_minimal() +
    ggtitle(str_remove_all(mtbl_name, '`')) + theme(plot.title = element_text(size = 30, hjust = 0.5, face = 'bold'),
                                                    # axis.title.x = element_text(size = 14, face = 'bold', vjust = -1),
                                                    axis.text.x = element_text(size = 11, face = 'bold'),
                                                    axis.title.x = element_text(size = 18, face = 'bold'),
                                                    axis.title.y = element_text(size = 18, face = 'bold'),
                                                    axis.text.y = element_text(size = 11, face = 'bold'),
                                                    legend.position = 'right') +
    ylab(y_label) + xlab(x_label) + 
    geom_line(aes(linetype = Category, color = Category), lwd = 1.3, position = pd) + 
    geom_point(aes(color = Category, shape = Category), size = 4.5, position = pd) + 
    geom_errorbar(aes(ymin=eval(parse(text = ymin)), ymax=eval(parse(text = ymax)), color = Category), width = 0.1, linewidth = 0.7, position = pd) +
    scale_color_manual(values = palette) + ## Better use a palette instead of annual adjustment every time
    scale_shape_manual(values = c(15, 16, 17, 18))
  
  
  
  return(p1)
}

#' Plots Single-Variable Line Graphs Across Time (of Different Categories) with error bars -- adapted from earlier lab code
#'
#' @param qtc Dataframe. Dataframe containing the values of the Variables (should be single/median values), their categorical classifications and time (metadata).
#' @param colName_toSortBy String. Column name containing the different categories of data to apply contrasts on.
#' @param xAxisName String. Column name for the X-axis. Generally should be related to Time analysis - indicating the different time points.
#' @param mtbl_name String. Column name of variable to be plotted in that graph.
#' @param y_label String. Y-label for graph. Default is 'Log10(Counts)'.
#' @param x_label String. X-label for graph. Default is 'Time'.
#'
#' @return A ggplot2 plot of the line graph across time.
#' @export
#'
#' @examples plot_graph_line2(qtc, 'Classification', 'Time in Minutes', 'Biomarker 33')
plot_graph_line2 <- function(qtc, colName_toSortBy, xAxisName, mtbl_name, ymin, ymax, y_label = 'Counts', x_label = 'Time after meal (min)', labels,
                             palette = c("olivedrab3", "gold2", "darkmagenta", "lightgreen", "skyblue2", "lightpink", "plum")){
  require(tidyverse)
  
  mtbl_name <- paste0('`', mtbl_name, '`')  ## We cannot eval-parse spaces or slashes!
  
  pd <- position_dodge(width = 0.4)
  
  p1 <- ggplot(data = qtc, aes(x = get(xAxisName), y = eval(parse(text = mtbl_name)),
                               group = get(colName_toSortBy))) +
    theme_classic() +
    ggtitle(str_remove_all(mtbl_name, '`')) + theme(plot.title = element_text(size = 30, hjust = 0.5, face = 'bold'),
                                                    # axis.title.x = element_text(size = 14, face = 'bold', vjust = -1),
                                                    axis.text.x = element_text(size = 11, face = 'bold'),
                                                    axis.title.x = element_text(size = 18, face = 'bold'),
                                                    axis.title.y = element_text(size = 18, face = 'bold'),
                                                    axis.text.y = element_text(size = 11, face = 'bold'),
                                                    legend.position = 'right') +
    ylab(y_label) + xlab(x_label) + 
    geom_line(aes(linetype = Category, color = Category), lwd = 1.3, position = pd) + 
    geom_point(aes(color = Category, shape = Category), size = 4.5, position = pd) + 
    geom_errorbar(aes(ymin=eval(parse(text = ymin)), ymax=eval(parse(text = ymax)), color = Category), width = 0.1, linewidth = 0.3, position = pd) +
    scale_color_manual(values = palette) + ## Better use a palette instead of annual adjustment every time
    scale_shape_manual(values = c(15, 16, 17, 18))
  
  
  
  return(p1)
}