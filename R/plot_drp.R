
#' @name get_reducedDimPlot.sce
#' @title Visualizing tSNE and PCA results
#' @description
#' 
#' This will generate a ggplot starting from a \code{DimRedPlot} object.
#' See \code{\link{get_reducedDimPlot.sce}} to start from a \code{SingleCellExperiment}
#' object.
#'
#' @param drp_object a \code{\link{DimRedPlot}} object
#' @param X name of the drp_object column that contains the values for the x axis
#' @param Y name of the drp_object column that contains the values for the y axis
#' @param which_quantile default: 1
#' @param color_by specify the column of \code{drp_object$plot_data} whose
#' entries will be used to assign either discrete or continous color schemes.
#' @param shape_by the factor to assign the shape of the points to; max. 10 levels
#' are accepted for \code{shape_by}. If you want to specify the (same) shape for
#' \emph{all} points, use an integer (default shape is 19).
#' @param size_by the factor to assign the size of the points to;
#' default: \code{NULL}. If you want to specify the (same) size for \emph{all}
#' points, use one integer (the default size is 3).
#' @param circle_by a factor (e.g. "condition") that will be used to draw circles
#' around the points belonging to the same factor instance.
#' @param ignore_drp_labels If set to TRUE, \emph{all} entries of
#' \code{drp_object$factors} will be ignored. If you want to ignore only selected
#' labels, indicate their name (e.g. c("color_by", "exprs_val_type")).
#' The default setting (\code{FALSE}) will make use of the information from 
#' the drp_object \emph{unless} specified here via \code{color_by}, \code{shape_by},
#' \code{size_by}, \code{circle_by}.
#' @param theme_size font size for the plot
#' @param legend choose wether there should be no legend ("none"), a legend for 
#' every factor or whether the legend will only show up for those factors that
#' have more than 1 level ("auto", default setting).
#' @param alpha define the opacity of the points, default: .65
#' @param remove_rug Define whether the bars at both axes that
#' represent the two 1d marginal distributions should be removed. Default: FALSE.
#' @param set_colors Set to \code{FALSE} if you want to add your own color scheme
#' to the resulting plot. The default behavior is to try to automatically assign
#' an optimized coloring scheme. Default: \code{TRUE}.
#' @param set_fill_colors Set to \code{FALSE} if you want to add your own color
#' scheme to the resulting plot \emph{for the background circles}. The default
#' behavior is to try to automatically assign an optimized coloring scheme.
#'
#' @seealso \code{\link{get_reducedDimPlot.sce}}, which generates the DimRedPlot
#' object that will be used here.
#' 
#' @return \code{ggplot2} object (a plot)
#' 
#' @examples
#' data("tiny_sce")
#' pdrp <- get_reducedDimPlot.sce(tiny_sce, which_reddim = "PCA",
#'                               which_pcs = c(2:3),
#'                               color_by = names(colData(tiny_sce))[1],
#'                               dim_red_type = "PCA",
#'                               add_cell_info = names(colData(tiny_sce)),
#'                               exprs_values = "log10medScaled")
#' # color based on the gene used during the first step
#' plt.DimRedPlot(pdrp)
#' 
#' # no color
#' plt.DimRedPlot(pdrp, ignore_drp_labels = "color_by")
#'   
#'# color based on condition, avoid "log10medScaled" being printed in the legend
#' plt.DimRedPlot(pdrp, color_by = "condition", ignore_drp_labels = "exprs_val_type")
#'
#' # color based on barcode, encircle points that belong to the same condition 
#' plt.DimRedPlot(pdrp, color_by = "barcode", ignore_drp_labels = "exprs_val_type", circle_by = "condition")
#' 
#' @export
#' 
plt.DimRedPlot <- function(drp_object, 
    X="x_axs", 
    Y="y_axs",
    color_by=NULL, 
    shape_by=NULL, 
    size_by=NULL, 
    circle_by = NULL, 
    which_quantile = 1,
    ignore_drp_labels = FALSE,
    theme_size = 14, 
    legend = "auto", 
    alpha = .75,
    remove_rug = FALSE,
    set_colors = TRUE, 
    set_fill_colors = TRUE) {
  
    ##--set the defaults for the plot---------------------------------------------
    if(is.numeric(shape_by) & length(shape_by) == 1){
        shape_val <- shape_by
        shape_by <- NULL
        ignore_drp_labels <- c(ignore_drp_labels, "shape_by")
    }else{
        shape_val <- 19
    }
  
    if(is.numeric(size_by) & length(size_by) == 1 ){
        size_val <- size_by
        size_by <- NULL
        ignore_drp_labels <- c(ignore_drp_labels, "size_by")
    }else{
        size_val <- 3
    }
  
    ggplot2::update_geom_defaults("point", list(shape = shape_val, size = size_val))
   
    ##--extract basic plot information from DRP object----------------------------
    color_by <- fx.set_factors(drp_object, fctr = color_by, 
        fctr_type = "color_by", ignore_drp_labels)
    shape_by <- fx.set_factors(drp_object, fctr = shape_by, 
        fctr_type = "shape_by", ignore_drp_labels)
    size_by <- fx.set_factors(drp_object, fctr = size_by, 
        fctr_type = "size_by", ignore_drp_labels)
    circle_by <- fx.set_factors(drp_object, fctr = circle_by, 
        fctr_type = "circle_by", ignore_drp_labels)
  
    ##--extract the data.frame needed for plotting--------------------------------
    df_to_plot <- drp_object$plot_data
    df_to_plot$key <- row.names(df_to_plot)
  
    ##--check legend argument ----------------------------------------------------
    legend <- match.arg(legend, c("auto", "none", "all"), several.ok = FALSE)

    ## if only one level for the variable, set legend to NULL
    if ( legend == "auto" ) {
        if (length(unique(df_to_plot$color_by)) == 1){
            color_by <- NULL
        }
        if (length(unique(df_to_plot$shape_by)) == 1){
            shape_by  <- NULL
        }
        if (length(unique(df_to_plot$size_by)) == 1){
            size_by <- NULL
        }
        if (length(unique(df_to_plot$circle_by)) == 1){
            circle_by <- NULL
        }
    }
  
    ##--generate base plot, apply colour_by, shape_by and size_by variables ------
    ## (NULL will simply be ignored by ggplot)
    if( !is.null(color_by) && !(color_by %in% names(df_to_plot)) ){color_by <- NULL}
    if( !is.null(shape_by) && !(shape_by %in% names(df_to_plot)) ){shape_by <- NULL}
    if( !is.null(size_by)  && !(size_by  %in% names(df_to_plot)) ){size_by <- NULL}
  
    if(is.numeric(df_to_plot[,shape_by])){
        df_to_plot[, shape_by] <- factor( df_to_plot[, shape_by])
    }
  
    ##==BASE PLOT ================================================================
    plot_out <- ggplot(df_to_plot,
        aes_string(x = X, y = Y, key = "key",
            colour = ABCutilities::fx.parse_column_names(color_by),
            shape = ABCutilities::fx.parse_column_names(shape_by),
            size = ABCutilities::fx.parse_column_names(size_by))
        )
  
    ## return(plot_out)
    ##  stop()
    if(!remove_rug){
        plot_out <- plot_out + geom_rug(colour = "gray20", alpha = 0.5)
    }
   
    ##--add x and y labels--------------------------------------------------------
    ld <- drp_object$label_data
  
    if( !is.null(ld$variance_pct) & !("variance_pct" %in% ignore_drp_labels) ){
        x_lab <- paste0(ld$x_lab, ": ", round(ld$variance_pct[1] * 100), "% variance")
        y_lab <- paste0(ld$y_lab, ": ", round(ld$variance_pct[2] * 100), "% variance")
    }else{
        x_lab <- ld$x_lab
        y_lab <- ld$y_lab
    }
  
    plot_out <- plot_out + xlab(x_lab) +  ylab(y_lab)
  
    ##--add background circles (do not move below geom_point!)--------------------
    if(!is.null(circle_by)){
    
        c.df <- ABCpca::fx.get_circle_df(df_to_plot, separate_by = circle_by,  keep_quant= which_quantile)

        plot_out <- plot_out + ggalt::geom_encircle(data = c.df,
            aes_string(x = X, y = Y, fill = ABCutilities::fx.parse_column_names(circle_by)),
            ## color = NULL, shape = NULL, size = NULL)
            inherit.aes = FALSE,
            s_shape = 0.9, 
            expand = 0.07, 
            alpha = .2)
    }
  
    ##--add points (leave this _after_ the background circles so that the points--
    ##--are plotted on top of the circles ----------------------------------------
    plot_out <- plot_out +  geom_point(alpha = alpha)
  
    ##--assign a sensible color scheme based on the features of color_by----------
    if(!is.null(color_by)){
        if(set_colors){
            plot_out <- ABCutilities::fx.resolve_plot_colors(plot_out, df_to_plot[,color_by],
                colour_by_name = color_by,
                fill = FALSE)
        }
   
        if(!is.null(ld$exprs_val_type) & !("exprs_val_type" %in% ignore_drp_labels)){
            color_label <- paste0(color_by, "\n(", ld$exprs_val_type, ")")
        }else{
            color_label <- color_by
        } 
        ## ignore alpha for the color legend
        plot_out <- plot_out + guides(color = guide_legend(override.aes = list(alpha = 1), title = color_label))
    }  
  
    ##--assign a sensible color scheme for the background circle------------------
    if(!is.null(circle_by) & set_fill_colors){
        plot_out <- ABCutilities::fx.resolve_plot_colors(plot_out,
            df_to_plot[,circle_by],
            colour_by_name = circle_by,
            fill = TRUE)
        ## ignore alpha for the fill legend
        plot_out <- plot_out + guides(fill = guide_legend(override.aes = list(alpha = 1)))
    }
  
    ##--add sensible shapes-------------------------------------------------------
    if(!is.null(shape_by)){
        plot_out <- ABCutilities::fx.resolve_plot_shapes(plot_out, df_to_plot[,shape_by])
        ## ignore alpha for the scale legend
        plot_out <- plot_out +  guides(shape = guide_legend(override.aes = list(alpha = 1)))
    }
  
    ##--define plotting theme ----------------------------------------------------
    if ( requireNamespace("cowplot", quietly = TRUE) ){
        plot_out <- plot_out + cowplot::theme_cowplot(theme_size)
    } else{
        plot_out <- plot_out + theme_bw(theme_size)
    }
  
    ##--remove legend if so desired ----------------------------------------------
    if ( legend == "none" ){
        plot_out <- plot_out + theme(legend.position = "none")
    }
  
    ##--return plot---------------------------------------------------------------
    return(plot_out)
}
