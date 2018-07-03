
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
            colour = fx.parse_column_names(color_by),
            shape = fx.parse_column_names(shape_by),
            size = fx.parse_column_names(size_by))
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
    
        c.df <- fx.get_circle_df(df_to_plot, separate_by = circle_by,  keep_quant= which_quantile)

        plot_out <- plot_out + ggalt::geom_encircle(data = c.df,
            aes_string(x = X, y = Y, fill = fx.parse_column_names(circle_by)),
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
            plot_out <- fx.resolve_plot_colors(plot_out, df_to_plot[,color_by],
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
        plot_out <- fx.resolve_plot_colors(plot_out,
            df_to_plot[,circle_by],
            colour_by_name = circle_by,
            fill = TRUE)
        ## ignore alpha for the fill legend
        plot_out <- plot_out + guides(fill = guide_legend(override.aes = list(alpha = 1)))
    }
  
    ##--add sensible shapes-------------------------------------------------------
    if(!is.null(shape_by)){
        plot_out <- fx.resolve_plot_shapes(plot_out, df_to_plot[,shape_by])
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



#' @name fx.get_circle_df 
#' @title Make a data.frame for encircling points in a reduced dimension plot
#' @description
#'
#' For the encircling, we will use the values of X and Y to determine
#' the coordinates for the circle that will either capture all values, or those
#' within the quantile defined by \code{which_quant}. 
#' Multiple circles may be drawn based on groups provided via \code{separate_by}.
#' 
#' @details  Sometimes, it is useful to remove outlier coordinates before drawing
#' the circles. If \code{keep_quant= 1}, the output will be exactly the content 
#' of \code{in.df}.
#' \code{keep_quant} should be a number between 1 and 0.5, usually you will want
#' to use something between 0.95 and 1. 
#' 
#' @param in.df data.frame that is used for making the original plot.
#' @param X name (or integer index) of the column of in.df that contains the 
#' values to be used for determining the x axis coordinates (default: "x_axs").
#' @param Y name (or integer index) of the column that contains the values to be 
#' used for determining the y axis coordinates (default: "y_axs".,
#' @param separate_by factor that will be used to split the values of X and Y
#' into distinct groups (e.g., by "condition"); if NULL (default), \code{which_quant}
#' @param keep_quant a number between 0.5 and 1 to define the max. quantile for
#' any value of \code{X} or \code{Y} that will be kept for defining the circle
#' coordinates. The default is to capture all points that belong to one instance
#' of \code{separate_by}, i.e. \code{keep_quant= 1}. It might be useful to set it
#' to less than 1 (e.g., 0.99) if you don't want the circle to be blown up by
#' individual outliers. 0.99 would ensure that only those points between the 99th
#' percentile and the 0.01st percentile (1 - 0.99) would be kept.
#' 
#' @details Integer indeces for X and Y must be given in proper integer notation
#' (e.g., 1L) so that \code{is.integer(X)} passes.
#'
#' @return data.table, most likely a subset of the in.df, depending on the setting
#' of \code{which_quant}. This data.table can be used with ggalt::geom_encircle.
#'
#' @seealso \code{\link{plot_reduced_dim.df}},
#' \code{\link{scABC2::generate.DimRedPlot}}
#' @importFrom stats quantile
#' 
#' @examples \dontrun{
#' c.dt <- fx.get_circle_df(df_to_plot, fill_by = circle_by, keep_quant= which_quantile)
#' plot_out <- plot_out +  ggalt::geom_encircle(data = c.dt,
#'                                              aes_string(fill = fx.parse_column_names(circle_by)),
#'                                              s_shape=0.9, expand=0.07, alpha = .1)
#' }
#' 
fx.get_circle_df <- function(in.df, X = "x_axs", Y = "y_axs", separate_by = NULL, keep_quant= 1){
  
    check_columns(separate_by, in.df,"data.frame for plotting (df_to_plot)", "fx.get_circle_df")
  
    comps <- fx.get_axes_val_names(in.df, X = X, Y = Y)
  
    if(keep_quant < 1 & keep_quant >= 0.5){
        pdt <- as.data.table(in.df)
        pdt[, x.high := quantile( get(comps[1]), probs = keep_quant), by = separate_by]
        pdt[, y.high := quantile( get(comps[2]), probs = keep_quant), by = separate_by]
        pdt[, x.low  := quantile( get(comps[1]), probs = 1 - keep_quant), by = separate_by]
        pdt[, y.low  := quantile( get(comps[2]), probs = 1 - keep_quant), by = separate_by]
    
        out.dt <- pdt[get(comps[1]) >= x.low & get(comps[1]) <= x.high & get(comps[2]) >= y.low & get(comps[2]) <= y.high, 
            colnames(in.df), with = FALSE]
        out.df <- as.data.frame(out.dt)
    }else{
        if(keep_quant != 1){
            message(paste0("The value you gave for keep_quant (", keep_quant, ") was not between 0.5 and 1 and was therefore ignored."))
        }
        out.df <- in.df
    }
  
    if(!is.null(separate_by)){
        if( !(all(is.character(out.df[,separate_by]))) | is.null(levels(out.df[,separate_by])) ){
            #warning(paste("The parameter you chose for the encircling is going to be turned into a factor."))
            out.df[, separate_by] <- factor(out.df[,separate_by])
        }
    }
  
    return(out.df)
  
}



#' @name fx.get_axes_val_names
#' @title Extract the column names that correspond to the X and Y axis
#' @description
#'
#' Very simple function that's basically just an error check for
#' whether the regex we think should be extracting the names of columns corresponding
#' to X and Y values work as expected: there should only be two values returned.
#' 
#' @param in.df Data.frame whose names() will be used to extract the instance that
#' matches \code{X} and \code{Y}
#' @param X name or integer index of the column that will be used for X axis
#' values. The default is "x_axs".
#' @param Y name or integer index of the column that will be used for Y axis
#' values. The default is "y_axs".
#' @return vector of names that correspond to the values where X_regex and Y_regex
#' match.
#' 
#' @examples \dontrun{
#' 
#' ## using column names
#' fx.get_axes_val_names(test.df, X = "PC1", Y = "PC2")
#' 
#' # using integer indeces
#' fx.get_axes_val_names(test.df, X = 1L, Y = 2L)
#' 
#' }
fx.get_axes_val_names <- function(in.df, X = "x_axs", Y = "y_axs"){
  
    ## > test.df <- data.frame(x_axs = c(1:4), y_axs = c(6:9), 
    ##                         value = c(11:14), group = c("a","b","c","d"))
    ## > fx.get_axes_val_names(test.df)
    ## [1] "x_axs" "y_axs"
    ## > fx.get_axes_val_names(test.df, X = 1L,Y = 2L)
    ## [1] "x_axs" "y_axs"
  
    if(is.integer(c(X, Y))){
        use_ints <- TRUE
    }else{
        use_ints <- FALSE
        check_columns(c(X, Y), in.df, "data.frame for plotting (df_to_plot)", "fx.get_circle_df")
    }
  
    if(use_ints){
        axs_names <- names(in.df)[c(X,Y)]
    }else{
        axs_names <- c(X, Y)
    } 
  
    return(axs_names)
}



#' @name check_columns
#' @title Check for the presence of named columns
#' @description
#'
#' Check for the presence of named columns
#'
#' @param which_names Column names to check, e.g. c("peaks", "genes")
#' @param input The names() of this will be checked.
#' @param input_name Name of the input, e.g. "start_params" just to make it
#' identifiable if an error message is returned.
#' @param function_name Name of function for which this test is carried out.
#' Again, just to make the returning error message a bit more readable. This
#' helps with debugging if you're wrapping many functions within each other.
#'
#' @return Returns an error and stop signal if entries of \code{which_names}
#'  are missing in the \code{input}.
#' @examples
#' \dontrun{
#' check_columns( c("cells", "sample", "condition"),
#'                long_df, "long_df", "plot_profile")
#' }
check_columns <- function(which_names, input, input_name, function_name){
    check <- which_names %in% names(input)
    if(!all(check)){
        stop(paste0("The input (", input_name , ") supplied to ", function_name, " is missing the following named columns: ",   paste(which_names[!check], collapse = ", ") ))
    }
}




#' @name fx.parse_column_names
#' @title Parse column names
#' @description
#'
#' Parse column names
#'
#' @details There are numerous ways to access the values
#' stored in the columns of a \code{data.frame}. The one
#' using a dollar sign (e.g. \code{df$this_column}) is
#' vulnerable to special characters, such as "-", spaces etc.
#' Since the default method of \code{ggplot2} is to use the
#' dollar sign approach, we need to ensure that special
#' characters, which are often part of gene names, are taken
#' care of. This function simply wraps those names into
#' single quotes.
#'
#' @param cn column name
#' @return column name enclosed in single quotes
#' 
fx.parse_column_names <- function(cn){
    if(is.null(cn)){
        cnp <- cn
    }else{
        cnp <- paste0("`", cn, "`")
    }
  
    return(cnp)
}




#' @name fx.resolve_plot_colors
#' @title Get nice plotting color schemes for very general color variables
#' @description
#'
#' Wrapper around \code{fx.get_palette_ABC} that checks the numbers
#' and types of coloring variables that are being used and tries to return a 
#' sensible color scheme. The colors are pre-defined in \code{fx.get_palette_ABC}.
#' 
#' @details This function is based on a very similar function in the scater package.
#' 
#' @param plot_out ggplot2 object
#' @param colour_by vector of values that determine the coloring of \code{plot_out}
#' @param colour_by_name string indicating the title/name for \code{colour_by}, e.g. the name of a gene
#' @param fill Boolean, default: \code{FALSE}
#' 
#' @return \code{ggplot2} object with adjusted coloring scheme
#' 
#' @seealso \code{fx.get_palette_ABC}
#'
fx.resolve_plot_colors <- function(plot_out, colour_by, colour_by_name, fill = FALSE) {
    ## if the colour_by object is NULL, return the plot_out object unchanged
    if ( is.null(colour_by) ){
        return(plot_out)
    }
    ## Otherwise, set a sensible colour scheme and return the plot_out object
    leg_title <- colour_by_name
  
    if (fill) { ## routine for fill
        if (is.numeric(colour_by)) {
            plot_out <- plot_out + viridis::scale_fill_viridis(name = leg_title)
        } else {
            nlevs_colour_by <- nlevels(as.factor(colour_by))
            if (nlevs_colour_by <= 14) {
                plot_out <- plot_out + scale_fill_manual(
                    values = fx.get_palette_ABC("paired_pal"),
                    name = leg_title)
            } else {
                if (nlevs_colour_by > 14 && nlevs_colour_by <= 20) {
                    plot_out <- plot_out + scale_fill_manual(
                        values = fx.get_palette_ABC("tableau20"),
                        name = leg_title)
                } else {
                    plot_out <- plot_out + viridis::scale_fill_viridis(name = leg_title, discrete = TRUE)
                }
            }
        }
    } else { ## routine for color
        if (is.numeric(colour_by)){
            plot_out <- plot_out + viridis::scale_color_viridis(name = leg_title)
        } else {
            nlevs_colour_by <- nlevels(as.factor(colour_by))
            if (nlevs_colour_by <= 14) {
                plot_out <- plot_out + scale_colour_manual(
                    values = fx.get_palette_ABC("paired_pal"),
                    name = leg_title)
            } else {
                if (nlevs_colour_by > 14 && nlevs_colour_by <= 20) {
                    plot_out <- plot_out + scale_colour_manual(
                        values = fx.get_palette_ABC("tableau20"),
                        name = leg_title)
                } else {
                    plot_out <- plot_out + viridis::scale_color_viridis(name = leg_title, discrete = TRUE)
                }
            }
        }
    }
    plot_out
}




#' @name fx.get_palette_ABC
#' @title Color palettes
#' @description 
#' 
#' This function simply supplies pre-defined color palettes.
#'
#' @details Based on \code{scater}'s defaults, but with significant changes to the
#' standard colors that were bing used
#' 
#' @return vector of color names
#' 
#' @seealso \code{fx.resolve_plot_colors}
fx.get_palette_ABC = function(palette_name) {
    switch(palette_name,
        tableau10medium = c("#34B20D","#FFAE18", "#be2f00","#73FFC3" , "#0D14B2",
            "#8EB20E",  "#FF81DE", "#FFFF00", "#0CCC9C", "#656BB2"),
        paired_pal = c("#A6CEE3","limegreen", "grey30", "grey80", "#1F78B4" , "#33A02C", "#FB9A99" ,"#E31A1C" ,"#FDBF6F", "#FF7F00" ,"#CAB2D6" ,"#6A3D9A" ,"#FFFF99" ,"#B15928"), # from RColorBrewer::brewer.pal(12,  "Paired")
        tableau20 = c("#34B20D","#FFAE18", "#be2f00","#73FFC3" , "#0D14B2",
            "#8EB20E",  "#FF81DE", "#FFFF00", "#0CCC9C", "#656BB2", 
            c( "#FF10FC" , "#3E68FF", "#8B9440", "#7F85C3", "#FF85C3", "#F52306" ,"#FFD2BD", "#25FFED", "black", "#FFEC83" )),
        colorblind10 = c("#006BA4", "#FF800E", "#ABABAB", "#595959",
            "#5F9ED1", "#C85200", "#898989", "#A2C8EC","#FFBC79", "#CFCFCF"),
        trafficlight = c("#B10318", "#DBA13A", "#309343", "#D82526",
            "#FFC156", "#69B764", "#F26C64", "#FFDD71", "#9FCD99"),
        purplegray12 = c("#7B66D2", "#A699E8", "#DC5FBD", "#FFC0DA",
            "#5F5A41", "#B4B19B", "#995688", "#D898BA","#AB6AD5", "#D098EE", "#8B7C6E", "#DBD4C5"),
        bluered12 = c("#2C69B0", "#B5C8E2", "#F02720", "#FFB6B0", "#AC613C",
            "#E9C39B", "#6BA3D6", "#B5DFFD", "#AC8763", "#DDC9B4","#BD0A36", "#F4737A"),
        greenorange12 = c("#32A251", "#ACD98D", "#FF7F0F", "#FFB977",
            "#3CB7CC", "#98D9E4", "#B85A0D", "#FFD94A","#39737C", "#86B4A9", "#82853B", "#CCC94D"),
        cyclic = c("#1F83B4", "#1696AC", "#18A188", "#29A03C", "#54A338",
            "#82A93F", "#ADB828", "#D8BD35", "#FFBD4C", "#FFB022",
            "#FF9C0E", "#FF810E", "#E75727", "#D23E4E", "#C94D8C","#C04AA7", "#B446B3", "#9658B1", "#8061B4", "#6F63BB"))
}







