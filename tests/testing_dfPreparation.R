## fx.get_reducedDim_df.sce --------------------------------------------------
test_that("S3 object (DimRedPlot) generated for the plot looks as expected",
          {
            ## run on test data
            gn <- rownames(tiny_sce)[1]
            cdt1 <- names(colData(tiny_sce))[1]
            cdt2 <- names(colData(tiny_sce))[2]
            cdt3 <- names(colData(tiny_sce))[3]
            df2plot_from_sce <- get_reducedDim_df.sce(tiny_sce, which_reddim = "PCA",
                                                      color_by = gn, shape_by = cdt1,
                                                      size_by = cdt2, circle_by = cdt3)
          
            expect_that(df2plot_from_sce, is_a("DimRedPlot"))
            
            ## check output 1
            x <- df2plot_from_sce$plot_data
            expect_that(x, is_a("data.frame"))
            expect_match(names(x)[1], "x_axs")
            expect_match(names(x)[2], "y_axs")
            expect_match(names(x), paste(c(gn,cdt1, cdt2, cdt3), collapse = "|"), all = FALSE)
            
            ## check output 2
            xx <- df2plot_from_sce$label_data
            expect_that(xx, is_a("list"))
            expect_match( names(xx), paste(c("x_lab","y_lab","exprs_val_type"),collapse="|"), all = FALSE)
            
            
            })

