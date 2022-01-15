
plot_genes_in_pseudotime2 <- function (cds_subset,module_mat, min_expr = NULL, scaled = F, cell_size = 0.75, nrow = NULL, 
          ncol = 1, panel_order = NULL, color_cells_by = "pseudotime", 
          trend_formula = "~ splines::ns(pseudotime, df=3)", label_by_short_name = TRUE, 
          vertical_jitter = NULL, horizontal_jitter = NULL, returncurve = F) 
{
  if(scaled){
    minagg = abs(min(module_mat))
    module_mat = module_mat + minagg +1 ##getting it above 0 for the model matching to work (i don't really know why)
    min_expr <- -minagg
  }
  
  cds_subset <- cds_subset[1,]
  
  assertthat::assert_that(methods::is(cds_subset, "cell_data_set"))
  tryCatch({
    pseudotime(cds_subset)
  }, error = function(x) {
    stop(paste("No pseudotime calculated. Must call order_cells first."))
  })
  colData(cds_subset)$pseudotime <- pseudotime(cds_subset)
  if (!is.null(min_expr)) {
    assertthat::assert_that(assertthat::is.number(min_expr))
  }
  assertthat::assert_that(assertthat::is.number(cell_size))
  if (!is.null(nrow)) {
    assertthat::assert_that(assertthat::is.count(nrow))
  }
  assertthat::assert_that(assertthat::is.count(ncol))
  assertthat::assert_that(is.logical(label_by_short_name))
  if (label_by_short_name) {
    assertthat::assert_that("gene_short_name" %in% names(rowData(cds_subset)), 
                            msg = paste("When label_by_short_name = TRUE,", 
                                        "rowData must have a column of gene", "names called gene_short_name."))
  }
  assertthat::assert_that(color_cells_by %in% c("cluster", 
                                                "partition") | color_cells_by %in% names(colData(cds_subset)), 
                          msg = paste("color_cells_by must be a column in the", 
                                      "colData table."))
  if (!is.null(panel_order)) {
    if (label_by_short_name) {
      assertthat::assert_that(all(panel_order %in% rowData(cds_subset)$gene_short_name))
    }
    else {
      assertthat::assert_that(all(panel_order %in% row.names(rowData(cds_subset))))
    }
  }
  
  assertthat::assert_that(nrow(rowData(cds_subset)) <= 100, 
                          msg = paste("cds_subset has more than 100 genes -", 
                                      "pass only the subset of the CDS to be", "plotted."))
  assertthat::assert_that(methods::is(cds_subset, "cell_data_set"))
  assertthat::assert_that("pseudotime" %in% names(colData(cds_subset)), 
                          msg = paste("pseudotime must be a column in", "colData. Please run order_cells", 
                                      "before running", "plot_genes_in_pseudotime."))
  if (!is.null(min_expr)) {
    assertthat::assert_that(assertthat::is.number(min_expr))
  }
  assertthat::assert_that(assertthat::is.number(cell_size))
  assertthat::assert_that(!is.null(size_factors(cds_subset)))
  if (!is.null(nrow)) {
    assertthat::assert_that(assertthat::is.count(nrow))
  }
  assertthat::assert_that(assertthat::is.count(ncol))
  assertthat::assert_that(is.logical(label_by_short_name))
  if (label_by_short_name) {
    assertthat::assert_that("gene_short_name" %in% names(rowData(cds_subset)), 
                            msg = paste("When label_by_short_name = TRUE,", 
                                        "rowData must have a column of gene", "names called gene_short_name."))
  }
  assertthat::assert_that(color_cells_by %in% c("cluster", 
                                                "partition") | color_cells_by %in% names(colData(cds_subset)), 
                          msg = paste("color_cells_by must be a column in the", 
                                      "colData table."))
  if (!is.null(panel_order)) {
    if (label_by_short_name) {
      assertthat::assert_that(all(panel_order %in% rowData(cds_subset)$gene_short_name))
    }
    else {
      assertthat::assert_that(all(panel_order %in% row.names(rowData(cds_subset))))
    }
  }
  assertthat::assert_that(nrow(rowData(cds_subset)) <= 100, 
                          msg = paste("cds_subset has more than 100 genes -", 
                                      "pass only the subset of the CDS to be", "plotted."))
  ##end of assertions
  
  f_id <- NA
  Cell <- NA
  
  cds_subset@colData@listData[["Size_Factor"]] <- rep(1.0,length(cds_subset@colData@listData[["Size_Factor"]]))
  names(cds_subset@colData@listData[["Size_Factor"]]) <- names(cds_subset@colData@listData[["pseudotime"]])
  cds_subset@assays@data@listData[["counts"]]@x <- module_mat
  cds_subset@assays@data@listData[["counts"]]@Dimnames[[1]] <- "module"
  cds_subset@rowRanges@elementMetadata@listData[["gene_short_name"]] <- "module"
  cds_subset@assays@data@listData[["counts"]]@i <- as.integer(rep(0, length(cds_subset@assays@data@listData[["counts"]]@x)))
  cds_subset@assays@data@listData[["counts"]]@p <- as.integer(0:length(cds_subset@assays@data@listData[["counts"]]@x))
  
  cds_subset = cds_subset[, is.finite(colData(cds_subset)$pseudotime)]
  cds_exprs <- SingleCellExperiment::counts(cds_subset)
  cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/size_factors(cds_subset))
  cds_exprs@x <- module_mat ##put here the data you need
  cds_exprs <- reshape2::melt((as.matrix(cds_exprs))) #melt(round(...))
  if (is.null(min_expr)) {
    min_expr <- 0
  }
  colnames(cds_exprs) <- c("f_id", "Cell", "expression")
  cds_colData <- colData(cds_subset)
  cds_rowData <- rowData(cds_subset)
  cds_exprs <- merge(cds_exprs, cds_rowData, by.x = "f_id", by.y = "row.names")
  cds_exprs <- as.data.frame(cds_exprs)
  cds_exprs <- merge(cds_exprs, cds_colData, by.x = "Cell", by.y = "row.names")
  cds_exprs <- as.data.frame(cds_exprs)
  cds_exprs$adjusted_expression <- cds_exprs$expression
  if (label_by_short_name == TRUE) {
    if (is.null(cds_exprs$gene_short_name) == FALSE) {
      cds_exprs$feature_label <- as.character(cds_exprs$gene_short_name)
      cds_exprs$feature_label[is.na(cds_exprs$feature_label)] <- cds_exprs$f_id
    }else {
      cds_exprs$feature_label <- cds_exprs$f_id
    }
  }else {
    cds_exprs$feature_label <- cds_exprs$f_id
  }
  cds_exprs$f_id <- as.character(cds_exprs$f_id)
  cds_exprs$feature_label <- factor(cds_exprs$feature_label)
  new_data <- data.frame(pseudotime = colData(cds_subset)$pseudotime)
  model_tbl = fit_models(cds_subset, model_formula_str = trend_formula)
  model_expectation <- model_predictions(model_tbl, new_data = colData(cds_subset))
  colnames(model_expectation) <- colnames(cds_subset)
  expectation <- plyr::ddply(cds_exprs, plyr::.(f_id, Cell), 
                             function(x) {
                               data.frame(expectation = model_expectation[x$f_id, 
                                                                          x$Cell])
                             })
  cds_exprs <- merge(cds_exprs, expectation)
  cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
  cds_exprs$expectation[cds_exprs$expectation < min_expr] <- min_expr
  if (!is.null(panel_order)) {
    cds_exprs$feature_label <- factor(cds_exprs$feature_label, 
                                      levels = panel_order)
  }
  if(scaled){
    cds_exprs$expression <- cds_exprs$expression - minagg-1
    cds_exprs$expectation <- cds_exprs$expectation - minagg-1
  }
  
  q <- ggplot(aes(pseudotime, expression), data = cds_exprs)
  if (!is.null(color_cells_by)) {
    q <- q + geom_point(aes_string(color = color_cells_by), 
                        size = I(cell_size), position = position_jitter(horizontal_jitter, 
                                                                        vertical_jitter))
    if (class(colData(cds_subset)[, color_cells_by]) == 
        "numeric") {
      q <- q + viridis::scale_color_viridis(option = "C")
    }
  }else {
    q <- q + geom_point(size = I(cell_size), position = position_jitter(horizontal_jitter, 
                                                                        vertical_jitter))
  }
  q <- q + geom_line(aes(x = pseudotime, y = expectation), data = cds_exprs)
  if(!scaled){
    q <- q + scale_y_log10() + facet_wrap(~feature_label, nrow = nrow, ncol = ncol, scales = "free_y")
  }
  if (min_expr < 1) {
    q <- q + expand_limits(y = c(min_expr, 1))
  }
  q <- q + ylab("Expression")
  q <- q + xlab("pseudotime")
  #q <- q + monocle_theme_opts()

  if(returncurve){
    return(cds_exprs[,c("pseudotime","expectation","expression")])
  }
  
  q
}





plot_genes_of_module = function(cds, genes, scaled = F, logscale = F, legend = F, trend_formula = "~ splines::ns(pseudotime, df=4)", min.expr = 1){
  plot <- ggplot()
  if(logscale & scaled){
    stop("'logscale' and 'scaled' cannot be True at the same time")
  }
  for(i in 1:length(genes)){
    
    genedata <- exprs(cds)[genes[i],]/size_factors(cds)
    if (var(genedata)==0){
      next
    }
    if(scaled){
      genedata <- scale(genedata)
      genedata <- as.numeric(genedata)
    }
    newcurve = plot_genes_in_pseudotime2(cds, min_expr = 0,module_mat = genedata, scaled = scaled,returncurve = T, trend_formula = trend_formula )
    if(i ==1){
      curves = newcurve
    }else{
      curves <- cbind(curves,newcurve$expectation)
    }
    #print(genes[i])
    #if(any(newcurve$expectation>12)){
    #  print(genes[i])
    #}
    if(logscale){
      newcurve$expectation <- newcurve$expectation +min.expr
    }
    newcurve$name <- genes[i]
    if(legend){
      plot <- plot + geom_line(aes(x = pseudotime, y = expectation, color = name), data = newcurve)
    }else{
      plot <- plot + geom_line(aes(x = pseudotime, y = expectation), data = newcurve)
    }
    
  }
  if(logscale){
    plot <- plot + scale_y_log10()
  }
  #curves
  plot + theme_minimal()
}

