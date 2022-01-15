library(ggplot2)
library(rlang)
library(ggrepel)
library(cowplot)

LayeredFeaturePlot <- function (object, features = NULL, assay = NULL, 
                                reduction = "tsne", cells.use = NULL, pt.size = 1, 
                                colors = NULL, min.cutoff = NULL, max.cutoff = NULL, 
                                fig.ratio = 6, label = FALSE, repel = FALSE, bold = FALSE, 
                                background.color = "white", pt.brightness = 10) {
  ##### Arguments ###########
  # object: A Seuratv3 object
  # features: features (genes) to plot
  # assay: The assay that you want to plot. If not set, you'll need to 
  # provide assay key in the feature name to prevent duplicated feature names in different assays to cause issues.
  # reduction: The dimensionality reduction to visualize expression on
  #            (e.g., tsne (default), umap, or pca)
  # cells.use: A vector of cell barcodes to be plotted (default: all cells)
  # pt.size: The size of each point on the resulting plot
  # colors: Color names to visualize each feature. Must be as many as features
  #         (e.g., c("red", "blue", "green", "purple"))
  # min.cutoff: Log-normalized expression level lower than this will be 
  #             seen as 0. An absolute value (e.g., 3) or a percentile (e.g., 
  #             "q70" means 70th percentile)
  # max.cutoff: Log-normalized expression level higher than this will be capped.
  #             (e.g., if max.cutoff = 5, anything above 5 will be plotted as 5.)
  #             This can be useful if you have some extreme high expression in 
  #             something that you don't really care. A percentile is accepted 
  #             like min.cutoff
  # fig.ratio: The ratio of figure itself to figure legend (e.g., fig.ratio = 6 
  #            makes (the width of figure):(the width of legend) = 6:1).
  # label: If TRUE, label the clusters with active.idents
  # repel: (Only useful when label = TRUE) If TRUE, overlapped labels will nudge to prevent overlapping
  # bold: (Only useful when label = TRUE) If TRUE, the labels will be in bold
  # background.color: The background of the figure panel (default is light grey (grey92)).
  
  
  ##### Internal func #######
  exp_encolor <- function(data, features, pt.brightness) {
    # Adjust expression to fit 0 - max per feature
    # In a scale of 0 - 100 (brightness)
    data <- data[ , features, drop = FALSE]
    
    norm_data <- apply(data, 2, function(x) {
      if (max(x) == 0) {
        return(x)
      }
      norm_to <- max(x)/(100 - pt.brightness)
      return(x/norm_to)
    })
    
    # Convert to long
    ldata <- reshape(data.frame(norm_data), 
                     times = features, 
                     varying = features, 
                     timevar = "type", 
                     idvar = "row.names",
                     v.names= "l", 
                     direction = "long", 
                     new.row.names = seq(nrow(data) * length(features)))
    # Determine hues each type corresponding to
    if (is.null(colors)) {
      type_hue <- rep(90, length(features))
      angle_step <- 360/length(features)
      times <- seq(length(features)) - 1
      type_hue <- (type_hue + (angle_step * times)) %% 360
    } else {
      type_hue <- sapply(colors, function(x) rgb2hsv(col2rgb(x))["h", ] * 360)
    }
    names(type_hue) <- features
    ldata$h <- type_hue[ldata$type]
    ldata$s <- ifelse(ldata$l > 0, 100, 0)
    ldata$l <- ldata$l + pt.brightness
    
    
    
    # Convert hsl to rgb
    ldata$color <- apply(ldata[ , c("h", "s", "l")], 1, function(x) {
      hvalue <- min(x[1]/360, 1)
      svalue <- min(x[2]/100, 1)
      lvalue <- min(x[3]/100, 1) 
      hsv(hvalue, svalue, lvalue)
    })            
    ldata$r <- strtoi(substr(ldata$color, 2, 3), base = 16)
    ldata$g <- strtoi(substr(ldata$color, 4, 5), base = 16)
    ldata$b <- strtoi(substr(ldata$color, 6, 7), base = 16)
    
    min.col <- hsv(1, 0, pt.brightness/100)
    ldata$r[ldata$color == min.col] <- ldata$r[ldata$color == min.col]/2
    ldata$g[ldata$color == min.col] <- ldata$g[ldata$color == min.col]/2
    ldata$b[ldata$color == min.col] <- ldata$b[ldata$color == min.col]/2
    
    
    palatte <- data.frame(
      "r" = sapply(tapply(ldata$r, ldata$row.names, sum), 
                   function(x) min(x, 255)), 
      "g" = sapply(tapply(ldata$g, ldata$row.names, sum), 
                   function(x) min(x, 255)), 
      "b" = sapply(tapply(ldata$b, ldata$row.names, sum), 
                   function(x) min(x, 255))
    )
    finalhex <- apply(palatte, 1, function(x) {rgb(x[1], x[2], x[3], 
                                                   maxColorValue = 255)})
    
    return(finalhex)
  }
  
  make_legend <- function(data, features, colors = NULL, min.col) {
    dummy <- data[ , features, drop = FALSE]
    dummy[nrow(dummy) + 1, ] <- 0
    dummy$x <- c(1:nrow(dummy))
    dummy$y <- c(1:nrow(dummy))
    
    if (is.null(colors)) {
      type_hue <- rep(90, length(features))
      angle_step <- 360/length(features)
      times <- seq(length(features)) - 1
      type_hue <- (type_hue + (angle_step * times)) %% 360
    } else {
      type_hue <- sapply(colors, function(x) rgb2hsv(col2rgb(x))["h", ] * 360)
    }
    
    names(type_hue) <- features
    type_col <- sapply(type_hue, function(x) hsv(x/360, s = 1, v = 1))
    legends <- list("ncol" = ceiling(length(features)/4))
    for (i in features) {
      gene <- sym(i)
      legends[[i]] <- get_legend(
        ggplot(dummy, aes(x = x, y = y, color = !! gene)) + 
          geom_point() +
          scale_color_gradient(low = min.col, 
                               high = type_col[i]) +
          labs(color = featurekey_ori[i])
      )
    }
    combined_legend <- do.call(plot_grid, args = legends)
    return(combined_legend)
  }
  
  # Deal with cutoffs
  cut_values <- function(data, min.cutoff, max.cutoff) {
    data_mod <- apply(data, 2, function(x) {
      # Initialize each column as a vector
      values <- x
      
      # Initialize upper and lower limits to infinity
      lcutoff <- -Inf
      hcutoff <- Inf
      
      # lower cut
      if (!is.null(min.cutoff) ) {
        if (substr(min.cutoff, 1, 1) == "q") {
          lcutoff <- quantile(x, 
                              probs = as.numeric(substr(min.cutoff, 2, 
                                                        nchar(min.cutoff)))/100)
          values <- ifelse(values < lcutoff, 0, values)
        } else {
          lcutoff <- min.cutoff
        }
      } 
      
      # Upper cut
      if (!is.null(max.cutoff)) {
        if (substr(max.cutoff, 1, 1) == "q") {
          hcutoff <- quantile(x, 
                              probs = as.numeric(substr(max.cutoff, 2, 
                                                        nchar(max.cutoff)))/100)
        } else {
          hcutoff <- max.cutoff
        }
      }
      if (hcutoff < lcutoff) {
        stop("min.cutoff cannot exceed max.cutoff")
      }
      
      # Cap and floor the values
      values <- ifelse(values > hcutoff, hcutoff, values)
      values <- ifelse(values < lcutoff, 0, values)
      return(values)
    })
    data_mod <- data.frame(data_mod)
    return(data_mod)
  }
  
  ########################################
  if (is.null(cells.use)) {
    cells.use <- row.names(slot(object, name = "meta.data"))
  }
  
  if (!is.null(colors) & length(colors) != length(features)) {
    stop(paste(
      length(features), 
      "colors are required, but", 
      length(colors), 
      "are given."
    ))
  }
  
  if (is.null(assay)) {
    non_meta <- features[!features %in% colnames(slot(object, name = "meta.data"))]
    if (!all(grepl("_", non_meta))) {
      warning(paste("You seem not to be providing key-specified feature names.", 
                    "This could cause troubles when there are same feature names in multiple assays", 
                    "Use [key_feature] name to be safe (e.g., 'rna_dpn' instead of 'dpn' alone)."))
    }
  }
  
  # Extract cell embedings
  reduction <- tolower(reduction)
  reductionobj <- slot(object, name = "reductions")[[reduction]]
  cellemb <- data.frame(slot(reductionobj, name = "cell.embeddings"))
  
  # Get keys for the assay
  if (!is.null(assay)) {
    assayobj <- slot(object, name = "assays")[[assay]]
    assaykey <- slot(assayobj, name = "key")
    featurekey <- paste0(assaykey, features)
  } else {
    featurekey <- features
  }
  
  # Extract expression matrix
  featuresexp <- FetchData(object, vars = featurekey, 
                           cells = cells.use)
  
  # To deal with special characters in gene symbols
  featurekey_ori <- featurekey
  featurekey <- make.names(featurekey)
  names(featurekey_ori) <- featurekey
  colnames(featuresexp) <- make.names(colnames(featuresexp))
  
  # Deal with cutoffs
  featuresexp <- cut_values(featuresexp, min.cutoff, max.cutoff)
  featuresexp$color <- exp_encolor(featuresexp, features = featurekey, pt.brightness = pt.brightness)
  
  
  # Create plotting df
  plotdf <- merge(cellemb, featuresexp, by = "row.names")
  # Plotting
  axis_1 <- sym(colnames(cellemb[1]))
  axis_2 <- sym(colnames(cellemb[2]))
  
  base_plot <- ggplot(plotdf, aes(x = !!axis_1, y = !!axis_2)) +
    geom_point(color = plotdf$color, size = pt.size) +
    labs(x = paste(reduction, "1", sep = "_"), 
         y = paste(reduction, "2", sep = "_")) +
    theme_dark() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = background.color))
  
  if (label) {
    cellemb$ident <- slot(object, name = "active.ident")
    label_pos <- data.frame(
      x = tapply(cellemb[ ,1], cellemb$ident, mean), 
      y = tapply(cellemb[ ,2], cellemb$ident, mean), 
      label = names(tapply(cellemb[ ,2], cellemb$ident, mean)), 
      stringsAsFactors = FALSE
    )
    if (repel) {
      if (bold) {
        label_layer <- geom_text_repel(data = label_pos, aes(x = x, y = y, label = label), 
                                       fontface = "bold")
      } else {
        label_layer <- geom_text_repel(data = label_pos, aes(x = x, y = y, label = label))
      }
    } else {
      if (bold) {
        label_layer <- geom_text(data = label_pos, aes(x = x, y = y, label = label), 
                                 fontface = "bold")
      } else {
        label_layer <- geom_text(data = label_pos, aes(x = x, y = y, label = label))
      }
      
    }
    base_plot <- base_plot + label_layer
  }
  min.col <- hsv(1, 0, pt.brightness/100)
  fig_legends <- make_legend(featuresexp, featurekey, colors, min.col = min.col)
  
  result <- plot_grid(base_plot, fig_legends, 
                      rel_widths = c(fig.ratio, 1), ncol = 2)
  return(result)
}
