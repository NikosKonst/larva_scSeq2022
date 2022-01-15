library(ranger)

SplitRename <- function(x, ident1, ident2) {
    idents <- x
    idents <- as.character(idents)
    idents[!idents %in% c(ident1, ident2)] <- "others"
    
    idents[idents %in% ident1] <- "1"
    idents[idents %in% ident2] <- "2"
    
    idents <- factor(idents)
    return(idents)
    
}

SplitByIdent <- function(object, ident1, ident2, slot.use = "data", assays = "RNA", 
						ratio = 0.1, min_cells = 5, split.test = FALSE, 
						keep.orig.ident = FALSE, features = NULL) {
    clade_to_test <- subset(object, idents = c(ident1, ident2))
	data_obj <- slot(clade_to_test, name = "assays")[[assays]]
    if(is.null(features)){
		features = data_obj@var.features
	}

    if (!split.test) {
            results <- list(exp = t(as.matrix(slot(data_obj, slot.use)[features, ])),
                            ident = clade_to_test@active.ident)
            if (keep.orig.ident) {
                results[["orig.ident"]] <- results[["ident"]]
            }

    # Recode idents
    results[["ident"]] <- SplitRename(x = results[["ident"]],
                                      ident1 = ident1, ident2 = ident2)

    return(results)
    }

    # Decide the size of test set (at least five cells)
    freq_tbl <- as.data.frame(table(clade_to_test@active.ident))
    freq_tbl$test_num <- sapply(round(freq_tbl$Freq * ratio), function(x) return(max(x, min_cells)))
    row.names(freq_tbl) <- freq_tbl$Var1

    # Take the cell list for each ident
    id_table <- as.data.frame(clade_to_test@active.ident)
    colnames(id_table) <- "id"
    cbclist <- lapply(c(ident1, ident2), function(x) {
        row.names(id_table)[id_table$id == x]
    })
    names(cbclist) <- c(ident1, ident2)

    # Define the test set
    test_cbc <- unlist(lapply(names(cbclist), function(x) {
        number <- freq_tbl[x, "test_num"]
        cell_bc <- sample(cbclist[[as.character(x)]], size = number, replace = FALSE)
        return(cell_bc)

    }))

    # Define train cbc
    train_cbc <- setdiff(row.names(id_table), test_cbc)



    # Get expression matrices for training and testing
	print(train_cbc[1:5])
    train_mat <- slot(data_obj, slot.use)[features, train_cbc]
    test_mat <- slot(data_obj, slot.use)[features, test_cbc]

    # Define a list to return
    results <- list(train = t(as.matrix(train_mat)), test = t(as.matrix(test_mat)),
                    train_id = id_table[train_cbc, ], test_id = id_table[test_cbc, ])
    if (keep.orig.ident) {
        results[["train_id_orig"]] <- results[["train_id"]]
        results[["test_id_orig"]] <- results[["test_id"]]
    }
    results[["train_id"]] <- SplitRename(x = results[["train_id"]],
                                  ident1 = ident1, ident2 = ident2)
    results[["test_id"]] <- SplitRename(x = results[["test_id"]],
                              ident1 = ident1, ident2 = ident2)
    return(results)
}

GetClusterBranch <- function(object, node) {
    edgelist <- object@tools$BuildClusterTree$edge

# Correction to address that node/tip number are 1-based while Seurat cluster identities are 0-based
    tiplist <- seq(length(object@tools$BuildClusterTree$tip.label))
    tipname <- object@tools$BuildClusterTree$tip.label
    names(tipname) <- tiplist

    init_branches <- edgelist[edgelist[, 1] %in% node , 2]

    metaclusters <- lapply(init_branches, function(x) {
        nodes <- x
        i <- 1
        while(any(nodes %in% edgelist[ , 1])) {
            nodes_new <- edgelist[edgelist[, 1] %in% nodes, 2]
            nodes <- c(intersect(nodes, tiplist), nodes_new)
        }
        nodes <- tipname[nodes]
        return(nodes)
    })
}


AssessNode <- function(object, node, assay) {
    traindata <- SplitByIdent(object, assays = assay,
                              ident1 = GetClusterBranch(object, node = node)[[1]],
                              ident2 = GetClusterBranch(object, node = node)[[2]])

    to_train <- as.data.frame(traindata$exp)
    to_train$ident <- traindata$ident

    rfmodel <- ranger(data = to_train, dependent.variable.name = "ident",
                  classification = TRUE, num.trees = 1500)
    oob <- rfmodel$prediction.error
    return(oob)
}



GetInternalNodes <- function(object) {
    nodelist <- slot(object, name = "tools")[["BuildClusterTree"]][["edge"]]
    terminal_nodes <- setdiff(nodelist[ , 2], nodelist[ , 1])
    internal_nodes <- setdiff(as.vector(unique(nodelist)), terminal_nodes)
    return(internal_nodes)
}



AssessAllNodes <- function(object, cut.off = 0.05, assay = "RNA") {
    nodelist <- GetInternalNodes(object)
    oob_error <- lapply(nodelist, function(x) {
        return(AssessNode(object, x, assay))
    })
    names(oob_error) <- nodelist
    results <- list(
        oobe = oob_error,
        low_confidence = names(oob_error)[oob_error > cut.off]
    )
    return(results)
}

IdentRename <- function(x, old.ident, new.ident) {
    idents <- x
    idents <- as.character(idents)
    idents[idents %in% old.ident] <- new.ident
    return(idents)
    
}

AssessPair <- function(object, ident1, ident2, assay = "RNA") {
    traindata <- SplitByIdent(object, assays = assay,
                              ident1 = ident1,
                              ident2 = ident2)

    to_train <- as.data.frame(traindata$exp)
    to_train$ident <- traindata$ident

    rfmodel <- ranger(data = to_train, dependent.variable.name = "ident",
                  classification = TRUE, num.trees = 1500)
    oob <- rfmodel$prediction.error
    return(oob)
}
