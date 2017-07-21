
my.facto_summarize <- function (X, element, node.level = 1, group.names, result = c("coord", 
    "cos2", "contrib"), axes = 1:2, select = NULL) 
{
    if (!element %in% c("row", "col", "var", "ind", "quanti.var", 
        "quali.var", "group", "partial.axes", "partial.node")) 
        stop("Te argument element should be one of \"row\", \"col\", \"var\", \"ind\", \"quanti.var\", \"quali.var\", \"group\", \"partial.axes\", \"partial.node\"")
    facto_class <- factoextra:::.get_facto_class(X)
    element <- element[1]
    if (facto_class == "CA") {
        if (element %in% c("ind", "row")) 
            elmt <- get_ca_row(X)
        else if (element %in% c("var", "col")) 
            elmt <- get_ca_col(X)
    }
    else if (facto_class == "PCA") {
        if (element %in% c("var", "col")) 
            elmt <- get_pca_var(X)
        else if (element %in% c("ind", "row")) 
            elmt <- get_pca_ind(X)
    }
    else if (facto_class == "MCA") {
        if (element %in% c("var", "col")) 
            elmt <- get_mca_var(X)
        else if (element %in% c("ind", "row")) 
            elmt <- get_mca_ind(X)
    }
    else if (facto_class == "MFA") {
        if (element %in% c("quanti.var", "col")) 
            elmt <- get_mfa_quanti_var(X)
        else if (element %in% c("quali.var", "col")) 
            elmt <- get_mfa_quali_var(X)
        else if (element %in% c("group", "col")) 
            elmt <- get_mfa_group(X)
        else if (element %in% c("partial.axes", "col")) 
            elmt <- get_mfa_partial_axes(X)
        else if (element %in% c("ind", "row")) 
            elmt <- get_mfa_ind(X)
    }
    else if (facto_class == "HMFA") {
        if (element %in% c("quanti.var", "col")) 
            elmt <- get_hmfa_quanti_var(X)
        else if (element %in% c("quali.var", "col")) 
            elmt <- get_hmfa_quali_var(X)
        else if (element %in% c("group", "col")) 
            elmt <- get_hmfa_group(X)
        else if (element %in% c("ind", "row")) 
            elmt <- get_hmfa_ind(X)
        else if (element %in% c("partial.node", "row")) 
            elmt <- get_hmfa_partial(X)
    }
    if (class(elmt)[[2]] == "hmfa_partial") 
        ndim <- ncol(elmt[[1]])
    else ndim <- ncol(elmt$coord)
    if (max(axes) > ndim) 
        stop("The value of the argument axes is incorrect. ", 
            "The number of axes in the data is: ", ncol(elmt$coord), 
            ". Please try again with axes between 1 - ", ncol(elmt$coord))
    res = NULL
    if ("coord" %in% result) {
        dd <- data.frame(elmt$coord[, axes, drop = FALSE])
        coord <- apply(dd^2, 1, sum)
        res = cbind(dd, coord = coord)
    }
    if ("cos2" %in% result) {
        cos2 <- elmt$cos2[, axes]
        if (length(axes) > 1) 
            cos2 <- apply(cos2, 1, sum, na.rm = TRUE)
        res <- cbind(res, cos2 = cos2)
    }
    if ("contrib" %in% result) {
        contrib <- elmt$contrib[, axes]
        if (length(axes) > 1) {
            eig <- get_eigenvalue(X)[axes, 1]
            contrib <- t(apply(contrib, 1, function(var.contrib, 
                pc.eig) {
                #var.contrib * pc.eig
                var.contrib * pc.eig/sum(pc.eig)
            }, eig))
            contrib <- apply(contrib, 1, sum)
        }
        res <- cbind(res, contrib = contrib)
    }
    if ("coord.partial" %in% result) {
        dd <- data.frame(elmt$coord.partiel[, axes, drop = FALSE])
        groupnames <- data.frame(do.call("rbind", strsplit(as.character(rownames(dd)), 
            ".", fixed = TRUE)))
        colnames(groupnames) <- c("name", "group.name")
        coord.partial <- apply(dd^2, 1, sum)
        res.partial <- data.frame(groupnames, dd, coord.partial)
    }
    if ("coord.node.partial" %in% result) {
        node <- as.data.frame(elmt[[node.level]])
        name <- rep(rownames(node), length(group.names))
        dim.group <- dim.names <- dd <- coord.partial <- dim.coord <- dim.name <- NULL
        for (i in axes[1]:length(axes)) {
            dim.group <- NULL
            for (j in 1:length(group.names)) {
                dim.name <- paste0("Dim", axes[i], ".", j)
                dim.coord <- abind::abind(dim.coord, node[, dim.name])
                dim.group <- c(dim.group, rep(group.names[j], 
                  length(node[, dim.name])))
            }
            dim.names <- c(dim.names, paste0("Dim", axes[i]))
            dd <- cbind(dd, dim.coord)
            dim.coord <- NULL
        }
        colnames(dd) <- dim.names
        coord.partial <- apply(dd^2, 1, sum)
        res.partial <- data.frame(group.name = dim.group, name, 
            dd, coord.partial)
    }
    if ("coord.node.partial" %in% result) 
        res <- res.partial
    else {
        name <- rownames(elmt$coord)
        if (is.null(name)) 
            name <- as.character(1:nrow(elmt$coord))
        res <- cbind.data.frame(name = name, res)
        if (!is.null(select)) 
            res <- factoextra:::.select(res, select)
        if ("coord.partial" %in% result) {
            res = list(res = res, res.partial = res.partial)
        }
    }
    res
}

my.fviz_contrib <-  function (X, choice = c("row", "col", "var", "ind", "quanti.var", 
    "quali.var", "group", "partial.axes"), axes = 1, fill = "steelblue", 
    color = "steelblue", sort.val = c("desc", "asc", "none"), 
    top = Inf) 
{
    #title <- .build_title(choice[1], "Contribution", axes)
    title <- factoextra:::.build_title(choice[1], "Contribution", axes)
    dd <- my.facto_summarize(X, element = choice, result = "contrib", 
        axes = axes)
    contrib <- dd$contrib
    names(contrib) <- rownames(dd)
    theo_contrib <- 100/length(contrib)
    if (length(axes) > 1) {
        eig <- get_eigenvalue(X)[axes, 1]
        eig <- eig/sum(eig)
        theo_contrib <- sum(theo_contrib * eig)
    }
    p <- factoextra:::.ggbarplot(contrib, fill = fill, color = color, sort.value = sort.val[1], 
        top = top, title = title, ylab = "Contributions (%)") + 
        geom_hline(yintercept = theo_contrib, linetype = 2, color = "red")
    p
}
