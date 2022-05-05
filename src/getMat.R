getmat <- function (data, nsets = 6, nintersects = 40, sets = NULL, keep.order = F, 
    set.metadata = NULL, intersections = NULL, matrix.color = "gray23", mat_col=NULL,
    main.bar.color = "gray23", mainbar.y.label = "Intersection Size", 
    mainbar.y.max = NULL, sets.bar.color = "gray23", sets.x.label = "Set Size", 
    point.size = 2.2, line.size = 0.7, mb.ratio = c(0.7, 0.3), 
    expression = NULL, att.pos = NULL, att.color = main.bar.color, 
    order.by = c("freq", "degree"), decreasing = c(T, F), show.numbers = "yes", 
    number.angles = 0, group.by = "degree", cutoff = NULL, queries = NULL, 
    query.legend = "none", shade.color = "gray88", shade.alpha = 0.25, 
    matrix.dot.alpha = 0.5, empty.intersections = NULL, color.pal = 1, 
    boxplot.summary = NULL, attribute.plots = NULL, scale.intersections = "identity", 
    scale.sets = "identity", text.scale = 1, set_size.angles = 0, 
    set_size.show = FALSE, set_size.numbers_size = NULL, set_size.scale_max = NULL)  {

    startend <- UpSetR:::FindStartEnd(data)
    first.col <- startend[1]
    last.col <- startend[2]
    if (color.pal == 1) {
        palette <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", 
            "#9467BD", "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", 
            "#17BECF")
    }
    else {
        palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
            "#0072B2", "#D55E00", "#CC79A7")
    }
    if (is.null(intersections) == F) {
        Set_names <- unique((unlist(intersections)))
        Sets_to_remove <- UpSetR:::Remove(data, first.col, last.col, Set_names)
        New_data <- UpSetR:::Wanted(data, Sets_to_remove)
        Num_of_set <-UpSetR:::Number_of_sets(Set_names)
        if (keep.order == F) {
            Set_names <- UpSetR:::order_sets(New_data, Set_names)
        }
        All_Freqs <- UpSetR:::specific_intersections(data, first.col, 
            last.col, intersections, order.by, group.by, decreasing, 
            cutoff, main.bar.color, Set_names)
    }
    else if (is.null(intersections) == T) {
        Set_names <- sets
        if (is.null(Set_names) == T || length(Set_names) == 0) {
            Set_names <- UpSetR:::FindMostFreq(data, first.col, last.col, 
                nsets)
        }
        Sets_to_remove <- UpSetR:::Remove(data, first.col, last.col, Set_names)
        New_data <- UpSetR:::Wanted(data, Sets_to_remove)
        Num_of_set <- UpSetR:::Number_of_sets(Set_names)
        if (keep.order == F) {
            Set_names <- UpSetR:::order_sets(New_data, Set_names)
        }
        All_Freqs <- UpSetR:::Counter(New_data, Num_of_set, first.col, 
            Set_names, nintersects, main.bar.color, order.by, 
            group.by, cutoff, empty.intersections, decreasing)
    }
    Matrix_setup <- UpSetR:::Create_matrix(All_Freqs)
    labels <- UpSetR:::Make_labels(Matrix_setup)
    att.x <- c()
    att.y <- c()
    if (is.null(attribute.plots) == F) {
        for (i in seq_along(attribute.plots$plots)) {
            if (length(attribute.plots$plots[[i]]$x) != 0) {
                att.x[i] <- attribute.plots$plots[[i]]$x
            }
            else if (length(attribute.plots$plots[[i]]$x) == 
                0) {
                att.x[i] <- NA
            }
            if (length(attribute.plots$plots[[i]]$y) != 0) {
                att.y[i] <- attribute.plots$plots[[i]]$y
            }
            else if (length(attribute.plots$plots[[i]]$y) == 
                0) {
                att.y[i] <- NA
            }
        }
    }
    BoxPlots <- NULL
    if (is.null(boxplot.summary) == F) {
        BoxData <- UpSetR:::IntersectionBoxPlot(All_Freqs, New_data, first.col, 
            Set_names)
        BoxPlots <- list()
        for (i in seq_along(boxplot.summary)) {
            BoxPlots[[i]] <- UpSetR:::BoxPlotsPlot(BoxData, boxplot.summary[i], 
                att.color)
        }
    }
    customAttDat <- NULL
    customQBar <- NULL
    Intersection <- NULL
    Element <- NULL
    legend <- NULL
    EBar_data <- NULL
    if (is.null(queries) == F) {
        custom.queries <- UpSetR:::SeperateQueries(queries, 2, palette)
        customDat <- UpSetR:::customQueries(New_data, custom.queries, 
            Set_names)
        legend <- UpSetR:::GuideGenerator(queries, palette)
        legend <- UpSetR:::Make_legend(legend)
        if (is.null(att.x) == F && is.null(customDat) == F) {
            customAttDat <- UpSetR:::CustomAttData(customDat, Set_names)
        }
        customQBar <- UpSetR:::customQueriesBar(customDat, Set_names, 
            All_Freqs, custom.queries)
    }
    if (is.null(queries) == F) {
        Intersection <- UpSetR:::SeperateQueries(queries, 1, palette)
        Matrix_col <- intersects(UpSetR:::QuerieInterData, Intersection, 
            New_data, first.col, Num_of_set, All_Freqs, expression, 
            Set_names, palette)
        Element <- UpSetR:::SeperateQueries(queries, 1, palette)
        EBar_data <- UpSetR:::ElemBarDat(Element, New_data, first.col, 
            expression, Set_names, palette, All_Freqs)
    }
    else {
        Matrix_col <- NULL
    }
    if (!is.null(mat_col)) {
      Matrix_col <- mat_col
    }
    Matrix_layout <- UpSetR:::Create_layout(Matrix_setup, matrix.color, 
        Matrix_col, matrix.dot.alpha)
    return(Matrix_layout)

}
