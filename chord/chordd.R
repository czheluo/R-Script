chordd <- function (data, type = "directional", width = NULL, height = NULL, 
          margin = 100, palette = "Dark2", palette2 = "Greys", showGroupnames = TRUE, 
          groupNames = NULL, groupColors = NULL, groupThickness = 0.1, 
          groupPadding = 2, groupnamePadding = 30, groupnameFontsize = 18, 
          groupedgeColor = NULL, chordedgeColor = "#808080", categoryNames = NULL, 
          categorynamePadding = 100, categorynameFontsize = 28, showTicks = TRUE, 
          tickInterval = NULL, ticklabelFontsize = 10, fadeLevel = 0.1, 
          showTooltips = TRUE, showZeroTooltips = TRUE, tooltipNames = NULL, 
          tooltipUnit = NULL, tooltipFontsize = 12, tooltipGroupConnector = " &#x25B6; ", 
          precision = NULL, clickAction = NULL, clickGroupAction = NULL) 
{
  if (!is.matrix(data)) 
    stop("'data' must be a matrix class object.")
  d <- dim(data)
  if (type == "bipartite") {
    g1 <- d[1]
    g2 <- d[2]
    n <- g1 + g2
    m <- matrix(0, nrow = n, ncol = n)
    m[1:g1, (g1 + 1):n] <- data
    m[(g1 + 1):n, 1:g1] <- t(data)
    g1.names <- row.names(data)
    g2.names <- colnames(data)
    m.names <- c(g1.names, g2.names)
    row.names(m) <- m.names
    colnames(m) <- m.names
    if (is.null(categoryNames)) {
      categoryNames <- names(dimnames(data))
    }
    data <- m
  }
  else if (type == "directional") {
    if (d[1] != d[2]) 
      stop("'data' must be a square matrix.")
    n <- d[1]
  }
  if (!is.null(categoryNames) & type != "bipartite") {
    warning("category names are only used for bipartite chord diagrams.")
  }
  if (!is.null(groupNames)) {
    g <- length(groupNames)
    if (g != n) 
      stop(paste0("length of 'groupNames' [", g, "] not equal to matrix extent [", 
                  n, "]."))
  }
  else {
    groupNames <- colnames(data)
  }
  rnames <- row.names(data)
  if (!is.null(rnames)) {
    if (!identical(rnames, groupNames)) 
      warning("row names of the 'data' matrix differ from its column names or the 'groupNames' argument.")
  }
  if (length(groupnamePadding) == 1) {
    groupnamePadding <- rep(groupnamePadding, n)
  }
  if (is.null(groupColors)) {
    if (type == "directional") {
      color <- grDevices::colors()[grep("gr(a|e)y", grDevices::colors(), invert = T)]
      rcolor <- color[sample(1:length(color))]
      groupColors <- rcolor#RColorBrewer::brewer.pal(n, palette)
    }
    else if (type == "bipartite") {
      color <- grDevices::colors()[grep("gr(a|e)y", grDevices::colors(), invert = T)]
      rcolor <- color[sample(1:length(color))]
      expcolor<-grDevices::rainbow(100)
      groupColors <- paste(color[c(1:g1)],expcolor[c(1:g2)])
      print(groupColors)
      print(g1)
      print(g2)
      #c(RColorBrewer::brewer.pal(g1,palette), 
       # RColorBrewer::brewer.pal(g2,palette2))
                       
    }
  }
  if (is.null(tickInterval)) {
    tickInterval <- 10^(floor(log10(max(data))) - 1)
  }
  if (is.null(tooltipNames)) {
    tooltipNames = groupNames
  }
  if (is.null(tooltipUnit)) {
    tooltipUnit <- ""
  }
  if (is.null(precision)) {
    precision <- "null"
  }
  params = list(matrix = data, options = list(type = type, 
                                              width = width, height = height, margin = margin, showGroupnames = showGroupnames, 
                                              groupNames = groupNames, groupColors = groupColors, 
                                              groupThickness = groupThickness, groupPadding = pi * 
                                                groupPadding/180, groupnamePadding = groupnamePadding, 
                                              groupnameFontsize = groupnameFontsize, groupedgeColor = groupedgeColor, 
                                              chordedgeColor = chordedgeColor, categoryNames = categoryNames, 
                                              categorynamePadding = categorynamePadding, categorynameFontsize = categorynameFontsize, 
                                              showTicks = showTicks, tickInterval = tickInterval, 
                                              ticklabelFontsize = ticklabelFontsize, fadeLevel = fadeLevel, 
                                              showTooltips = showTooltips, showZeroTooltips = showZeroTooltips, 
                                              tooltipNames = tooltipNames, tooltipFontsize = tooltipFontsize, 
                                              tooltipUnit = tooltipUnit, tooltipGroupConnector = tooltipGroupConnector, 
                                              precision = precision, clickAction = clickAction, clickGroupAction = clickGroupAction))
  params = Filter(Negate(is.null), params)
  htmlwidgets::createWidget(name = "chorddiag", params, width = width, 
                            height = height, htmlwidgets::sizingPolicy(padding = 0, 
                                                                       browser.fill = TRUE), package = "chorddiag")
}
