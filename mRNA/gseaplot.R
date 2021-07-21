library(RColorBrewer)
library(ggplot2)
library(enrichplot)
pp <- ggplot(gse, aes(x = Rank_at_MAX)) + xlab(NULL) + theme_classic(base_size) +
        theme(panel.grid.major = element_line(colour = "grey92"),
            panel.grid.minor = element_line(colour = "grey92"),
            panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
        scale_x_continuous(expand = c(0, 0))



esline <- geom_line(aes(y = ES, color = "green"),size = 1)


pes <- pp + esline + theme(legend.position = c(0.8, 0.8),
        legend.title = element_blank(), legend.background = element_rect(fill = "transparent"))
pes <- pes + ylab("Running Enrichment Score") + theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.line.x = element_blank(),
        plot.margin = margin(t = 0.2, r = 0.2, b = 0, l = 0.2,
            unit = "cm"))



gsInfo <- function(object, geneSetID) {
    geneList <- object@geneList

    if (is.numeric(geneSetID))
        geneSetID <- object@result[geneSetID, "ID"]

    geneSet <- object@geneSets[[geneSetID]]
    exponent <- object@params[["exponent"]]
    df <- gseaScores(geneList, geneSet, exponent, fortify=TRUE)
    df$ymin <- 0
    df$ymax <- 0
    pos <- df$position == 1
    h <- diff(range(df$runningScore))/20
    df$ymin[pos] <- -h
    df$ymax[pos] <- h
    df$geneList <- geneList

    df$Description <- object@result[geneSetID, "Description"]
    return(df)
}

 h <- diff(range(gse$ES))/20
 gse$ymin[pos] <- -h
 gse$ymax[pos] <- h

gse$ymin <- 0
gse$ymax <- 1
gse$position <- 1
gse$x<-seq(1:length(gse$ES))
p2 <- ggplot(gse, aes(x = x)) + geom_linerange(aes(ymin = ymin,
        ymax = ymax, color = "red")) + xlab(NULL) + ylab(NULL) +
        theme_classic(base_size) + theme(legend.position = "none",
        plot.margin = margin(t = -0.1, b = 0, unit = "cm"), axis.ticks = element_blank(),
        axis.text = element_blank(), axis.line.x = element_blank()) +
        scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0,
        0))
v <- seq(1, sum(gse$Rank_at_MAX), length.out = 9)
inv <- findInterval(rev(cumsum(gse$position)), v)
if (min(inv) == 0)
    inv <- inv + 1
col = c(rev(brewer.pal(5, "Blues")), brewer.pal(5, "Reds"))
ymin <- min(p2$data$ymin)
yy <- max(p2$data$ymax - p2$data$ymin) * 0.3
xmin <- which(!duplicated(inv))
xmax <- xmin + as.numeric(table(inv)[unique(inv)])
d <- data.frame(ymin = ymin, ymax = yy, xmin = xmin,
            xmax = xmax, col = col[unique(inv)])
p2 <- p2 + geom_rect(aes_(xmin = ~xmin, xmax = ~xmax,
          ymin = ~ymin, ymax = ~ymax, fill = ~I(col)), data = d,
          alpha = 0.9, inherit.aes = FALSE)

df2 <- p$data
df2$y <- p$data$geneList[df2$x]
p.pos <- p + geom_segment(data = df2, aes_(x = ~x, xend = ~x,
y = ~y, yend = 0), color = "grey")
p.pos <- p.pos + ylab("Ranked list metric") + xlab("Rank in Ordered Dataset") +
        theme(plot.margin = margin(t = -0.1, r = 0.2, b = 0.2,
            l = 0.2, unit = "cm"))

plotlist <- list(pes, p2)
n <- length(plotlist)
plotlist[[n]] <- plotlist[[n]] + theme(axis.line.x = element_line(),
        axis.ticks.x = element_line(), axis.text.x = element_text())

    plot_grid(plotlist = plotlist, ncol = 1, align = "v", rel_heights = rel_heights)    
 
