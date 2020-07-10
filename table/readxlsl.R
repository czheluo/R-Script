

install.packages("XLConnect")


setwd("C:\\Users\\meng.luo\\Desktop")
library("XLConnect")
library(xlsx)
#install.packages("xlsx")

wb <- loadWorkbook("result.xlsx")

sheet1 <- getSheets(wb)[[1]]

# get all rows
rows  <- getRows(sheet1)
cells <- getCells(rows)

styles <- sapply(cells, getCellStyle) #This will get the styles
cellColor <- function(style) {
  fg  <- style$getFillForegroundXSSFColor()
  rgb <- tryCatch(fg$getRgb(), error = function(e) NULL)
  rgb <- paste(rgb, collapse = "")
  return(rgb)
}

color_info <- sapply(styles, cellColor)


