setwd("C:/Users/arctu/Dropbox/School/Thesis/qPCR Data/CFX96 Data")

data1 <- read.csv(file="Fold Change Full.csv")

panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...){                    ## Define a new function that codes for scatterplot matrix with r,
  usr <- par("usr"); on.exit(par(usr))                                    ## extra information added with code by Hossam Abdel-Moniem.
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

pairs(data1, lower.panel = panel.smooth, upper.panel = panel.cor)
