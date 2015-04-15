setwd("/home/tim/git/AgeTimeProp/TimeLab")

# this object is produced by R/DataPrep.R ...
DAT     <- local(get(load("Data/HMDall.Rdata")))
library(reshape2)

lx      <- acast(DAT, Age ~ Year + Sex + Code, value.var = "lx")
lx      <- lx / 1e5      # divide out radix
n       <- 1000          # haphazard selection worked ok!
colnames(lx)[n]          # 1875, females, FRATNP 
lxp     <- lx[56:111, n] # 2nd 1/2 of lx makes a good hourglass building block!
lxp     <- c(rev(lxp), lxp[ -1]) # symmetric

bluecol <- colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"), space = "Lab")(9)[3]
ylcol   <- colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"), space = "Lab")(9)[3]

graphics.off()
border  <- gray(.3)
mai     <- rep(.03, 4) # narrow borders
height  <- 2 + mai[1] * 2
width   <- 1 + mai[1] * 2
#dev.new(height=height,width=width)
pdf("Figures/logo.pdf", height = height, width = width)
par(xaxs = "i", yaxs = "i", mai = mai, xpd = TRUE)
x <- seq(0, 2, length = 111)
# open empty device 
plot(NULL, type = 'n', ylim = c(0, 2), xlim = c(0, 1), asp = 1, axes = FALSE, xlab = "", ylab = "")
# colored triangles (underneath)
polygon(c(0, 1, 1), c(0, 0, 1), col = bluecol, border = NA)
polygon(c(0, 1, 0), c(1, 2, 2), col = ylcol, border = NA)
# hourglass
polygon(c(lxp, 1 - lxp), c(x, rev(x)), border = border, lwd = 2, col = gray(.9))
# triangle outlines
polygon(c(0, 1, 1), c(0, 0, 1), border = border, lwd = 2)
polygon(c(0, 1, 0), c(1, 2, 2), border = border, lwd = 2)
dev.off()


system('cd /home/tim/git/AgeTimeProp/TimeLab/CoverLetter | pdflatex Riffe_Cover.tex')