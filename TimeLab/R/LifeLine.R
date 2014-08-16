setwd("/home/triffe/git/TimeProposal/TimeLab")
# get plotting functions
source("R/Functions.R")

a <- 5

graphics.off()
#dev.new(height = 1, width = 6)
pdf("Figures/LifeLine.pdf", height = 1.25, width = 6)
par(mai = c(0, .5, .25, .5), xaxs = "i", yaxs = "i")
plot(NULL, 
        type = 'n', 
        xlim = c(0, 10), 
        ylim = c(0, 2), 
        axes = FALSE, 
        xlab = "", 
        ylab = "",
        asp = 1)
curlyBrace1(1, 1.5, length = a, radius1 = .25, radius2 = .15, col = gray(.5) )
curlyBrace1(a + 1, 1.5, length = 8 - a, radius1 = .25, radius2 = .15, col = gray(.5))
curlyBrace1(1, .5, length = 8, radius1 = .25, radius2 = .15, top = FALSE, col = gray(.5))

segments(1, 1, a + 1, 1, lwd = 2, col = "forestgreen")
segments(a + 1, 1, 9, 1, lwd = 2, col = "blue")
points(1, 1, pch = 1, cex = 1.5)
points(9, 1, pch = 19, cex = 1.5)

segments(a + 1, .9, a + 1, 1.1, lwd = 1)
text(c(1 + a / 2, 1 + a + (8 - a) / 2), 1.8, c("years lived", "years left"), col = c("forestgreen", "blue"))
text(5, .2, "lifespan")
text(0, 1, "Birth\n(state entry)", xpd = TRUE)
text(10, 1, "Death\n(state exit)", xpd = TRUE)
dev.off()



