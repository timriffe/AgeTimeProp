
setwd("/home/triffe/git/TimeProposal/TimeLab")
library(RColorBrewer)

Greens <- colorRampPalette(brewer.pal(9,"Greens")[1:7],space="Lab")
green <- Greens(7)[7]

makeFade <- function(x1, x2, ramp, rev = FALSE, n = 10){
    x <- seq(x1,x2,length=(n+1))
    cols <- ramp(n)
    if (rev){
        cols <- rev(cols)
    }
    rect(x[1:n],0,x[2:(n+1)],1,col=cols,border=NA,lwd=.3)
}
makeSolid <- function(x1,x2,col){
    rect(x1,0,x2,1,col=col,border=NA)
}
openPlot <- function(){
    par(mai=c(.05,.05,.05,.05),xpd=TRUE,xaxs="i",yaxs="i")
    plot(NULL, type = "n",axes=FALSE,xlab="",ylab="",xlim=c(0,5),ylim=c(0,1))
}
makeBorder <- function(x1,x2){
    rect(x1,0,x2,1,border=gray(.2),lwd=.5,lend='square')
}
graphics.off()

height <- .2
width  <- 3

pdf("Figures/Gantt/0_Years.pdf",width=width,height=height * 1.5)
openPlot()
text(.5:4.5,.5,1:5)
segments(0:5,0,0:5,1,lwd=.5,col=gray(.2))
segments(0,1,5,1,col=gray(.2))
dev.off()

#dev.new(width=5,height=.5)
pdf("Figures/Gantt/1_Math.pdf",width=width,height=height)
openPlot()
makeSolid(0,1.5,green)
makeFade(1.5,5,Greens,TRUE,21)
makeBorder(0,5)
dev.off()

pdf("Figures/Gantt/1_Methods.pdf",width=width,height=height)
openPlot()
makeSolid(0,2.5,green)
makeFade(2.5,5,Greens,TRUE,15)
makeBorder(0,5)
dev.off()

pdf("Figures/Gantt/2_Data.pdf",width=width,height=height)
openPlot()
makeSolid(1,2,green)
makeFade(0,1,Greens,FALSE,6)
makeFade(2,3,Greens,TRUE,6)
makeBorder(0,3)
dev.off()

pdf("Figures/Gantt/3_Empirical.pdf",width=width,height=height)
openPlot()
makeSolid(3,5,green)
makeFade(0,3,Greens,FALSE,18)
#makeFade(2,3,Greens,TRUE,12)
makeBorder(0,5)
dev.off()

