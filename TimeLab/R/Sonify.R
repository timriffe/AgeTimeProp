# written 24 July, 2014 for fun

# first part copied from Fischer.R, sonifaction code after

# change this as needed:
setwd("/home/triffe/git/TimeProposal/TimeLab")

# -----------------------------------------------------
# Utility functions.
Mna0 <- cmpfun(function(M){
            M[is.na(M)]  <- 0
            M[is.nan(M)] <- 0
            M
        })

Minf0 <- cmpfun(function(M){
            M[is.infinite(M)]  <- 0
            M
        })

Thano <- cmpfun(function(Px, dx, stagger = TRUE){
            Np <- length(Px)
            Nd <- length(dx)
            if (Np != Nd){
                N <- max(Np, Nd)
                Px <- c(Px, rep(0, N - Np))
                dx <- c(dx, rep(0, N - Nd))
            } else {
                N <- Np
            }
            
            ay      <- 1:N - 1
            
            dx      <- Mna0(dx)   # remove NAs if any       
            dx      <- c(dx, dx * 0) / sum(dx) # pad out with 0s
            EDx     <- matrix(dx[col(matrix(nrow = N, 
                                            ncol = N)) + ay], 
                    nrow = N, 
                    ncol = N, 
                    dimnames = list(Ex = ay, 
                            Age = ay)
            )
            if (stagger){
                EDx <- (EDx + cbind(EDx[, 2:ncol(EDx)], 0)) / 2
            }
            t(Px * Minf0(Mna0(EDx / rowSums(EDx))))
        })

# this construction is rough, as things need to be staggered still,
# as with the Leslie matrix.
ThanoProjMatrix <- cmpfun(function(Fy, da, lambda = .9){
            N       <- length(Fy)
            # discount for part of infant mortality not surviving until end of year
            da[1]   <- da[1] * (1-lambda)
            
            # NxN matrix
            # fertility component
            Y       <- outer(da, Fy, "*") + 
                    # add survival element-wise
                    rbind(cbind(0,diag(N - 1)), 0)
            
            # reduce e0 fertility by 1/2, as only exposed for part of year
            Y[, 1]  <- Y[, 1] / 2
            # do not allow for Inf or NA values: impute 0
            Y       <- Mna0(Minf0(Y))
            # return projection matrix
            Y
        })

Leslie <- function(Fx, Lx){
    NN     <- length(Lx)
    N      <- NN - 1
    lambda <- -diff(Lx)/2
    Fx     <- (Fx[1:N] + Fx[-1])/2 * (1-lambda) 
    Sx     <- Lx[2:NN] / Lx[1:N]
    Mna0(cbind(rbind(Fx, diag(Sx)), 0))
}

# get data from web. Uncomment below code, and only run it once. Saves three big data.frames
# to a /Data/ folder in your working directory (you'll need to make that folder probably)
# -------------------------------------------------------------------------------------------#
# run only once, take a while due to HFD web structure and excessive parsing that is needed. #
# -------------------------------------------------------------------------------------------#
# HMD grabbing takes less time, but is slow due to size.                                     #
# -------------------------------------------------------------------------------------------#
# this code is way superfluous, as the example only uses data from 1950 US females, so skip  #
# this and adapt if you don' want the whole fricking databases. Commented out because you    #
# should only do it once!                                                                    #
# -------------------------------------------------------------------------------------------#
{
#library(devtools)
#install_github("DemogBerkeley", subdir = "DemogBerkeley", username = "UCBdemography")
#library(DemogBerkeley)
#library(reshape2)
#countries <- getHFDcountries()
#xyz <- countries[1]
#DAT <- do.call(rbind, lapply(countries,function(xyz,pw,us){ 
#            DAT    <- readHFDweb(xyz,"birthsRR",password=pw,username=us) 
#            BM     <- acast(DAT,Age~Year,value.var="Total")
#            YA     <- matrix(0,nrow=12,ncol=ncol(BM),dimnames=list(0:11,colnames(BM)))
#            OA     <- matrix(0,nrow=55,ncol=ncol(BM),dimnames=list(56:110,colnames(BM)))
#            BM     <- rbind(YA,BM,OA)
#            DAT    <- melt(BM, varnames=c("Age","Year"),value.name="Births")
#            DAT$Code <- xyz
#            DAT
#        },pw=pw,us=us))
    ## save out for later use
#save(DAT,file="Data/HFDall.Rdata")
    
# -------------------------------------------
# now get HMD exposures (differ from HFD exposures sometimes)
#library(DemogBerkeley)
#countries <- getHMDcountries()
#DAT <- do.call(rbind, lapply(countries,function(xyz,pw,us){
#            DAT             <- readHMDweb(xyz,"Exposures_1x1",password=pw,username=us) 
#            Fem             <- DAT[,c("Year","Age","Female")]
#            Mal             <- DAT[,c("Year","Age","Male")]
#            colnames(Fem)   <- colnames(Mal) <- c("Year","Age","Exp")
#            Fem$Sex         <- "f"
#            Mal$Sex         <- "m"
#            DAT             <- rbind(Fem,Mal)
#            DAT$Code        <- xyz
#            DAT
#        },pw=pw,us=us))
#save(DAT,file="Data/HMDExpall.Rdata")
# #and again for HMD lifetables (males grabbed superfluously)    
#DAT <- do.call(rbind, lapply(countries,function(xyz,pw,us){
#            DATm    <- readHMDweb(xyz,"mltper_1x1",password=pw,username=us) 
#            DATf    <- readHMDweb(xyz,"fltper_1x1",password=pw,username=us) 
#            DATm$Sex <- "m"
#            DATf$Sex <- "f"
#            DATi <- rbind(DATf,DATm)
#            DATi$Code <- xyz
#            DATi
#        },pw=pw,us=us))
#save(DAT,file="Data/HMDall.Rdata")
}

# assuming you have the above data.frames:
# Read in largeish data objects
HMDE    <- local(get(load("Data/HMDExpall.Rdata")))
HFD     <- local(get(load("Data/HFDall.Rdata")))
HMD     <- local(get(load("Data/HMDall.Rdata"))) # produced in Logo.R
# sort before binding
HMD     <- HMD[with(HMD, order(Code,Year,Sex,Age)),]
HMDE    <- HMDE[with(HMDE, order(Code,Year,Sex,Age)),]
HMD     <- cbind(HMD, Exp = HMDE$Exp)
# Deaths, not needed, but just to round things out
HMD$Deaths <- HMD$mx * HMD$Exp
# unisex, no need for males
HMD     <- HMD[HMD$Sex == "f", ]
# 'HMD' includes total births of both sexes.
# we can get crude sex ratios at birth
# from the HMD
# in a hurry, so assume 106 SRB
PF      <- 100 / 206  # PF is prop female, could have used K. Wachter's assumption too
HFD$Bf  <- HFD$Births * PF
# Code, Year:
.Code   <- "USA"
.Year   <- 1950
# --------------------------------------------------#
# grab data chunks based on year, country selection #
# --------------------------------------------------#
Bx <- with(HFD, Bf[Year == .Year & Code == .Code])
Ex <- with(HMD, Exp[Year == .Year & Code == .Code])
Fx <- Bx / Ex
Lx <- with(HMD, Lx[Year == .Year & Code == .Code]) / 1e5
# make Leslie Matrix
A  <- Leslie(Fx, Lx)
dx <- st(with(HMD, dx[Year == .Year & Code == .Code]))
By <- rowSums(Thano(Bx,dx))
Ey <- rowSums(Thano(Ex,dx))
Fy <- By / Ey
Y  <- ThanoProjMatrix(Fy,dx)
Av <- Re(eigen(A,symmetric=FALSE)$vectors)
Yv <- Re(eigen(Y,symmetric=FALSE)$vectors)
# ---------------------------------------------- #
# begin sonification chunk                       #
# ---------------------------------------------- #
# can't resist trying to sonify the eigenvectors...
# see playitbyr.org for better installation instructions.
# y synopsis:
# this requires installation of Csound. Google it to figure out how.
# First install Csound, THEN install these two R packages:

#install.packages("csound")
#install.packages("playitbyr")
library(csound)
library(playitbyr) 
library(reshape2)
#example(createPerformance) 
#createPerformance(sndcheck)

rownames(Yv) <- 0:110
colnames(Yv) <- 0:110
YvDF <- melt(Yv, varnames=c("y1","y2"),value.name="vec")
AvDF <- melt(Av, varnames=c("a1","a2"),value.name="vec")

# sonification code. Plenty that could be done to tweak it, but it requires
# either lots of trial and error or else someone who know their way around
# a synthesizer/Csound
Yson <- sonify(YvDF[YvDF$y2 < 10 & YvDF$y1 < 60,], sonaes(time = y1, pitch = vec)) + shape_scatter() + 
        scale_time_continuous(c(0, 2)) + scale_pitch_continuous(c(3, 12)) + sonfacet(y2)
Ason <- sonify(AvDF[AvDF$a2 < 10 & AvDF$a1 > 50,], sonaes(time = a1, pitch = vec)) + shape_scatter() + 
        scale_time_continuous(c(0, 2)) + scale_pitch_continuous(c(3, 12)) + sonfacet(a2)

# wave files on soundcloud, embedded in blog
sonsave(Yson, "Data/Yson.wav")
sonsave(Ason, "Data/Ason.wav")

# images in soundcloud:
png("Figures/Av.png")
matplot(0:110,Av,type='l',lty=1,col = "#00000020",ylim=c(-1,1))
dev.off()
png("Figures/Yv.png")
matplot(0:110,Yv,type='l',lty=1,col = "#00000020",ylim=c(-1,1))
dev.off()