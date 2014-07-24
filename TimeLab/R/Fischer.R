# written 23 July, 2014 in order to verify stuff I write before I write it!

setwd("/home/triffe/git/TimeProposal/TimeLab")
source("R/Functions.R")
# -------------------------------------------------------------------------------------------#
# run only once, take a while due to HFD web structure and excessive parsing that is needed. #
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
}

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
PF <- 100 / 206  # PF is prop female, could have used K. Wachter's assumption too
HFD$Bf <- HFD$Births * PF

# Code, Year:

.Code <- "USA"
.Year <- 1950

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
# close, but not same because Fy didn't come from stable pop
# had it come from stable pop, it'd still not be same due
# to discrete approx and rounding, but it'd be closer still
log(max(Re(eigen(A)$values)))
log(max(Re(eigen(Y)$values)))

# eigenvectors
#library(Matrix)
#A <- Matrix(A, sparse = TRUE)
Av <- Re(eigen(A,symmetric=FALSE)$vectors)
Yv <- Re(eigen(Y,symmetric=FALSE)$vectors)

graphics.off()
png("Figures/Eigens.png",width=900,height=450)
par(mfrow=c(1,2), mai = c(.5,.5,.5,.5))
matplot(0:110,Av,type='l',lty=1,col = "#00000020",ylim=c(-1,1), main = "chronological age")
lines(0:110,abs(Av[,1]),col="red")
matplot(0:110,Yv,type='l',lty=1,col = "#00000020",ylim=c(-1,1), main = "thanatological age")
lines(0:110,abs(Yv[,1]),col="red")
dev.off()
EA <- eigen.analysis(A)
EY <- eigen.analysis(Y)

# this is necessary because solve(eigen(A)$vectors)) returns singularity error
# eigen.analysis frequently won't work for these particular Leslie matrices
getRv <- function(A,leslie=TRUE){
    v1 <- eigen( t(A) )$vectors[,1]
    if (leslie){
        v1 <- Re(v1)/Re(v1[1])
    }
    v1
}
divst <- function(x){
    x/sd(x)
}
# the usual normalization
#plot(0:110,divst(getRv(A)),type="l", ylim=c(0,3))
#lines(0:110,divst(getRv(Y,FALSE)), col = "green")
#
## stuff to add, out of curiosity
#abline(v = sum(.5:110.5*dx),col="red")
#abline(v=which.max(dx)-.5,col="red",lty=2)
#lines(0:110,divst(Fx),col = "blue")
#lines(0:110,divst(dx),col = "magenta")

