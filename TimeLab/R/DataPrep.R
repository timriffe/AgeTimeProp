# all data is read live from the web using the DemogBerkeley package.
# this code takes a long time to run and should only be run once!
# Note all code relying on these data actually only uses a single example
# year/sex of data, so you can actually save some time by adapting the code
# to a single data file, rather than compiling the whole HMD/HFD. I prefer
# to have the whole thing in one chunk.

# change as necessary
setwd("/home/triffe/git/TimeProposal/TimeLab")
#library(devtools)
#install_github("DemogBerkeley", subdir = "DemogBerkeley", username = "UCBdemography")
library(DemogBerkeley)
library(reshape2)
countries <- getHFDcountries()
DAT <- do.call(rbind, lapply(countries,function(xyz,pw,us){ 
            DAT    <- readHFDweb(xyz,"birthsRR",password=pw,username=us) 
            BM     <- acast(DAT,Age~Year,value.var="Total")
            YA     <- matrix(0,nrow=12,ncol=ncol(BM),dimnames=list(0:11,colnames(BM)))
            OA     <- matrix(0,nrow=55,ncol=ncol(BM),dimnames=list(56:110,colnames(BM)))
            BM     <- rbind(YA,BM,OA)
            DAT    <- melt(BM, varnames=c("Age","Year"),value.name="Births")
            DAT$Code <- xyz
            DAT
        },pw=pw,us=us))
   ## save out for later use
save(DAT,file="Data/HFDall.Rdata")
   
# -------------------------------------------
# now get HMD exposures (differ from HFD exposures sometimes)

countries <- getHMDcountries()
DAT <- do.call(rbind, lapply(countries,function(xyz,pw,us){
            DAT             <- readHMDweb(xyz,"Exposures_1x1",password=pw,username=us) 
            Fem             <- DAT[,c("Year","Age","Female")]
            Mal             <- DAT[,c("Year","Age","Male")]
            colnames(Fem)   <- colnames(Mal) <- c("Year","Age","Exp")
            Fem$Sex         <- "f"
            Mal$Sex         <- "m"
            DAT             <- rbind(Fem,Mal)
            DAT$Code        <- xyz
            DAT
        },pw=pw,us=us))
save(DAT,file="Data/HMDExpall.Rdata")
 #and again for HMD lifetables (males grabbed superfluously)    
DAT <- do.call(rbind, lapply(countries,function(xyz,pw,us){
            DATm    <- readHMDweb(xyz,"mltper_1x1",password=pw,username=us) 
            DATf    <- readHMDweb(xyz,"fltper_1x1",password=pw,username=us) 
            DATm$Sex <- "m"
            DATf$Sex <- "f"
            DATi <- rbind(DATf,DATm)
            DATi$Code <- xyz
            DATi
        },pw=pw,us=us))
save(DAT,file="Data/HMDall.Rdata")
