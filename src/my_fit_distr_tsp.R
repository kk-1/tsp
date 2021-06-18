

########################################################################################################
#Clean the environment
remove(list = ls())
########################################################################################################

library(fitdistrplus)
library(logspline)
library(TSP)
library(psych)
library(DescTools)

########################################################################################################

myWorkDir <- "/my_hd1/my_dir/my_prg/Camerino/tsp"
myDataDir <- "/my_hd1/my_dir/my_prg/Camerino/tsp/data"
myPlotDir <- "/my_hd1/my_dir/my_prg/Camerino/tsp/fit-distr-plot"
setwd(myWorkDir)


if (!(dir.exists(myPlotDir ))){
  dir.create(myPlotDir , recursive = TRUE)
}

########################################################################################################
#Some Globals
#########################################################################################################



tspDataSets <- c(
  "att48"                    ,
  "berlin52"                 ,
  "pr76"                     ,
  "kroA100"                  ,
  "lin105"                   ,
  "ch130"                    ,
  "ch150"                    ,
  "a280"                     ,
  "pcb442"                   ,
  "Tnm52"                    ,
  "Tnm76"                    ,
  "Tnm100"                   ,
  "Tnm127"                   ,
  "Tnm154"                   ,
  "Tnm178"                   ,
  "Tnm199"                   ,
  "myRND-100"                ,
  "myRND-200"                ,
  "myRND-300"                ,
  "myRND-400"                ,
  "myLattice-10x10-100"          ,
  "myLattice-10x20-200"          ,
  "myLattice-15x20-300"          ,
  "myLattice-20x20-400"          ,
  "myRNDLattice-12x12-100"   ,
  "myRNDLattice-12x23-200"   ,
  "myRNDLattice-18x23-300"   ,
  "myRNDLattice-23x23-400"   ,
  "myHexLattice-10x10-100"       ,
  "myHexLattice-10x20-200"       ,
  "myHexLattice-15x20-300"       ,
  "myHexLattice-20x20-400"       ,
  "myRNDHexLattice-12x12-100",
  "myRNDHexLattice-12x23-200",
  "myRNDHexLattice-18x23-300",
  "myRNDHexLattice-23x23-400"
)

bigTSPDataSets <- c(
  "pr1002",
  "pr2392",
  "rl5915",
  "usa13509",
  "pla33810E"
)

bigCustomTSPDataSets <- c(
  "myHexLattice-100x100-10000",
  "myHexLattice-100x200-20000",
  "myHexLattice-150x200-30000",
  "myLattice-100x100-10000",
  "myLattice-100x200-20000",
  "myLattice-150x200-30000",
  "myRNDHexLattice-105x105-10000",
  "myRNDHexLattice-105x210-20000",
  "myRNDHexLattice-158x210-30000",
  "myRNDLattice-105x105-10000",
  "myRNDLattice-105x210-20000",
  "myRNDLattice-158x210-30000"
)

#########################################################################################################
#Some Globals
#########################################################################################################

plotWidth <- 8
plotHeight <- 5

#Select all or certain dataset group
#dsetname <- c(tspDataSets, bigTSPDataSets, bigCustomTSPDataSets )
#dsetname <- c(tspDataSets)
#dsetname <- c(bigCustomTSPDataSets)

dsetname <- c(bigCustomTSPDataSets)



tspNo <- 1


for (d in dsetname) {
  
  
  #DEBUG
  ########################################################################################################
  #Set the tsp file name here without .tsp!!!
  #theTSP <- "pr2392"  #OK
  ########################################################################################################
  
  
  #This is sily but helps for debugging!
  d <- dsetname[tspNo]
  theTSP <- d 
  
  cat(tspNo, " --- Distr Fit for", theTSP,"\n")
  
  theTSPplotRoot <- paste(myPlotDir,theTSP, sep = "/")
  theTSPFileName <- paste0(theTSP,".tsp")
  
  cat("Opening the TSP:", theTSP, "\n")
  theTSPObject <- read_TSPLIB(paste(myDataDir,theTSPFileName, sep = "/"))
  theTSPcoordDF <- as.data.frame(theTSPObject)
  #insert vertex ids
  theTSPcoordDF[,3] <- 1:nrow(theTSPcoordDF)
  colnames(theTSPcoordDF) <- c("x", "y", "vno")
  
  theTSPcoordMTX <- as.matrix(theTSPcoordDF)
  
  datasetLen <- nrow(theTSPcoordMTX)
  vecLen <- datasetLen * datasetLen
  
  #For plot axis
  xmin <- min(theTSPcoordMTX[,1])
  xmax <- max(theTSPcoordMTX[,1])
  
  ymin <- min(theTSPcoordMTX[,2])
  ymax <- max(theTSPcoordMTX[,2])
  
  
  theDistMTX <- as.matrix(dist(theTSPcoordMTX[,1:2], diag = TRUE, upper = TRUE))
  theDistVector <- as.vector(theDistMTX)
  theDistVectorNoZero <- theDistVector[theDistVector != 0] #Remove zeros
  
  outfileName <- (paste0(theTSPplotRoot,"-fitdistr-results.txt"))
  
  ########################################################################################################
  
  #describe(theDistVector)
  
  #t <- describe(theDistVector)
  #cat("Vec Stat:\n",paste(print(t))[0],"\n",file=outfileName,append=TRUE)
  
  
  #describe(theDistVectorNoZero)
  
  
  if (vecLen > 1000000) {
    #Big Data use sample
    theDistSampleNoZero <- sample(theDistVectorNoZero, 1000000, replace=FALSE)
   
    summary(theDistVectorNoZero)
    summary(theDistSampleNoZero)
    
    #Clean these big arrays from memory to avoid memory problems
    rm(theDistMTX)
    rm(theDistVector)
    rm(theDistVectorNoZero)
  }else{
    #Small Data use the distvector directly
    theDistSampleNoZero <- theDistVectorNoZero
    summary(theDistVectorNoZero)
    summary(theDistSampleNoZero)
  }
  
  
  
  #describe(theDistSampleNoZero)
  
  
  
  cairo_ps(filename = (paste0(theTSPplotRoot,"-descr-stat-plot.eps")), width = plotWidth, height = plotHeight, pointsize = 12, fallback_resolution = 600)
  # Desc(theDistVectorNoZero, plotit=TRUE)
  print(Desc(theDistSampleNoZero, plotit=TRUE, main = paste0(theTSP," - Descriptive Stats") ))
  dev.off()
  
  
  
 
  
  
  #May take too much time
  #par(mfrow = c(1, 1))
  #plotdist(theDistSampleNoZero, histo = TRUE, demp = TRUE)
  #descdist(theDistSampleNoZero, boot = 1000)
  
  
  
  # xnzp : theDistSampleNoZero proportion
  xnzp <- theDistSampleNoZero / max(theDistSampleNoZero)
  xnzp <- xnzp[xnzp > 0 & xnzp < 1] # Since Gamma assumes this range
  summary(xnzp)
  
  # xnzs : theDistSampleNoZero scaled
  xnzs <- (theDistSampleNoZero - min(theDistSampleNoZero)) / max(theDistSampleNoZero)
  xnzs <- xnzs[xnzs > 0 & xnzs < 1]
  summary(xnzs)
  
  summary(fbetp <- fitdist(xnzp, "beta", control=list(trace=1, REPORT=1)))
  summary(fbets <- fitdist(xnzs, "beta", control=list(trace=1, REPORT=1)))
 
  
  summary(fws <- fitdist(xnzs, "weibull"))
  summary(flns <- fitdist(xnzs, "lnorm"))
  summary(fgs <- fitdist(xnzs, "gamma"))
  summary(fns <- fitdist(xnzs, "norm"))
  
  plot.legend <- c("Weibull", "Lognormal", "Gamma", "Beta", "Normal")
  
  myCol <- c("red", "blue", "orange", "green","pink")
  
  cairo_ps(filename = (paste0(theTSPplotRoot,"-density-plot.eps")), width = plotWidth, height = plotHeight, pointsize = 12, fallback_resolution = 600)
  
  print(denscomp(list(fws, flns, fgs, fbets, fns), main=paste0(theTSP,"-density-plot"), legendtext = plot.legend, plotstyle = "ggplot", fitcol = myCol))
  
  dev.off()
  
  
  
  
  
  cairo_ps(filename = (paste0(theTSPplotRoot,"-qq-plot.eps")), width = plotWidth, height = plotHeight, pointsize = 12, fallback_resolution = 600)
  
  print(qqcomp(list(fws, flns, fgs, fbets, fns), main=paste0(theTSP,"-qq-plot"), legendtext = plot.legend, plotstyle = "ggplot", fitcol = myCol))
  
  dev.off()
  
  
  
  
  cairo_ps(filename = (paste0(theTSPplotRoot,"-cdf-plot.eps")), width = plotWidth, height = plotHeight, pointsize = 12, fallback_resolution = 600)
  
  print(cdfcomp(list(fws, flns, fgs, fbets, fns), main=paste0(theTSP,"-cdf-plot"), legendtext = plot.legend, plotstyle = "ggplot", fitcol = myCol))
  
  dev.off()
  
  
  
  cairo_ps(filename = (paste0(theTSPplotRoot,"-pp-plot.eps")), width = plotWidth, height = plotHeight, pointsize = 12, fallback_resolution = 600)
  
  print(ppcomp(list(fws, flns, fgs, fbets, fns), main=paste0(theTSP,"-pp-plot"), legendtext = plot.legend, plotstyle = "ggplot", fitcol = myCol))
  
  dev.off()
  
  
  fitResults <- gofstat(list(fws, flns, fgs, fbets, fns))
  # Gamma is the best distribution for xnzs series
  
  
  
  cat("------------------------------------------------------------------\n",file=outfileName,append=TRUE)
  cat("Dist Fit for", theTSP,"\n",file=outfileName,append=TRUE)
  cat("------------------------------------------------------------------\n",file=outfileName,append=TRUE)
  cat("Dists  :", plot.legend,"\n",file=outfileName,append=TRUE)
  cat("------------------------------------------------------------------\n",file=outfileName,append=TRUE)
  cat("KS-dec :", fitResults$kstest,"\n",file=outfileName,append=TRUE)
  cat("KS-val :", fitResults$ks,"\n",file=outfileName,append=TRUE)
  cat("CM-dec :", fitResults$cvmtest,"\n",file=outfileName,append=TRUE)
  cat("CM-val :", fitResults$cvm,"\n",file=outfileName,append=TRUE)
  cat("AD-dec :", fitResults$adtest,"\n",file=outfileName,append=TRUE)
  cat("AD-val :", fitResults$ad,"\n",file=outfileName,append=TRUE)
  cat("------------------------------------------------------------------\n",file=outfileName,append=TRUE)
  cat("AIC    :", fitResults$aic,"\n",file=outfileName,append=TRUE)
  cat("------------------------------------------------------------------\n",file=outfileName,append=TRUE)
  cat("Min KS  :", plot.legend[which.min(fitResults$ks)],"\n",file=outfileName,append=TRUE)
  cat("Min CM  :", plot.legend[which.min(fitResults$cvm)],"\n",file=outfileName,append=TRUE)
  cat("Min AD  :", plot.legend[which.min(fitResults$ad)],"\n",file=outfileName,append=TRUE)
  cat("Min AIC :", plot.legend[which.min(fitResults$aic)],"\n",file=outfileName,append=TRUE)
  cat("------------------------------------------------------------------\n",file=outfileName,append=TRUE)
  
  tspNo <- tspNo + 1
  
}









