
########################################################################################################
#Benchmark the algo and the data set n times and saves results into a text file
########################################################################################################
#Clean the environment
remove(list = ls())
########################################################################################################

myWorkDir <- "/my_hd1/my_dir/my_prg/Camerino/tsp"
myResultsDir <- "/my_hd1/my_dir/my_prg/Camerino/tsp/results-big-custom-cr=0.9"
myDataDir <- "/my_hd1/my_dir/my_prg/Camerino/tsp/data-big"
#########################################################################################################
library(TSP)
library(stringr)  
library(dplyr)
library(tictoc)
library(sf)
library(spatstat.geom)
library(gissr)
library(concaveman)
library(colorspace)
library(scales)
library(readr)
library(deldir)
library(tripack)
#########################################################################################################

setwd(myWorkDir)
source("tsp_functions-V3.R")


#########################################################################################################


########################################################################################################
#Arg processing
########################################################################################################
#Print the arguments for the Rscript
args <- commandArgs()
cat("Args:\n")
print(args)


########################################################################################################

########################################################################################################
#DEBUG
#quit()
##Set the values for interactive run here (uncomment the following 2 lines!):
#Arg Values start from the 7th element!!!
#Arguments comes as text!!!
#Comment them for script runs
########################################################################################################


########################################################################################################
#For interactive trial uncomment the next 2 lines
########################################################################################################
#args <- matrix(seq(1:11),nrow=1,ncol=11)
#args[7:10] <- c("repetitive_nn", "Tnm52", 1, 551609,1)
########################################################################################################

#Get the necessary values from args
theMethod <- args[7]
theTSP <- args[8]
nTrial <- as.numeric(args[9])
theOptimTourCost <- as.numeric(args[10]) #-1 if there is no tour file or optim tour cost is unknown
minAWDFlag <-  as.numeric(args[11]) #1: find minAWD or 0 : Do not find minAWD so time consuming




cat("--------------------------------------------------------------------------------------------\n")
cat("Testing",theMethod,"with",theTSP,nTrial,"times.\n")

########################################################################################################
#Create dir for each dataset

dataSetResultDir <- paste(myResultsDir,theTSP, sep = "/")

if (!(dir.exists(dataSetResultDir))){
  dir.create(dataSetResultDir, recursive = TRUE)
}

outputFileName <- paste0(theMethod,"-",nTrial,".csv")
outputFile <- paste(dataSetResultDir,outputFileName, sep = "/")
cat("Saving benchmark results to", outputFile,"\n")
########################################################################################################


coordDigitPrec <- 4

theTSPFileName <- paste0(theTSP,".tsp")
theTSPtourFileName <- paste0(theTSP,".opt.tour")
cat("The TSP is:", theTSP, "\n")
#########################################################################################################
#Some Globals
#########################################################################################################
theTSPObject <- read_TSPLIB(paste(myDataDir,theTSPFileName, sep = "/"))
theTSPcoordDF <- as.data.frame(theTSPObject)
#insert vertex ids
theTSPcoordDF[,3] <- 1:nrow(theTSPcoordDF)
colnames(theTSPcoordDF) <- c("x", "y", "vno")


tic.clearlog()


if (!file.exists(paste(myDataDir,theTSPtourFileName, sep = "/"))) {
  # Tour file is not there !!!
  cat("No tour file is found for", theTSP, "\n")
  theOptimTour <- NULL
  cat("The optim tour has cost:",theOptimTourCost,"\n")
  #theOptimTourCost <- -1
}else{
  cat("Tour file is found for", theTSP, "---",theTSPtourFileName,"\n")
  theOptimTour <- my_read_tsplib_tour(paste(myDataDir,theTSPtourFileName, sep = "/"))
  theOptimTourCost <- -1 * my_get_cost_perm(theOptimTour)
  cat("The optim tour has", length(theOptimTour), "vertices with total length:",theOptimTourCost,"\n")
}




ptsToGoDF <- theTSPcoordDF[,1:2]
ringPolyList <- list()
mergedRingPtsMTX <- NULL
theTSPcoordMTX <- as.matrix(theTSPcoordDF[,1:2])

theTSPcoordMTX <- as.matrix(theTSPcoordDF)






#############################################################################
#Delaunay edges stuff
#############################################################################


#############################################################################
#From tripack
#############################################################################
dxy <- tri.mesh(theTSPcoordMTX[,1], theTSPcoordMTX[,2])
DTEdgeList <- neighbours(dxy)

nList <- length(DTEdgeList)
nEdges <- length(unlist(DTEdgeList))
edgesDelaunayDF <- as.data.frame(matrix(0,nrow = nEdges, ncol = 2))
colnames(edgesDelaunayDF) <- c("src","dst")

addedEdge <- 0
for (i in 1:nList) {
  src <- i
  
  nNeighbor <- length(DTEdgeList[[i]])
  if (nNeighbor > 0){
    for (z in 1:nNeighbor) {
      dst <- DTEdgeList[[i]][z]
      addedEdge <- addedEdge + 1
      edgesDelaunayDF[addedEdge,] <- c(src,dst)
    }
  }
}

#DEBUG plot
#plot(dxy)

#DEBUG
# plot(theTSPObject, theOptimTour, tour_lty = 1, tour_col = "green")
# text(theTSPcoordMTX[,1], theTSPcoordMTX[,2], cex = 0.5, col = "red", pos = 2)
# plot(dxy, wlines = "triang", wpoints = "none", lty=2, lwd=0.2, add=T)

#############################################################################
#From tripack
#############################################################################


#############################################################################
#Delaunay edges stuff
#############################################################################












######################################################
#For plotting 
######################################################
xmin <- min(theTSPcoordDF[,1])
xmax <- max(theTSPcoordDF[,1])

ymin <- min(theTSPcoordDF[,2])
ymax <- max(theTSPcoordDF[,2])

myWidth <- 800
myHeight <- 600





##########################################################################################################
#Init - END
##########################################################################################################
algoNameData <- NULL
tourCostData <- NULL
aRatioData <- NULL
awdData <- NULL
minAWDData <- NULL
delEdgePercentData <- NULL

#Init stuff

theDistMTX <- as.matrix(dist(theTSPcoordMTX[,1:2], diag = TRUE, upper = TRUE))
diag(theDistMTX) <- Inf #So that min dist will not be the self loop











#The main benchmarking loop

for (i in 1:nTrial) {
  cat("---------------------------------------------------------------------------------------------\n")
  cat("Trial",i,"of",nTrial,"for",theMethod,"with",theTSP,"\n")
  cat("---------------------------------------------------------------------------------------------\n")
  
  
    if (theMethod == "nearest_insertion"){
    #standard TSP algos from R TSP
    
    #########################################################################################################
    #nearest_insertion
    #########################################################################################################
    tic("nearest_insertion")
    nearestInsOptimTour <- solve_TSP(theTSPObject, method = "nearest_insertion")
    toc(log = TRUE, quiet = TRUE)
    
    nearestInsOptimTourCost <- tour_length(nearestInsOptimTour)
    nearestInsApproxRatio <- round(nearestInsOptimTourCost/theOptimTourCost, 3)
    cat("The approximation ratio of nearest_insertion is: (",nearestInsOptimTourCost,"/", theOptimTourCost,") =",nearestInsApproxRatio,"\n")
    
    tourCostData <- c(tourCostData, nearestInsOptimTourCost)
    aRatioData <- c( aRatioData, nearestInsApproxRatio)
    
    
    nearestInsAWT <- my_get_awt_from_vtx1(nearestInsOptimTour)
    awdData <- c(awdData, nearestInsAWT)
    cat(theMethod,"has awt of:", nearestInsAWT,"\n")
    
    
    if (minAWDFlag == 1){
      nearestInsMinAWD <- my_get_min_awt_perm_fast(nearestInsOptimTour)
      minAWDData <- c(minAWDData,  nearestInsMinAWD[[1]])
    }else{
      minAWDData <- c(minAWDData, NA)
    }
    
    
    p <- my_get_del_percent_fast(nearestInsOptimTour, edgesDelaunayDF)
    cat("The Delaunay percentage is:",p,"\n")
    delEdgePercentData <- c(delEdgePercentData, p)
    #########################################################################################################
  } else if (theMethod == "farthest_insertion"){
    #########################################################################################################
    #farthest_insertion
    #########################################################################################################
    tic("farthest_insertion")
    farthestInsOptimTour <- solve_TSP(theTSPObject, method = "farthest_insertion")
    toc(log = TRUE, quiet = TRUE)
    
    farthestInsOptimTourCost <- tour_length(farthestInsOptimTour)
    farthestInsApproxRatio <-  round(farthestInsOptimTourCost/theOptimTourCost, 3)
    cat("The approximation ratio of farthest_insertion is: (",farthestInsOptimTourCost,"/", theOptimTourCost,") =",farthestInsApproxRatio,"\n")
    
    tourCostData <- c(tourCostData,  farthestInsOptimTourCost)
    aRatioData <- c( aRatioData,  farthestInsApproxRatio)
    
    farthestInsAWT <- my_get_awt_from_vtx1(farthestInsOptimTour)
    awdData <- c(awdData, farthestInsAWT)
    cat(theMethod,"has awt of:", farthestInsAWT,"\n")
    
    if (minAWDFlag == 1){
      farthestInsMinAWD <- my_get_min_awt_perm_fast(farthestInsOptimTour)
      minAWDData <- c(minAWDData,  farthestInsMinAWD[[1]])
    }else{
      minAWDData <- c(minAWDData, NA)
    }

    p <- my_get_del_percent_fast(farthestInsOptimTour, edgesDelaunayDF)
    cat("The Delaunay percentage is:",p,"\n")
    delEdgePercentData <- c(delEdgePercentData, p)
    #########################################################################################################
  } else if (theMethod == "cheapest_insertion"){
    #########################################################################################################
    #cheapest_insertion
    #########################################################################################################
    tic("cheapest_insertion")
    cheapestInsOptimTour <- solve_TSP(theTSPObject, method = "cheapest_insertion")
    toc(log = TRUE, quiet = TRUE)
    
    cheapestInsOptimTourCost <- tour_length(cheapestInsOptimTour)
    cheapestInsApproxRatio <-  round(cheapestInsOptimTourCost/theOptimTourCost, 3)
    cat("The approximation ratio of cheapest_insertion is: (",cheapestInsOptimTourCost,"/", theOptimTourCost,") =",cheapestInsApproxRatio,"\n")
    
    tourCostData <- c(tourCostData,  cheapestInsOptimTourCost)
    aRatioData <- c( aRatioData,  cheapestInsApproxRatio)
    
    cheapestInsAWT <- my_get_awt_from_vtx1(cheapestInsOptimTour)
    awdData <- c(awdData, cheapestInsAWT)
    cat(theMethod,"has awt of:", cheapestInsAWT,"\n")
    
    if (minAWDFlag == 1){
      cheapestInsMinAWD <- my_get_min_awt_perm_fast(cheapestInsOptimTour)
      minAWDData <- c(minAWDData, cheapestInsMinAWD[[1]])
    }else{
      minAWDData <- c(minAWDData, NA)
    }
    
        
    p <- my_get_del_percent_fast(cheapestInsOptimTour, edgesDelaunayDF)
    cat("The Delaunay percentage is:",p,"\n")
    delEdgePercentData <- c(delEdgePercentData, p)
    #########################################################################################################
  } else if (theMethod == "arbitrary_insertion"){
    #########################################################################################################
    #arbitrary_insertion
    #########################################################################################################
    tic("arbitrary_insertion")
    arbitraryInsOptimTour <- solve_TSP(theTSPObject, method = "arbitrary_insertion")
    toc(log = TRUE, quiet = TRUE)
    
    arbitraryInsOptimTourCost <- tour_length(arbitraryInsOptimTour)
    arbitraryInsApproxRatio <- round(arbitraryInsOptimTourCost/theOptimTourCost, 3)
    cat("The approximation ratio of arbitrary_insertion is: (",arbitraryInsOptimTourCost,"/", theOptimTourCost,") =",arbitraryInsApproxRatio,"\n")
    
    tourCostData <- c(tourCostData, arbitraryInsOptimTourCost)
    aRatioData <- c( aRatioData,  arbitraryInsApproxRatio)
    
    arbitraryInsAWT <- my_get_awt_from_vtx1(arbitraryInsOptimTour)
    awdData <- c(awdData, arbitraryInsAWT)
    cat(theMethod,"has awt of:", arbitraryInsAWT,"\n")
    
    if (minAWDFlag == 1){
      arbitraryInsMinAWD <- my_get_min_awt_perm_fast(arbitraryInsOptimTour)
      minAWDData <- c(minAWDData,  arbitraryInsMinAWD[[1]])
    }else{
      minAWDData <- c(minAWDData, NA)
    }
 
    
    p <- my_get_del_percent_fast(arbitraryInsOptimTour, edgesDelaunayDF)
    cat("The Delaunay percentage is:",p,"\n")
    delEdgePercentData <- c(delEdgePercentData, p)
    #########################################################################################################
  } else if (theMethod == "nn"){
    #########################################################################################################
    #nn
    #########################################################################################################
    tic("nn")
    nnOptimTour <- solve_TSP(theTSPObject, method = "nn")
    toc(log = TRUE, quiet = TRUE)
    
    nnOptimTourCost <- tour_length(nnOptimTour)
    nnApproxRatio <- round(nnOptimTourCost/theOptimTourCost, 3)
    cat("The approximation ratio of nn is: (",nnOptimTourCost,"/", theOptimTourCost,") =",nnApproxRatio,"\n")
    
    tourCostData <- c(tourCostData, nnOptimTourCost)
    aRatioData <- c( aRatioData,  nnApproxRatio)
    
    nnAWT <- my_get_awt_from_vtx1(nnOptimTour)
    awdData <- c(awdData,  nnAWT)
    cat(theMethod,"has awt of:", nnAWT,"\n")
    
    if (minAWDFlag == 1){
      nnMinAWD <- my_get_min_awt_perm_fast(nnOptimTour)
      minAWDData <- c(minAWDData, nnMinAWD[[1]])
    }else{
      minAWDData <- c(minAWDData, NA)
    }
 
    
    p <- my_get_del_percent_fast(nnOptimTour, edgesDelaunayDF)
    cat("The Delaunay percentage is:",p,"\n")
    delEdgePercentData <- c(delEdgePercentData, p)
    #########################################################################################################
  } else if (theMethod == "repetitive_nn"){
    #########################################################################################################
    #repetitive_nn
    #########################################################################################################
    tic("repetitive_nn")
    repetitive_nnOptimTour <- solve_TSP(theTSPObject, method = "repetitive_nn")
    toc(log = TRUE, quiet = TRUE)
    
    repetitive_nnOptimTourCost <- tour_length(repetitive_nnOptimTour)
    repetitive_nnApproxRatio <- round(repetitive_nnOptimTourCost/theOptimTourCost, 3)
    cat("The approximation ratio of repetitive_nn is: (",repetitive_nnOptimTourCost,"/", theOptimTourCost,") =",repetitive_nnApproxRatio,"\n")
    
    tourCostData <- c(tourCostData, repetitive_nnOptimTourCost)
    aRatioData <- c( aRatioData,  repetitive_nnApproxRatio)
    
    
    repetitive_nnAWT <- my_get_awt_from_vtx1(repetitive_nnOptimTour)
    awdData <- c(awdData,  repetitive_nnAWT)
    cat(theMethod,"has awt of:", repetitive_nnAWT,"\n")
    
    if (minAWDFlag == 1){
      repetitive_nnMinAWD <- my_get_min_awt_perm_fast(repetitive_nnOptimTour)
      minAWDData <- c(minAWDData, repetitive_nnMinAWD[[1]])
    }else{
      minAWDData <- c(minAWDData, NA)
    }
   
    
    p <- my_get_del_percent_fast(repetitive_nnOptimTour, edgesDelaunayDF)
    cat("The Delaunay percentage is:",p,"\n")
    delEdgePercentData <- c(delEdgePercentData, p)
    #########################################################################################################
  } else if (theMethod == "two_opt"){
    #########################################################################################################
    #two_opt
    #########################################################################################################
    tic("two_opt")
    two_optOptimTour <- solve_TSP(theTSPObject, method = "two_opt")
    toc(log = TRUE, quiet = TRUE)
    
    two_optOptimTourCost <- tour_length(two_optOptimTour)
    two_optApproxRatio <- round(two_optOptimTourCost/theOptimTourCost, 3)
    cat("The approximation ratio of two_opt is: (",two_optOptimTourCost,"/", theOptimTourCost,") =",two_optApproxRatio,"\n")
    
    tourCostData <- c(tourCostData, two_optOptimTourCost)
    aRatioData <- c( aRatioData,  two_optApproxRatio)
    
    
    two_optAWT <- my_get_awt_from_vtx1(two_optOptimTour)
    awdData <- c(awdData,  two_optAWT)
    cat(theMethod,"has awt of:", two_optAWT,"\n")
    
    if (minAWDFlag == 1){
      two_optMinAWD <- my_get_min_awt_perm_fast(two_optOptimTour)
      minAWDData <- c(minAWDData, two_optMinAWD[[1]])
    }else{
      minAWDData <- c(minAWDData, NA)
    }
    

    
    p <- my_get_del_percent_fast(two_optOptimTour, edgesDelaunayDF)
    cat("The Delaunay percentage is:",p,"\n")
    delEdgePercentData <- c(delEdgePercentData, p)
    #########################################################################################################
  }
  #standard TSP algos from R TSP
  
  
}






#########################################################################################################
#Timings
#########################################################################################################

log.txt <- tic.log(format = TRUE)
log.lst <- tic.log(format = FALSE)
tic.clearlog()
timings <- unlist(lapply(log.lst, function(x) x$toc - x$tic))
writeLines(unlist(log.txt))




#########################################################################################################
#ALL RESULTS TOGETHER
#########################################################################################################


timingData <- as.vector(timings)

if (is.null(theOptimTour)) {
  theOptimTourLen <- nrow(theTSPcoordDF)
}else{
  theOptimTourLen <- length(theOptimTour)
}

#On the last row put averages
timingData <- c(timingData, mean(timingData), sd(timingData) )
aRatioData <- c(aRatioData, mean(aRatioData), sd(aRatioData) )
tourCostData <- c(tourCostData, mean(tourCostData), sd(tourCostData))
awdData <- c(awdData, mean(awdData), sd(awdData))
minAWDData <- c(minAWDData, mean(minAWDData), sd(minAWDData))
delEdgePercentData <- c(delEdgePercentData, mean(delEdgePercentData), sd(delEdgePercentData))


theIdxCol <- seq(1:nTrial)
theIdxCol <- c(theIdxCol, "AVG")
theIdxCol <- c(theIdxCol, "STD")

algoNameData <- rep(theMethod, nTrial+2)


if (theOptimTourCost == -1){
  resultsDF <- data.frame(cbind(algoNameData,
                                theIdxCol,
                                round(timingData, 3), 
                                round(tourCostData, 3), 
                                round(awdData, 3),
                                round(minAWDData, 3),
                                round(delEdgePercentData, 3)
  )
  )
  
  colnames(resultsDF) <- c("Algo", "ResultId", "Time(sec)", "TourCost", "AWDfromVTX1", "MinAWD", "DelEdgeP")
}else{
  resultsDF <- data.frame(cbind(algoNameData,
                                theIdxCol,
                                round(timingData, 3), 
                                round(aRatioData,3),
                                round(tourCostData, 3), 
                                round(awdData, 3),
                                round(minAWDData, 3),
                                round(delEdgePercentData, 3)
  )
  )
  
  colnames(resultsDF) <- c("Algo", "ResultId", "Time(sec)", "ApproxRatio", "TourCost", "AWDfromVTX1", "MinAWD", "DelEdgeP")
}

write_excel_csv(resultsDF,outputFile)
cat("--------------------------------------------------------------------------------------------\n")



#########################################################################################################










