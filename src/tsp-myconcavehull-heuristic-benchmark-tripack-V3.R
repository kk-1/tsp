########################################################################################################
#concaveRingMergeNearestPtV2 - This is concave ring merging algo WITH heuristics 
#for the vertex insertion. The nearest point is used to select the best ring pt to merge
########################################################################################################
########################################################################################################
#Clean the environment
remove(list = ls())
########################################################################################################

myWorkDir <- "/my_hd1/my_dir/my_prg/Camerino/tsp"
myResultsDir <- "/my_hd1/my_dir/my_prg/Camerino/tsp/results-big-custom-cr=0.9"
myDataDir <- "/my_hd1/my_dir/my_prg/Camerino/tsp/data-big"


myPlotDir <- "/my_hd1/my_dir/my_prg/Camerino/tsp/plot"

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
library(readr)
library(deldir)
library(tripack)
#########################################################################################################

setwd(myWorkDir)
source("tsp_functions-V3.R")


#########################################################################################################
##########################
#Tested tsps
theTSP <- "att48"  #OK
theTSP <- "berlin52"  #OK
theTSP <- "ch150"  #OK
theTSP <- "pcb442" #OK
theTSP <- "pr2392" #OK
theTSP <- "rl5915" ; theOptimTourCost <- 565530#OK
theTSP <- "usa13509";  theOptimTourCost <- 19982859#OK
theTSP <- "pla33810" ; theOptimTourCost <- 66043099 #OK
theTSP <- "usa115475"  #OK But theDistMTX not enough
##########################


##########################
theTSP <- "Tnm52"  ; theOptimTourCost <- 551609#OK
theTSP <- "Tnm199"  ; theOptimTourCost <- 3139778 #OK
theTSP <- "ch150"  #OK
theTSP <- "myHexLattice-3434"; theOptimTourCost <- -1 #OK
theTSP <- "myHexLattice-672"; theOptimTourCost <- -1 #OK
theTSP <- "myHexLattice-1512"; theOptimTourCost <- -1 #OK
theTSP <- "myHexLattice-6163"; theOptimTourCost <- -1 #OK
theTSP <- "att48"  #OK
########################################################################################################






########################################################################################################
##########################
#Copy one of the previous line as next line
theTSP <- "myLattice-20x20"; theOptimTourCost <- -1 #OK
theTSP <- "myHexLattice-672"; theOptimTourCost <- -1 #OK
theTSP <- "myHexLattice-1512"; theOptimTourCost <- -1 #OK
theTSP <- "myRND-100"; theOptimTourCost <- -1 #OK
theTSP <- "myRND-200"; theOptimTourCost <- -1 #OK
theTSP <- "pcb442" #OK



#######################################
#Benchmark order for the paper results
######################################
theTSP <- "att48"  #OK
theTSP <- "berlin52"  #OK
theTSP <- "pr76"  #OK
theTSP <- "ch130"  #OK
theTSP <- "ch150"  #OK
theTSP <- "a280"  #OK
theTSP <- "Tnm52"  ; theOptimTourCost <- 551609#OK
theTSP <- "Tnm100"  ; theOptimTourCost <- 1398070 #OK
theTSP <- "Tnm199"  ; theOptimTourCost <- 3139778 #OK
theTSP <- "myRND-50"; theOptimTourCost <- -1 #OK
theTSP <- "myRND-150"; theOptimTourCost <- -1 #OK
theTSP <- "myRND-250"; theOptimTourCost <- -1 #OK
theTSP <- "myLattice-10x10"; theOptimTourCost <- -1 #OK
theTSP <- "myLattice-10x20"; theOptimTourCost <- -1 #OK
theTSP <- "myLattice-20x20"; theOptimTourCost <- -1 #OK
theTSP <- "myLattice-20x30"; theOptimTourCost <- -1 #OK
theTSP <- "myHexLattice-10x10"; theOptimTourCost <- -1 #OK
theTSP <- "myHexLattice-10x20"; theOptimTourCost <- -1 #OK
theTSP <- "myHexLattice-20x20"; theOptimTourCost <- -1 #OK
theTSP <- "myHexLattice-20x30"; theOptimTourCost <- -1 #OK

theTSP <- "att48"  #OK
theTSP <- "Tnm52"  ; theOptimTourCost <- 551609#OK
theTSP <- "myLattice-20x20"; theOptimTourCost <- -1 #OK
theTSP <- "myHexLattice-20x20"; theOptimTourCost <- -1 #OK


theTSP <- "mydantzig42"; theOptimTourCost <- 699 #OK
cr <- 0.7
##########################


#########################################################################################################


########################################################################################################
#Arg processing
########################################################################################################
#Print the arguments for the Rscript
#cat("Here is the usage:\n")
#cat("Rscript --vanilla tsp-benchmark-algo-V3.0.R <algo> <tspdatafile> <n>\n")
args <- commandArgs()
cat("Args:\n")
print(args)


########################################################################################################

########################################################################################################
#DEBUG
#quit()
##Set the values for interactive run here (uncomment the following 2 lines!):
#Arg Values start from the 7th element!!!
# #Arguments comes as text!!!
#Comment them for script runs
########################################################################################################
#args <- matrix(seq(1:9),nrow=1,ncol=9)
#args[7:9] <- c("concaveRingMergeNearestPtD", "myLattice-20x20", 1)
########################################################################################################


########################################################################################################
#For interactive trial uncomment the next 2 lines
########################################################################################################
#args <- matrix(seq(1:10),nrow=1,ncol=10)
#args[7:10] <- c("repetitive_nn", "Tnm52", 1, 551609)
#args <- matrix(seq(1:10),nrow=1,ncol=10)
#args[7:10] <- c("Tnm199", 3, 3139778, 0.9)
# 
# args <- matrix(seq(1:11),nrow=1,ncol=11)
# args[7:11] <- c("pr76", 1, -1, 0.9, 1)
# 
# args <- matrix(seq(1:11),nrow=1,ncol=11)
# args[7:11] <- c("myLattice-20x20-400", 1, -1, 0.9, 0)
########################################################################################################

#Get the necessary values from args
theTSP <- args[7]
nTrial <- as.numeric(args[8])
theOptimTourCost <- as.numeric(args[9]) #-1 if there is no tour file or optim tour cost is unknown
cr <- as.numeric(args[10])
minAWDFlag <-  as.numeric(args[11]) #1: find minAWD or 0 : Do not find minAWD so time consuming
theMethod <- "concaveRingMergeNearestPtV2"


cat("--------------------------------------------------------------------------------------------\n")
cat("Testing",theMethod,"with",theTSP,nTrial,"times cr=",cr,"\n")



########################################################################################################
#Create dir for each dataset

dataSetResultDir <- paste(myResultsDir,theTSP, sep = "/")

if (!(dir.exists(dataSetResultDir))){
  dir.create(dataSetResultDir, recursive = TRUE)
}

outputFileName <- paste0(theMethod,"-",nTrial,"-cr-",cr,".csv")
outputFile <- paste(dataSetResultDir,outputFileName, sep = "/")
cat("Saving benchmark results to", outputFile,"\n")
########################################################################################################

#ATTENTION concaveman rounds numbers to 4 prec digits
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


if (!file.exists(paste(myDataDir,theTSPtourFileName, sep = "/"))) {
  # Tour file is not there !!!
  cat("No tour file is found for", theTSP, "\n")
  theOptimTour <-  NULL
}else{
  cat("Tour file is found for", theTSP, "---",theTSPtourFileName,"\n")
  theOptimTour <- my_read_tsplib_tour(paste(myDataDir,theTSPtourFileName, sep = "/"))
  theOptimTourCost <- -1 * my_get_cost_perm(theOptimTour)
  cat("The optim tour has", length(theOptimTour), "vertices with total length:",theOptimTourCost,"\n")
}



theTSPcoordMTX <- as.matrix(theTSPcoordDF)


#############################################################################################################
#DEBUG plot
#############################################################################################################
# #############################################################################################################
# #Uncomment the following to see the plot of the optimtour
# #############################################################################################################
# #DEBUG
# 
# 
# 
# 
# xmin <- min(theTSPcoordMTX[,1])
# xmax <- max(theTSPcoordMTX[,1])
# 
# ymin <- min(theTSPcoordMTX[,2])
# ymax <- max(theTSPcoordMTX[,2])
# 
# #############################################################################################################
# #For EPS output
# cairo_ps(filename = (paste0(theTSP,"-vtx.eps")),
#          width = 12,
#          height = 7,
#          pointsize = 12,
#          fallback_resolution = 600)
# 
# #Draw the tsp
# plot(theTSPObject,
#      xlab = "x-coordinate",
#      ylab = "y-coordinate",
#      lty = 3, lwd=0.7, cex=1.5,  col = "black", pch= 19, asp=1,
#      main = paste0(theTSP," --- Vertices"),
#      xlim=c(xmin, xmax),
#      ylim=c(ymin, ymax)
# )
# 
# #Specific vertex label positions for pr76
# text(theTSPcoordMTX[c(1:71,73:75),1], theTSPcoordMTX[c(1:71,73:75),2], labels = c(1:71,73:75), cex = 0.8, col = "red", pos = 3)
# text(theTSPcoordMTX[c(72,76),1], theTSPcoordMTX[c(72,76),2], labels = c(72,76), cex = 0.8, col = "red", pos = 4)
# 
# # pos	
# # a position specifier for the text. If specified this overrides any adj value given. 
# # Values of 1, 2, 3 and 4, respectively indicate positions 
# # below, to the left of, above and to the right of the specified (x,y) coordinates.
# 
# dev.off()
#############################################################################################################


# #############################################################################################################
# 














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









#############################################################################################################
#DEBUG plot
#############################################################################################################

## Transparent colors
## Mark Gardener 2015
## www.dataanalytics.org.uk

# t_col <- function(color, percent = 50, name = NULL) {
#   #      color = color name
#   #    percent = % transparency
#   #       name = an optional name for the color
#   
#   ## Get RGB values for named color
#   rgb.val <- col2rgb(color)
#   
#   ## Make new color using input color as base and alpha set by transparency
#   t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
#                max = 255,
#                alpha = (100 - percent) * 255 / 100,
#                names = name)
#   
#   return(t.col)
#   ## Save the color
#   #invisible(t.col)
# }
# ## END
# myGreenTrans <- t_col("green", perc = 25)
# #############################################################################################################
# #For EPS output
# cairo_ps(filename = (paste0(theTSP,"-dt.eps")),
#          width = 12,
#          height = 7,
#          pointsize = 12,
#          fallback_resolution = 600)
# #Draw the tsp
# plot(theTSPObject, theOptimTour, 
#      xlab = "x-coordinate",
#      ylab = "y-coordinate",
#      lty = 3, lwd=1, cex=1.5,  col = "black", pch= 19, asp=1,
#      main = paste0(theTSP," --- Optimum tour and Delaunay Triangulation edges"),
#      xlim=c(xmin, xmax),
#      ylim=c(ymin, ymax)
# )
# plot(dxy, wlines = "triang", wpoints = "none", lty=1, lwd=2, col = myGreenTrans, asp=1, add=T)
#     
# points(theTSPcoordMTX[,1:2], lty = 3, lwd=0.7, cex=1.5,  col = "black", pch= 19, asp=1)
# #Specific vertex label positions for pr76
# text(theTSPcoordMTX[c(1:71,73:75),1], theTSPcoordMTX[c(1:71,73:75),2], labels = c(1:71,73:75), cex = 0.8, col = "red", pos = 3)
# text(theTSPcoordMTX[c(72,76),1], theTSPcoordMTX[c(72,76),2], labels = c(72,76), cex = 0.8, col = "red", pos = 4)
# 
# 
# dev.off()
#############################################################################################################








#############################################################################
#Delaunay edges stuff
#############################################################################








#############################################################################


algoNameData <- NULL
tourCostData <- NULL
aRatioData <- NULL
awdData <- NULL
minAWDData <- NULL
delEdgePercentData <- NULL

tic.clearlog()


#########################################################################################################
#BEGIN MY ALGO
#########################################################################################################

#Init stuff
tic("DistMtx")
theDistMTX <- as.matrix(dist(theTSPcoordMTX[,1:2], diag = TRUE, upper = TRUE))
diag(theDistMTX) <- Inf #So that min dist will not be the self loop
t <- toc()
DistMTXSecs <- t$toc - t$tic

cat(DistMTXSecs,"secs for dmtx creation.\n")

tic.clearlog()


for (i in 1:nTrial) {
  
  #DEBUG i <- 1
  
  cat("---------------------------------------------------------------------------------------------\n")
  cat("Trial",i,"of",nTrial,"for",theMethod,"with",theTSP,"cr=",cr,"\n")
  cat("---------------------------------------------------------------------------------------------\n")
  tic(theMethod)
  mergedRingPtsMTX <- NULL
  ptsToGoMTX <- theTSPcoordMTX
  #in polygons we will not have vertex ids but in PtsDF!!!
  ringMTXList <- list()
  #########################################################################################################
  #The loop for deleting chull repeatedly till all the points are covered
  #########################################################################################################
  
  
  iter <- 1
  while (nrow(ptsToGoMTX) > 2) {
    
    #Attention concaveman returns closed-cyclical mtx
    #All depends on that cr value
    #As it goes less than 1 more points and less number of rings 
    
    #chullMTX is cyclical!!!
    chullMTX <- concaveman(ptsToGoMTX, concavity = cr)
    nChullPts <- nrow(chullMTX) - 1
    if (nChullPts <= 2){
      #It is possible there are hundreds of points linear
      #In this case the longest line with the two end pointas will be the concavehull
      #So we will have 2 points
      #Time to break
      cat("Iter",iter,"breaking the loop!\n")
      break
    }
    
    
    
    #RingPtsMTXcontains all the points on ring no repetition of the first vtx at the end
    RingPtsMTX <- chullMTX[1:(nrow(chullMTX)-1),]
    ringMTXList[[iter]] <- chullMTX
    #Precision figures can be problem in join so lets normalize them
    RingPtsMTX <- round(RingPtsMTX, coordDigitPrec)
    ptsToGoMTX <- round(ptsToGoMTX, coordDigitPrec)
    #Now subtract the ring  points from the TSP pts and find new ring to do the same
    
    #Both  ptsToGoMTX and RingPtsMTXare not cyclical!!!
    
    #ptsToGoMTX <- anti_join(ptsToGoMTX, RingPtsMTX, by = c("x" = "x", "y" = "y", "vno" = "vno"))
    ptsToGoDF <- as.data.frame(ptsToGoMTX)
    colnames(ptsToGoDF) <-  c("x", "y", "vno")
    
    RingPtsDF <- as.data.frame(RingPtsMTX)
    colnames(RingPtsDF) <-  c("x", "y", "vno")
    
    ptsToGoMTX <- as.matrix(anti_join(ptsToGoDF, RingPtsDF, by = c("x" = "x", "y" = "y", "vno" = "vno")))
    
    iter <- iter + 1
  }
  
  ####################################################################
  #DO NOT FORGET TO ADD THEM!!!
  ptsNotInAnyRingMTX <- ptsToGoMTX
  ####################################################################
  lenRemained <- nrow(ptsToGoMTX)
  
  if (lenRemained > 0) { 
    cat(lenRemained,"pts remained that are not in any chull.\n")
  }else{
    cat("All pts are in chulls.\n")
  }
  
  nRings <- length(ringMTXList)
  
  cat(nRings, "rings generated.\n")
  
  
  
  #############################################################################################################
  #DEBUG plot
  #########################################################################################################
  #BEGIN plot all rings together uncomment the following lines
  #########################################################################################################  
  
  # ringCol <- qualitative_hcl(nRings, palette = "Dark 3")
  # 
  # ####################################################################
  # #Plot all rings on the same plot
  # ####################################################################
  # 
  # #For EPS output
  # cairo_ps(filename = (paste0(theTSP,"-rings.eps")),
  #          width = 12,
  #          height = 7,
  #          pointsize = 12,
  #          fallback_resolution = 600)
  # 
  # #Draw the tsp
  # plot(theTSPObject,
  #      xlab = "x-coordinate",
  #      ylab = "y-coordinate",
  #      lty = 3, lwd=0.7, cex=1.5,  col = "black", pch= 19, asp=1,
  #      main = paste0(theTSP," --- ",nRings," Rings"),
  #      xlim=c(xmin, xmax),
  #      ylim=c(ymin, ymax)
  # )
  # 
  # #Specific vertex label positions for pr76
  # text(theTSPcoordMTX[c(1:71,73:75),1], theTSPcoordMTX[c(1:71,73:75),2], labels = c(1:71,73:75), cex = 0.8, col = "blue", pos = 3)
  # text(theTSPcoordMTX[c(72,76),1], theTSPcoordMTX[c(72,76),2], labels = c(72,76), cex = 0.8, col = "blue", pos = 4)
  # 
  # 
  # 
  # for (z in 1:nRings) {
  #   cat("Drawing ring:",z,"\n")
  #   theRingMTX <- ringMTXList[[z]]
  #   
  #   
  #   
  #   np <- nrow(theRingMTX)
  #   theRingPtsMTX <-  theRingMTX[1:(np-1),]
  #   
  #   ringTour <-  theRingPtsMTX[,3]
  #   p <- my_get_del_percent_fast(ringTour, edgesDelaunayDF)
  #   cat("Ring",z, "has Delaunay percentage of",p,"\n")
  #   
  #   #Draw the merged ring 
  #   #colx =  rainbow(nRings)[sample(nRings,1)]
  #   colx = ringCol[z]
  #   points(theRingPtsMTX[,1], theRingPtsMTX[,2], col = colx, pch = 19, cex=2)
  #   lines(theRingMTX[,1], theRingMTX[,2], col = colx, pch = 19, cex=2, asp=1,)
  #   text(theRingPtsMTX[,1]+1, theRingPtsMTX[,2]+1, labels=1:nrow(theRingPtsMTX), cex=0.75, font=2)
  # }
  # 
  # #Add also points that are remained as black dots
  # if (lenRemained > 0) { 
  #   points(ptsNotInAnyRingMTX[,1], ptsNotInAnyRingMTX[,2],   col = "black", pch=19, cex=2)
  # }
  # 
  # dev.off()
  ####################################################################
  #########################################################################################################
  #END plot all rings together 
  #########################################################################################################   
  
  
  #########################################################################################################
  #Now unify the rings to single one: One ring to rule them all!!!
  #########################################################################################################
  # Watch the boundary cases!!!
  # When one of the pair is first or the last vertex you need to erase the closing duplicate vertex stuff 
  # and update the rings (polygons or chulls)!!!
  #Loop over rings
  # The first ring is the outmost one
  # All others are inside of it
  # In this loop slowly merge them to the outermost ring
  #ATTENTION:
  #The last and the first pts are same for polygons
  #For all other structures throw the last pts away
  
  
  #Put the vertices of the outmost ring to the "merged Ring" (the ruler ring!!!)
  mergedRingMTX <- ringMTXList[[1]]
  nr <- nrow(mergedRingMTX)
  mergedRingPtsMTX <<- mergedRingMTX[1:(nr-1),]
  #Pts of the mergedRingPtsMTX so far
  nPts <- nrow(mergedRingPtsMTX)
  #Save the vertex ids and update this vector as more pts are added
  mergedRingPtsVnoVec <- c(mergedRingPtsMTX[, 3])  
  
  
  
  #Start from the second ring, if any ring left after concave hull, and add them into the merged one
  
  if (nRings >= 2){
    for (r in 2:nRings) {
      
      #DEBUG
      #r <- 2
      
      
      #Put ring points into a buffer as they  will be updated later and the id numbers will change
      ringToBeMergedPtsMTX <- ringMTXList[[r]]
      #Take out the last points as it is the same with the first one
      ringToBeMergedPtsMTX <- ringToBeMergedPtsMTX[1:(nrow(ringToBeMergedPtsMTX)-1), ]
      nPtsToBeMerged <- nrow(ringToBeMergedPtsMTX)
      
      
      #Loop over points on the next -- inner ring
      #Some ideas to think:
      #The problem is where to start and should we connect every point individually or in chunks of chain?
      #Starting from the closest points?
      #If the point is furthest away from a threshold skip it?
      for (p in (1:nPtsToBeMerged)) {
        
        #DEBUG
        #p <- 5
        
        #DEBUG
        #p <- 1
        
        
        ################################################################################## 
        #DEBUG plot
        # mergedPts <- matrix(ncol=3)
        ################################################################################## 
        
        
        
        ##################################################################################  
        #Merging criteria: Closest point from the MergedRing
        ##################################################################################
        #Try to find the nearest pt from the merged ring
        #ATTENTION: The merged ring changes after each point addition
        #Number of distinct points in the ring,
        #excluding the last one since it is same with the first one (polygons and rings are closed)
        n <- nrow(mergedRingPtsMTX)
        
        PtToBeMerged <- ringMTXList[[r]][p,]
        PtToBeMergedVno <- ringMTXList[[r]][p,3]
        
        #The min pt should also be on the merged ring
        minPtidx <- which.min(theDistMTX[PtToBeMergedVno, mergedRingPtsVnoVec])
        minPtVno <- mergedRingPtsVnoVec[minPtidx]
        
        nearestPtsIdx <- minPtidx
        nearestPtsVno <- minPtVno
        
        
        
        ##################################################################################
        #Now we find the nearest pt on the mergedRing to connect the pt.
        
        
        if (nearestPtsIdx == 1) {
          #The nearest pt is the first pt in mergedRing
          #The previous point is the penultimate point of the merged ring
          #Attention the last one and the first one is same!!!
          preNearestIdx <- n
          preNearestVno <- mergedRingPtsMTX[n, 3]
          postNearestIdx <- 2 
          postNearestVno <- mergedRingPtsMTX[2, 3]
          
          #Get the distances to see which way is the cheapest to connect!!!
          #Before the nearest pt or after the nearest pt!!!
          #Distance between pre and the nearest
          preDist <- theDistMTX[nearestPtsVno,  preNearestVno]
          #Distance between post and the nearest
          postDist <- theDistMTX[nearestPtsVno,  postNearestVno]
          #Distance between the point to be merged and the nearest
          nearest2ptDist <- theDistMTX[PtToBeMergedVno, nearestPtsVno]
          #Distance between the pre of the nearest and pt to be merged
          pre2ptDist <- theDistMTX[PtToBeMergedVno, preNearestVno]
          #Distance between the post of the nearest and pt to be merged
          post2ptDist <- theDistMTX[PtToBeMergedVno, postNearestVno]
          
          if ((preDist + post2ptDist) < (postDist + pre2ptDist)){
            #after the nearest
            #DEBUG
            #cat("1: nearest first after nearest\n")
            mergedRingPtsMTX <<- rbind(mergedRingPtsMTX[1:nearestPtsIdx,], PtToBeMerged, mergedRingPtsMTX[(nearestPtsIdx + 1):n,] )
            #Update vnovector too
            mergedRingPtsVnoVec <- c(mergedRingPtsVnoVec[1:nearestPtsIdx], PtToBeMergedVno, mergedRingPtsVnoVec[(nearestPtsIdx + 1):n])
          }else{
            #before the nearest
            #DEBUG
            #cat("2: nearest first before nearest\n")
            mergedRingPtsMTX <<- rbind(PtToBeMerged, mergedRingPtsMTX[1:n,])
            #Update vnovector too
            mergedRingPtsVnoVec <- c(PtToBeMergedVno, mergedRingPtsVnoVec)
          }
          
          
          
        }else if (nearestPtsIdx == n) {#The nearest pt is the last pt in mergedRing or
          
          preNearestIdx <- nearestPtsIdx - 1 
          preNearestVno <- mergedRingPtsMTX[nearestPtsIdx - 1, 3]
          
          #The next point of the last one is the first point!!!
          
          postNearestIdx <- 1
          postNearestVno <- mergedRingPtsMTX[1, 3]
          
          #Get the distances to see which way is the cheapest to connect!!!
          #Before the nearest pt or after the nearest pt!!!
          #Distance between pre and the nearest
          preDist <- theDistMTX[nearestPtsVno,  preNearestVno]
          #Distance between post and the nearest
          postDist <- theDistMTX[nearestPtsVno,  postNearestVno]
          #Distance between the point to be merged and the nearest
          nearest2ptDist <- theDistMTX[PtToBeMergedVno, nearestPtsVno]
          #Distance between the pre of the nearest and pt to be merged
          pre2ptDist <- theDistMTX[PtToBeMergedVno, preNearestVno]
          #Distance between the post of the nearest and pt to be merged
          post2ptDist <- theDistMTX[PtToBeMergedVno, postNearestVno]
          if ((preDist + post2ptDist) < (postDist + pre2ptDist)){
            #after the nearest
            #DEBUG
            #cat("3: nearest last after nearest\n")
            mergedRingPtsMTX <<- rbind(mergedRingPtsMTX[1:nearestPtsIdx,], PtToBeMerged)
            #Update vnovector too
            
            #PtToBeMerged will be the last element
            mergedRingPtsVnoVec <- c(mergedRingPtsVnoVec, PtToBeMergedVno)
            
          }else{
            #DEBUG
            #cat("4: nearest last before nearest\n")
            #before the nearest
            #update the last one since it will be the first point in the ring/polygon
            mergedRingPtsMTX <<- rbind(mergedRingPtsMTX[1:(nearestPtsIdx-1),], PtToBeMerged, mergedRingPtsMTX[nearestPtsIdx:n,])
            #Update vnovector too too
            mergedRingPtsVnoVec <- c(mergedRingPtsVnoVec[1:(nearestPtsIdx-1)], PtToBeMergedVno, mergedRingPtsVnoVec[nearestPtsIdx:n])
          }
          
          
          
          
          
        } else{ #The nearest (from the mergedRing) pt is in between
          #The nearest pt is the last pt in mergedRing or
          #The nearest (from the mergedRing) pt is in between
          preNearestIdx <- nearestPtsIdx - 1 
          preNearestVno <- mergedRingPtsMTX[nearestPtsIdx - 1, 3]
          
          
          postNearestIdx <- nearestPtsIdx + 1
          postNearestVno <- mergedRingPtsMTX[nearestPtsIdx + 1, 3]
          
          
          #Get the distances to see which way is the cheapest to connect!!!
          #Before the nearest pt or after the nearest pt!!!
          #Distance between pre and the nearest
          preDist <- theDistMTX[nearestPtsVno,  preNearestVno]
          #Distance between post and the nearest
          postDist <- theDistMTX[nearestPtsVno,  postNearestVno]
          #Distance between the point to be merged and the nearest
          nearest2ptDist <- theDistMTX[PtToBeMergedVno, nearestPtsVno]
          #Distance between the pre of the nearest and pt to be merged
          pre2ptDist <- theDistMTX[PtToBeMergedVno, preNearestVno]
          #Distance between the post of the nearest and pt to be merged
          post2ptDist <- theDistMTX[PtToBeMergedVno, postNearestVno]
          if ((preDist + post2ptDist) < (postDist + pre2ptDist)){
            #after the nearest
            #DEBUG
            #cat("5: nearest between after nearest\n")
            #IF THE NEAREST IS THE LAST THAN (nearestPtsIdx + 1) WILL GIVE ERROR!!!!
            mergedRingPtsMTX <<- rbind(mergedRingPtsMTX[1:nearestPtsIdx,], PtToBeMerged, mergedRingPtsMTX[(nearestPtsIdx + 1):n,])
            #Update vnovector too
            
            mergedRingPtsVnoVec <- c(mergedRingPtsVnoVec[1:nearestPtsIdx], PtToBeMergedVno, mergedRingPtsVnoVec[(nearestPtsIdx+1):n])
            
          }else{
            #DEBUG
            #cat("6: nearest between before nearest\n")
            #before the nearest
            #update the last one since it will be the first point in the ring/polygon
            mergedRingPtsMTX <<- rbind(mergedRingPtsMTX[1:(nearestPtsIdx-1),], PtToBeMerged, mergedRingPtsMTX[nearestPtsIdx:n,])
            #Update vnovector too too
            mergedRingPtsVnoVec <- c(mergedRingPtsVnoVec[1:(nearestPtsIdx-1)], PtToBeMergedVno, mergedRingPtsVnoVec[nearestPtsIdx:n])
          }
        } #else#else
        
        
        
        
        
        
        
        
        
        
        
        
        
        #############################################################################################################
        #DEBUG plot for merging
        #########################################################################################################
        #################################################################################################### 
        # mergedPts <- rbind(mergedPts, PtToBeMerged)
        # 
        # #For EPS output
        # cairo_ps(filename = (paste0(theTSP,"-merge-ring-", r, "-pt-", p,".eps")),
        #          width = 12,
        #          height = 7,
        #          pointsize = 12,
        #          fallback_resolution = 600)
        # 
        # #Draw the tsp
        # plot(theTSPObject,
        #      xlab = "x-coordinate",
        #      ylab = "y-coordinate",
        #      lty = 3, lwd=0.7, cex=1.5,  col = "black", pch= 19, asp=1,
        #      main = paste0(theTSP," --- Merging ring: ", r, " point: ", p),
        #      xlim=c(xmin, xmax),
        #      ylim=c(ymin, ymax)
        # )
        # 
        # #Specific vertex label positions for pr76
        # text(theTSPcoordMTX[c(1:71,73:75),1], theTSPcoordMTX[c(1:71,73:75),2], labels = c(1:71,73:75), cex = 0.8, col = "blue", pos = 3)
        # text(theTSPcoordMTX[c(72,76),1], theTSPcoordMTX[c(72,76),2], labels = c(72,76), cex = 0.8, col = "blue", pos = 4)
        # 
        # 
        # 
        # 
        # 
        # #mark next ring
        # if (nrow(ringToBeMergedPtsMTX ) > 0){
        #   points(ringToBeMergedPtsMTX[,1], ringToBeMergedPtsMTX[,2], col = "pink", pch = 19, cex=2)
        #   lines(ringToBeMergedPtsMTX[,1], ringToBeMergedPtsMTX[,2], col = "pink", pch = 19)
        #   lines(ringToBeMergedPtsMTX[c(1, nrow(ringToBeMergedPtsMTX)),1], 
        #         ringToBeMergedPtsMTX[c(1, nrow(ringToBeMergedPtsMTX)),2], col = "pink", pch = 19)
        #   text(ringToBeMergedPtsMTX[,1]+1, 
        #        ringToBeMergedPtsMTX[,2]+1,labels=1:nrow(ringToBeMergedPtsMTX), cex=0.75, font=2)
        # }
        # 
        # #Draw the merged ring 
        # colx =  "green" #rainbow(7)[sample(7,1)]
        # points(mergedRingPtsMTX[,1], mergedRingPtsMTX[,2], col = colx, pch = 19, cex=2)
        # lines(mergedRingPtsMTX[,1], mergedRingPtsMTX[,2], col = colx, pch = 19)
        # lines(mergedRingPtsMTX[c(1, nrow(mergedRingPtsMTX)),1], mergedRingPtsMTX[c(1, nrow(mergedRingPtsMTX)),2], col = colx, pch = 19)
        # text(mergedRingPtsMTX[,1]+1, mergedRingPtsMTX[,2]+1, labels=1:nrow(mergedRingPtsMTX), cex=0.75, font=2)
        # 
        # #mark the newly merged next ring vertices if any
        # if (nrow(mergedPts) >= 1){
        #   mergedPts <- matrix(mergedPts[2:nrow(mergedPts),], ncol=3) #remove the first row as it is NAs
        #   points(mergedPts[,1], mergedPts[,2], col = "blue", pch = 19, cex=2)
        # }
        # 
        # #mark the next point?
        # #PtToBeMerged <- st_point(ringToBeMergedPtsMTX[p,])
        # if (p < nPtsToBeMerged) {
        #   points(ringToBeMergedPtsMTX[p+1,1], ringToBeMergedPtsMTX[p+1,2], col="red", pch=19, cex=2)
        # }
        # 
        # dev.off()
        # 
        #############################################################################################################
        
        
        
        
        
        
        
      }#end for loop over ring points
    }#end for loop over rings
    
  }else {
    #Single ring
    cat("Single ring! Checking if there are pts remained.\n")
    
  }
  
  ####################################################################
  #DO NOT FORGET TO ADD THEM!!!
  #ptsNotInAnyRingMTX 
  ####################################################################
  
  lenRemained <- nrow(ptsNotInAnyRingMTX)
  
  if (lenRemained > 0) { 
    cat(lenRemained,"pts remained that are not in any ring.\n")
    ringToBeMergedPtsMTX <- ptsNotInAnyRingMTX
    nPtsToBeMerged <- nrow(ringToBeMergedPtsMTX)
    
    for (k in 1:lenRemained) {
      
      #DEBUG
      #k <- 1
      
      
      ##################################################################################  
      #Merging criteria: Closest point from the MergedRing
      ##################################################################################
      #Try to find the nearest pt from the merged ring
      #ATTENTION: The merged ring changes after each point addition
      #Number of distinct points in the ring,
      #excluding the last one since it is same with the first one (polygons and rings are closed)
      n <- nrow(mergedRingPtsMTX)
      
      
      PtToBeMerged <- ringToBeMergedPtsMTX[k,]
      PtToBeMergedVno <- ringToBeMergedPtsMTX[k,3]
      
      #The min pt should also be on the merged ring
      minPtidx <- which.min(theDistMTX[PtToBeMergedVno, mergedRingPtsVnoVec])
      minPtVno <- mergedRingPtsVnoVec[minPtidx]
      
      nearestPtsIdx <- minPtidx
      nearestPtsVno <- minPtVno
      
      
      
      ##################################################################################
      #Now we find the nearest pt on the mergedRing to connect the pt.
      
      
      if (nearestPtsIdx == 1) {
        #The nearest pt is the first pt in mergedRing
        #The previous point is the penultimate point of the merged ring
        #Attention the last one and the first one is same!!!
        preNearestIdx <- n
        preNearestVno <- mergedRingPtsMTX[n, 3]
        postNearestIdx <- 2 
        postNearestVno <- mergedRingPtsMTX[2, 3]
        
        #Get the distances to see which way is the cheapest to connect!!!
        #Before the nearest pt or after the nearest pt!!!
        #Distance between pre and the nearest
        preDist <- theDistMTX[nearestPtsVno,  preNearestVno]
        #Distance between post and the nearest
        postDist <- theDistMTX[nearestPtsVno,  postNearestVno]
        #Distance between the point to be merged and the nearest
        nearest2ptDist <- theDistMTX[PtToBeMergedVno, nearestPtsVno]
        #Distance between the pre of the nearest and pt to be merged
        pre2ptDist <- theDistMTX[PtToBeMergedVno, preNearestVno]
        #Distance between the post of the nearest and pt to be merged
        post2ptDist <- theDistMTX[PtToBeMergedVno, postNearestVno]
        
        if ((preDist + post2ptDist) < (postDist + pre2ptDist)){
          #after the nearest
          #DEBUG
          #cat("1: nearest first after nearest\n")
          mergedRingPtsMTX <<- rbind(mergedRingPtsMTX[1:nearestPtsIdx,], PtToBeMerged, mergedRingPtsMTX[(nearestPtsIdx + 1):n,] )
          #Update vnovector too
          mergedRingPtsVnoVec <- c(mergedRingPtsVnoVec[1:nearestPtsIdx], PtToBeMergedVno, mergedRingPtsVnoVec[(nearestPtsIdx + 1):n])
        }else{
          #before the nearest
          #DEBUG
          #cat("2: nearest first before nearest\n")
          mergedRingPtsMTX <<- rbind(PtToBeMerged, mergedRingPtsMTX[1:n,])
          #Update vnovector too
          mergedRingPtsVnoVec <- c(PtToBeMergedVno, mergedRingPtsVnoVec)
        }
        
        
        
      }else if (nearestPtsIdx == n) {#The nearest pt is the last pt in mergedRing or
        
        preNearestIdx <- nearestPtsIdx - 1 
        preNearestVno <- mergedRingPtsMTX[nearestPtsIdx - 1, 3]
        
        #The next point of the last one is the first point!!!
        
        postNearestIdx <- 1
        postNearestVno <- mergedRingPtsMTX[1, 3]
        
        #Get the distances to see which way is the cheapest to connect!!!
        #Before the nearest pt or after the nearest pt!!!
        #Distance between pre and the nearest
        preDist <- theDistMTX[nearestPtsVno,  preNearestVno]
        #Distance between post and the nearest
        postDist <- theDistMTX[nearestPtsVno,  postNearestVno]
        #Distance between the point to be merged and the nearest
        nearest2ptDist <- theDistMTX[PtToBeMergedVno, nearestPtsVno]
        #Distance between the pre of the nearest and pt to be merged
        pre2ptDist <- theDistMTX[PtToBeMergedVno, preNearestVno]
        #Distance between the post of the nearest and pt to be merged
        post2ptDist <- theDistMTX[PtToBeMergedVno, postNearestVno]
        if ((preDist + post2ptDist) < (postDist + pre2ptDist)){
          #after the nearest
          #DEBUG
          #cat("3: nearest last after nearest\n")
          mergedRingPtsMTX <<- rbind(mergedRingPtsMTX[1:nearestPtsIdx,], PtToBeMerged)
          #Update vnovector too
          
          #PtToBeMerged will be the last element
          mergedRingPtsVnoVec <- c(mergedRingPtsVnoVec, PtToBeMergedVno)
          
        }else{
          #DEBUG
          #cat("4: nearest last before nearest\n")
          #before the nearest
          #update the last one since it will be the first point in the ring/polygon
          mergedRingPtsMTX <<- rbind(mergedRingPtsMTX[1:(nearestPtsIdx-1),], PtToBeMerged, mergedRingPtsMTX[nearestPtsIdx:n,])
          #Update vnovector too too
          mergedRingPtsVnoVec <- c(mergedRingPtsVnoVec[1:(nearestPtsIdx-1)], PtToBeMergedVno, mergedRingPtsVnoVec[nearestPtsIdx:n])
        }
        
        
        
        
        
      } else{ #The nearest (from the mergedRing) pt is in between
        #The nearest pt is the last pt in mergedRing or
        #The nearest (from the mergedRing) pt is in between
        preNearestIdx <- nearestPtsIdx - 1 
        preNearestVno <- mergedRingPtsMTX[nearestPtsIdx - 1, 3]
        
        
        postNearestIdx <- nearestPtsIdx + 1
        postNearestVno <- mergedRingPtsMTX[nearestPtsIdx + 1, 3]
        
        
        #Get the distances to see which way is the cheapest to connect!!!
        #Before the nearest pt or after the nearest pt!!!
        #Distance between pre and the nearest
        preDist <- theDistMTX[nearestPtsVno,  preNearestVno]
        #Distance between post and the nearest
        postDist <- theDistMTX[nearestPtsVno,  postNearestVno]
        #Distance between the point to be merged and the nearest
        nearest2ptDist <- theDistMTX[PtToBeMergedVno, nearestPtsVno]
        #Distance between the pre of the nearest and pt to be merged
        pre2ptDist <- theDistMTX[PtToBeMergedVno, preNearestVno]
        #Distance between the post of the nearest and pt to be merged
        post2ptDist <- theDistMTX[PtToBeMergedVno, postNearestVno]
        if ((preDist + post2ptDist) < (postDist + pre2ptDist)){
          #after the nearest
          #DEBUG
          #cat("5: nearest between after nearest\n")
          #IF THE NEAREST IS THE LAST THAN (nearestPtsIdx + 1) WILL GIVE ERROR!!!!
          mergedRingPtsMTX <<- rbind(mergedRingPtsMTX[1:nearestPtsIdx,], PtToBeMerged, mergedRingPtsMTX[(nearestPtsIdx + 1):n,])
          #Update vnovector too
          
          mergedRingPtsVnoVec <- c(mergedRingPtsVnoVec[1:nearestPtsIdx], PtToBeMergedVno, mergedRingPtsVnoVec[(nearestPtsIdx+1):n])
          
        }else{
          #DEBUG
          #cat("6: nearest between before nearest\n")
          #before the nearest
          #update the last one since it will be the first point in the ring/polygon
          mergedRingPtsMTX <<- rbind(mergedRingPtsMTX[1:(nearestPtsIdx-1),], PtToBeMerged, mergedRingPtsMTX[nearestPtsIdx:n,])
          #Update vnovector too too
          mergedRingPtsVnoVec <- c(mergedRingPtsVnoVec[1:(nearestPtsIdx-1)], PtToBeMergedVno, mergedRingPtsVnoVec[nearestPtsIdx:n])
        }
      } #else#else
      
      
      
      
      
      
      
      
    }
    
  }else{
    cat("All pts were in rings.\n")
  }
  
  
  
  #########################################################################################################
  #END MY ALGO
  #########################################################################################################
  
  toc(log = TRUE, quiet = TRUE)
  
  
  
  #############################################################################################################
  #DEBUG plot for merged ring
  #########################################################################################################
  #################################################################################################### 
  
  
  # #For EPS output
  # cairo_ps(filename = (paste0(theTSP,"-final-merged-ring.eps")),
  #          width = 12,
  #          height = 7,
  #          pointsize = 12,
  #          fallback_resolution = 600)
  # 
  # #Draw the tsp
  # plot(theTSPObject, theOptimTour, 
  #      xlab = "x-coordinate",
  #      ylab = "y-coordinate",
  #      lty = 3, lwd=0.7, cex=1.5,  col = "black", pch= 19, asp=1,
  #      main = paste0(theTSP," --- Final merged ring"),
  #      xlim=c(xmin, xmax),
  #      ylim=c(ymin, ymax)
  # )
  # 
  # #Specific vertex label positions for pr76
  # text(theTSPcoordMTX[c(1:71,73:75),1], theTSPcoordMTX[c(1:71,73:75),2], labels = c(1:71,73:75), cex = 0.8, col = "red", pos = 3)
  # text(theTSPcoordMTX[c(72,76),1], theTSPcoordMTX[c(72,76),2], labels = c(72,76), cex = 0.8, col = "red", pos = 4)
  # 
  # 
  # 
  # 
  # #Draw the merged ring 
  # colx = myGreenTrans
  # points(mergedRingPtsMTX[,1], mergedRingPtsMTX[,2], col = "green", pch = 19, cex=2)
  # lines(mergedRingPtsMTX[,1], mergedRingPtsMTX[,2],  lwd=2, col = colx, pch = 19)
  # lines(mergedRingPtsMTX[c(1, nrow(mergedRingPtsMTX)),1], mergedRingPtsMTX[c(1, nrow(mergedRingPtsMTX)),2], lwd=2, col = colx, pch = 19)
  # text(mergedRingPtsMTX[,1]+1, mergedRingPtsMTX[,2]+1, labels=1:nrow(mergedRingPtsMTX), cex=0.75, font=2)
  # 
  # 
  # 
  # 
  # 
  # dev.off()
  # 
  
  
  #############################################################################################################
  #DEBUG plot for merged ring
  #########################################################################################################
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ############################### GET TOUR COST AND APPROXIMATION RATIO ############################### 
  #Add the first point to the end to complete the cycle!!!!
  m2 <- rbind(mergedRingPtsMTX[2:nrow(mergedRingPtsMTX),], mergedRingPtsMTX[1,])
  m3 <- mergedRingPtsMTX[1:nrow(mergedRingPtsMTX),]
  
  #Preventing NAs and integer overflow
  storage.mode(m2) <- "numeric"
  storage.mode(m3) <- "numeric"
  
  
  distMTX <- as.double(sqrt((m3[,1] - m2[,1]) * (m3[,1] - m2[,1]) + 
                              (m3[,2] - m2[,2]) * (m3[,2] - m2[,2])))
  concaveRingMergeNearestPtV2TourCost <- sum(distMTX)
  
  cat("The tour has", nrow(mergedRingPtsMTX), "vertices with total cost:", concaveRingMergeNearestPtV2TourCost ,"\n")
  
  
  if (theOptimTourCost == -1){
    cat("The concaveRingMergeNearestPtV2 with cr =", cr,"\n")
    aRatioData <- c(aRatioData, NA)
    
  }else{
    concaveRingMergeNearestPtV2Ratio <- round(concaveRingMergeNearestPtV2TourCost/theOptimTourCost, 3)
    cat("The approximation ratio of concaveRingMergeNearestPtV2 with cr =", cr,
        "is: (",concaveRingMergeNearestPtV2TourCost,"/", theOptimTourCost,") =",
        concaveRingMergeNearestPtV2Ratio,"\n")
    aRatioData <- c(aRatioData, concaveRingMergeNearestPtV2Ratio)
  }
  
  tourCostData <- c(tourCostData, concaveRingMergeNearestPtV2TourCost)
  
  
  #Lets find the vertex ids on the tour
  concaveRingMergeNearestPtV2TourLength <-  nrow(mergedRingPtsMTX)
  vtxIds <- vector()
  #Rounding is necessary for comparison
  theTSPcoordDF <- round(theTSPcoordDF, coordDigitPrec)
  mergedRingPtsMTX <- round(mergedRingPtsMTX, coordDigitPrec)
  
  for (i in 1:concaveRingMergeNearestPtV2TourLength) {
    idx <- which(theTSPcoordDF$x ==  mergedRingPtsMTX[i,1] & theTSPcoordDF$y ==  mergedRingPtsMTX[i,2])
    
    if (length(is.na(idx)) == 0) {
      cat("ERROR: The point", i ,"(", mergedRingPtsMTX[i,1],  mergedRingPtsMTX[i,2],") was not found in coords!!!\n")
    }else{
      vtxIds[i] <- idx 
    }
  }
  
  #Also
  #my_get_cost_perm(mergedRingPtsVnoVec)
  
  concaveRingMergeNearestPtV2Tour <- vtxIds
  if ((round(concaveRingMergeNearestPtV2TourCost)) != (round(-1 * my_get_cost_perm( concaveRingMergeNearestPtV2Tour)))){
    cat("ERROR:((round(concaveRingMergeNearestPtV2TourCost)) != (round(-1 * my_get_cost_perm( concaveRingMergeNearestPtV2Tour))))\n")
  }

  
  
  
  concaveRingMergeNearestPtV2AWT <- my_get_awt_from_vtx1(concaveRingMergeNearestPtV2Tour)
  awdData <- c(awdData, concaveRingMergeNearestPtV2AWT)
  cat("concaveRingMergeNearestPtV2 has awt of (From vtx no 1):", concaveRingMergeNearestPtV2AWT,"\n")
  
  
  
  if (minAWDFlag == 1){
    concaveRingMergeNearestPtV2MinAWD <- my_get_min_awt_perm_fast(concaveRingMergeNearestPtV2Tour)
    minAWDData <- c(minAWDData, concaveRingMergeNearestPtV2MinAWD[[1]])
  }else{
    minAWDData <- c(minAWDData, NA)
  }
  
  #Sometimes DT fails!!!
  #Just skip it
  #Uncomment the following then
  #delEdgePercentData <- c(delEdgePercentData, NA)
  
  p <- my_get_del_percent_fast(concaveRingMergeNearestPtV2Tour, edgesDelaunayDF)
  cat("The Delaunay percentage is:",p,"\n")
  delEdgePercentData <- c(delEdgePercentData, p)
  
  
  ############################### GET TOUR COST AND APPROXIMATION RATIO ###############################
  cat("cr =",cr,"--",nRings, "rings generated.\n")
  
  
  rSum <- my_sum_cost_rings(ringMTXList)
  cat("The sum of", nRings, "rings costs:", rSum, "the final merged ring cost:",concaveRingMergeNearestPtV2TourCost,"\n")
  
  
  cat("---------------------------------------------------------------------------------------------\n")
  
}


log.txt <- tic.log(format = TRUE)
log.lst <- tic.log(format = FALSE)
tic.clearlog()
timings <- unlist(lapply(log.lst, function(x) x$toc - x$tic))
writeLines(unlist(log.txt))


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
  theIdxCol <- c(theIdxCol, "OPT")
  algoNameData <- c(algoNameData, theMethod)
  timingData <- c(timingData, 0)
  aRatioData <- c(aRatioData, 1)
  tourCostData <- c(tourCostData, theOptimTourCost)
  
  #Some data set they have the cost data for optimtour but not the tour itself
  if (is.null(theOptimTour) ){
    awdData  <- c(awdData, NA)
    minAWDData  <- c(minAWDData, NA)
    delEdgePercentData <- c(delEdgePercentData, NA)
  }else {
    awdData  <- c(awdData, my_get_awt_perm(theOptimTour))
    
    if (minAWDFlag == 1){
      minOptimAWD <- my_get_min_awt_perm_fast(theOptimTour)
      minAWDData  <- c(minAWDData, minOptimAWD[[1]])
    }else{
      minAWDData <- c(minAWDData, NA)
    }
    
    
    
    
    delEdgePercentData <- c(delEdgePercentData, my_get_del_percent_fast(theOptimTour, edgesDelaunayDF))
  }
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


cat(DistMTXSecs,"secs for dmtx creation.\n")
cat("cr =",cr,"--",nRings, "rings generated.\n")
rSum <- my_sum_cost_rings(ringMTXList)
cat("The sum of", nRings, "rings costs:", rSum, "the final merged ring cost:",concaveRingMergeNearestPtV2TourCost,"\n")

write_excel_csv(resultsDF,outputFile)
cat("--------------------------------------------------------------------------------------------\n")




# 
# 
# #############################################################################################################
# #Uncomment the following to see the plot of the optimtour
# #############################################################################################################
# #DEBUG
# xmin <- min(theTSPcoordMTX[,1])
# xmax <- max(theTSPcoordMTX[,1])
# 
# ymin <- min(theTSPcoordMTX[,2])
# ymax <- max(theTSPcoordMTX[,2])
# 
# 
# myWidth <- 800
# myHeight <- 600
# 
# #Draw the optimal tour
# plot(mergedRingPtsMTX,
#      xlab = "x-coordinate", 
#      ylab = "y-coordinate", 
#      lty = 3, lwd=0.7, cex=1,  col = "black", pch= 4, 
#      main = paste0(theTSP," - ",theMethod," - cost=",round(concaveRingMergeNearestPtV3TourCost,3)," - cr=",cr),
#      xlim=c(xmin, xmax), 
#      ylim=c(ymin, ymax)
# )
# 
# 
# #Draw the merged ring 
# colx =  "green"
# points(mergedRingPtsMTX[,1], mergedRingPtsMTX[,2], col = colx, pch = 19, cex=2)
# lines(mergedRingPtsMTX[,1], mergedRingPtsMTX[,2],  lwd=2.0, col = colx, pch = 19)
# lines(mergedRingPtsMTX[c(1, nrow(mergedRingPtsMTX)),1], mergedRingPtsMTX[c(1, nrow(mergedRingPtsMTX)),2], lwd=2.0, col = colx, pch = 19)
# text(mergedRingPtsMTX[,1]+1, mergedRingPtsMTX[,2]+1, labels=1:nrow(mergedRingPtsMTX), cex=0.75, font=2)
# #############################################################################################################
# 




