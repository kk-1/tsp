
#######################################################################################################
#Clean the environment
remove(list = ls())
########################################################################################################

myWorkDir <- "/my_hd1/my_dir/my_prg/Camerino/tsp" #The work dir
myResultsDir <- "/my_hd1/my_dir/my_prg/Camerino/tsp/results" #Results dir


setwd(myWorkDir)
########################################################################################################

library(tidyverse)
library(ggplot2)
library(reshape2)
########################################################################################################




#############################################################################
# Just set these parallel arrays below and good to go 
#############################################################################

#List algorithms 
Algos <- c(
  "concaveRingMergeNearestPtV2",
  "concaveRingMergeNearestPtV3",
  "nearest_insertion",
  "farthest_insertion",
  "cheapest_insertion",
  "arbitrary_insertion",
  "nn",
  "repetitive_nn",
  "two_opt")


#Header for algo names  
myAlgoNames <- c("CH" , "CNH",  "NI" , "FI", "CI" , "AI" , "NN" , "RNN", "2-Opt")



#Result files for the selected algos
#Incase there are garbage csv files list the actual result files to be processed
csvFileNames <- c(
  "concaveRingMergeNearestPtV2-3-cr-0.9.csv",
  "concaveRingMergeNearestPtV3-3-cr-0.9.csv",
  "nearest_insertion-20.csv",
  "farthest_insertion-20.csv",
  "cheapest_insertion-20.csv",
  "arbitrary_insertion-20.csv",
  "nn-20.csv",
  "repetitive_nn-20.csv",
  "two_opt-20.csv")


#############################################################################
#Big datasets
#############################################################################


#Set it to 1 if there is a tour file 0 ow
dataSetsOptimTour <- c(
1,  # "pr1002",
1,  # "pr2392",
0,  # "rl5915",
0,  # "usa13509",
1,  # "pla33810E",
0,  # "myLattice-100x100-10000",
0,  # "myHexLattice-100x100-10000",
0,  # "myRNDHexLattice-105x105-10000",
0,  # "myRNDLattice-105x105-10000",
0,  # "myLattice-100x200-20000",
0,  # "myHexLattice-100x200-20000",
0,  # "myRNDHexLattice-105x210-20000",
0,  # "myRNDLattice-105x210-20000",
0,  # "myLattice-150x200-30000",
0,  # "myHexLattice-150x200-30000",
0,  # "myRNDLattice-158x210-30000",
0  # "myRNDHexLattice-158x210-30000"
)


#List big datasets
dataSets <- c(
"pr1002",
"pr2392",
"rl5915",
"usa13509",
"pla33810E",
"myLattice-100x100-10000",
"myHexLattice-100x100-10000",
"myRNDHexLattice-105x105-10000",
"myRNDLattice-105x105-10000",
"myLattice-100x200-20000",
"myHexLattice-100x200-20000",
"myRNDHexLattice-105x210-20000",
"myRNDLattice-105x210-20000",
"myLattice-150x200-30000",
"myHexLattice-150x200-30000",
"myRNDLattice-158x210-30000",
"myRNDHexLattice-158x210-30000"
)

#List their optim tour costs if they are known -1 ow. FROM THE LOG FILE OR FROM CSV FILES!!!
dataSetsOptimCost <- c(
259066,  # "pr1002",
378062,  # "pr2392",
565530,  # "rl5915",
19982859,  # "usa13509",
66043099,  # "pla33810E",
-1,  # "myLattice-100x100-10000",
-1,  # "myHexLattice-100x100-10000",
-1,  # "myRNDHexLattice-105x105-10000",
-1,  # "myRNDLattice-105x105-10000",
-1,  # "myLattice-100x200-20000",
-1,  # "myHexLattice-100x200-20000",
-1,  # "myRNDHexLattice-105x210-20000",
-1,  # "myRNDLattice-105x210-20000",
-1,  # "myLattice-150x200-30000",
-1,  # "myHexLattice-150x200-30000",
-1,  # "myRNDLattice-158x210-30000",
-1   # "myRNDHexLattice-158x210-30000"
)





# #############################################################################
# #Small datasets
# #############################################################################
# #Check if the tsplib dataset has opt.tour or not
# dataSetsOptimTour <- c(
#   1,   #   "att48",
#   1,   #   "berlin52",
#   1,   #   "pr76",
#   1,   #   "kroA100",
#   1,   #   "lin105",
#   1,   #   "ch130",
#   1,   #   "ch150",
#   1,   #   "a280",
#   1,   #   "pcb442",
#   0,   #   "Tnm52",
#   0,   #   "Tnm76",
#   0,   #   "Tnm100",
#   0,   #   "Tnm127",
#   0,   #   "Tnm154",
#   0,   #   "Tnm178",
#   0,   #   "Tnm199",
#   0,   #   "myRND-100",
#   0,   #   "myRND-200",
#   0,   #   "myRND-300",
#   0,   #   "myRND-400",
#   0,   #   "myLattice-10x10-100",
#   0,   #   "myLattice-10x20-200",
#   0,   #   "myLattice-15x20-300",
#   0,   #   "myLattice-20x20-400",
#   0,   #   "myRNDLattice-12x12-100"
#   0,   #   "myRNDLattice-12x23-200"
#   0,   #   "myRNDLattice-18x23-300"
#   0,   #   "myRNDLattice-23x23-400"
#   0,   #   "myHexLattice-10x10-100",
#   0,   #   "myHexLattice-10x20-200",
#   0,   #   "myHexLattice-15x20-300",
#   0,   #   "myHexLattice-20x20-400",
#   0,   #   "myRNDHexLattice-12x12-100"
#   0,   #   "myRNDHexLattice-12x23-200"
#   0,   #   "myRNDHexLattice-18x23-300"
#   0    #   "myRNDHexLattice-23x23-400"
# )
# 
# 
# #List datasets
# dataSets <- c(
#   "att48",
#   "berlin52",
#   "pr76",
#   "kroA100",
#   "lin105",
#   "ch130",
#   "ch150",
#   "a280",
#   "pcb442",
#   "Tnm52",
#   "Tnm76",
#   "Tnm100",
#   "Tnm127",
#   "Tnm154",
#   "Tnm178",
#   "Tnm199",
#   "myRND-100",
#   "myRND-200",
#   "myRND-300",
#   "myRND-400",
#   "myLattice-10x10-100",
#   "myLattice-10x20-200",
#   "myLattice-15x20-300",
#   "myLattice-20x20-400",
#   "myRNDLattice-12x12-100",
#   "myRNDLattice-12x23-200",
#   "myRNDLattice-18x23-300",
#   "myRNDLattice-23x23-400",
#   "myHexLattice-10x10-100",
#   "myHexLattice-10x20-200",
#   "myHexLattice-15x20-300",
#   "myHexLattice-20x20-400",
#   "myRNDHexLattice-12x12-100",
#   "myRNDHexLattice-12x23-200",
#   "myRNDHexLattice-18x23-300",
#   "myRNDHexLattice-23x23-400"
# )
# 
# 
# # 259066.663,  #   "pr1002",
# 
# #List their optim tour costs if they are known
# dataSetsOptimCost <- c(
#   33523.709,   #   "att48",
#   7544.366,    #   "berlin52",
#   108159.438,  #   "pr76",
#   21285.443,   #   "kroA100",
#   14382.996,   #   "lin105",
#   6110.861,    #   "ch130",
#   6532.281,    #   "ch150",
#   2586.77,     #   "a280",
#   50783.548,   #   "pcb442",
#   551609,      #   "Tnm52",
#   949961,      #   "Tnm76",
#   1398070,     #   "Tnm100",
#   1871162,     #   "Tnm127",
#   2350345,     #   "Tnm154",
#   2771953,     #   "Tnm178",
#   3139778,     #   "Tnm199",
#   -1,   #   "myRND-100",
#   -1,   #   "myRND-200",
#   -1,   #   "myRND-300",
#   -1,   #   "myRND-400",
#   -1,   #   "myLattice-10x10-100",
#   -1,   #   "myLattice-10x20-200",
#   -1,   #   "myLattice-15x20-300",
#   -1,   #   "myLattice-20x20-400",
#   -1,   #   "myRNDLattice-12x12-100"
#   -1,   #   "myRNDLattice-12x23-200"
#   -1,   #   "myRNDLattice-18x23-300"
#   -1,   #   "myRNDLattice-23x23-400"
#   -1,   #   "myHexLattice-10x10-100",
#   -1,   #   "myHexLattice-10x20-200",
#   -1,   #   "myHexLattice-15x20-300",
#   -1,   #   "myHexLattice-20x20-400",
#   -1,   #   "myRNDHexLattice-12x12-100"
#   -1,   #   "myRNDHexLattice-12x23-200"
#   -1,   #   "myRNDHexLattice-18x23-300"
#   -1    #   "myRNDHexLattice-23x23-400"
# )
# #############################################################################



########################################################################################################

#Read files

cat("These should be equal to eachother:",length(dataSets),"data sets", length(dataSetsOptimTour),"optim tour data", length(dataSetsOptimCost),"optim tour cost data!\n")

dataSetIDX <- 1
for (d in dataSets) {
  
  # #DEBUG
  #d <- "att48"; dataSetIDX <- 1
  #d <- "berlin52"; dataSetIDX <- 2
  # d <-  "Tnm52"; dataSetIDX <- 6
  # d <-  "myHexLattice-20x30"; dataSetIDX <- 19
  
  
 
  
  
  csvDir <- paste(myResultsDir, d, sep = "/")
  
  cat("Processing for",d, "data in:", csvDir,"\n")
  
  csvFiles <-  paste(csvDir, csvFileNames , sep = "/")
  dataSetDF <- csvFiles %>% map_df(~read_csv(.))
  avgDF <- dataSetDF[which(dataSetDF$ResultId == "AVG"),]
  if (dataSetsOptimCost[dataSetIDX] != -1){
    avgDF2 <- avgDF %>% select(1,3:8)
  }else{
    avgDF2 <- avgDF %>% select(1,3:7)
  }
  
  if (dataSetIDX == 1){
    #For the first data set populate the mega dataframe 
    dsCol <- rep(d,nrow(avgDF))
    dsColDF <- as.data.frame(dsCol)
    colnames(dsColDF) <- c("DataSet")
    allDF <- cbind(dsColDF,avgDF)
    
    #The row for the optimum tour
    if (dataSetsOptimCost[dataSetIDX] == -1){
      #OPT tour cost is unknown so tour is also not known
      optRow <- c(d,"OPT","OPT", NA,NA, -1,NA,NA,NA)
      allDF <- rbind(allDF,optRow)
    }else{
      #OPT tour cost is known, check if tour is known
      if (dataSetsOptimTour[dataSetIDX] == 1){
        #The optim tour is known so get metrics for AWD and DelEdgeP
        #WE ASSSUME THAT ALL DATASETS WITH KNOWN TOURS GIVE DEL EDGES!!!!
        optimAWDfromVTX1 <- dataSetDF[which(dataSetDF$ResultId == "OPT"),]$AWDfromVTX1[1]
        optimMinAWD <- dataSetDF[which(dataSetDF$ResultId == "OPT"),]$MinAWD[1]
        optimDelEdgeP <- dataSetDF[which(dataSetDF$ResultId == "OPT"),]$DelEdgeP[1]
        
        optRow <- c(d,"OPT","OPT", NA,1,dataSetsOptimCost[dataSetIDX], optimAWDfromVTX1, optimMinAWD, optimDelEdgeP)
      }else{
        #The optim tour is not known so get metrics only for the cost
        optRow <- c(d,"OPT","OPT", NA,1,dataSetsOptimCost[dataSetIDX], NA, NA, NA)
      }
      allDF <- rbind(allDF,optRow)
    }
    
    
    
    
  }else{
    #For the datasets other than the first one just append their data to mega DF
    dsCol <- rep(d,nrow(avgDF))
    dsColDF <- as.data.frame(dsCol)
    colnames(dsColDF) <- c("DataSet")
    tempDF <- cbind(dsColDF,avgDF)
    
    
    #The row for the optimum tour
    if (dataSetsOptimCost[dataSetIDX] == -1){
      
      #Non-standard row without approx ratio there is problem to bind
      #Insert approxratio col with NAa
      arCol <- rep(NA,nrow(tempDF))
      arColDF <- as.data.frame(arCol)
      colnames(arColDF) <- c("ApproxRatio")
      tempDF <- cbind((tempDF %>% select(1:4)), arColDF, (tempDF %>% select(5:8)) )
      
      allDF <- rbind(allDF, tempDF)
      
      #OPT tour cost is unknown so tour is also not known
      optRow <- c(d,"OPT","OPT", NA,NA, -1,NA,NA,NA)
      allDF <- rbind(allDF,optRow)
    }else{
      #Standard row with approx ratio no problem to bind
      allDF <- rbind(allDF, tempDF)
      #OPT tour cost is known, check if tour is known
      if (dataSetsOptimTour[dataSetIDX] == 1){
        #The optim tour is known so get metrics for AWD and DelEdgeP
        #WE ASSSUME THAT ALL DATASETS WITH KNOWN TOURS GIVE DEL EDGES!!!!
        optimAWDfromVTX1 <- dataSetDF[which(dataSetDF$ResultId == "OPT"),]$AWDfromVTX1[1]
        optimMinAWD <- dataSetDF[which(dataSetDF$ResultId == "OPT"),]$MinAWD[1]
        optimDelEdgeP <- dataSetDF[which(dataSetDF$ResultId == "OPT"),]$DelEdgeP[1]
        
        optRow <- c(d,"OPT","OPT", NA,1,dataSetsOptimCost[dataSetIDX], optimAWDfromVTX1, optimMinAWD, optimDelEdgeP)
      }else{
        #The optim tour is not known so get metrics only for the cost
        optRow <- c(d,"OPT","OPT", NA,1,dataSetsOptimCost[dataSetIDX], NA, NA, NA)
      }
      allDF <- rbind(allDF,optRow)
    }
    
    
    
    
  }
  
  
  
  
  dataSetIDX <- dataSetIDX + 1
  
  
}
########################################################################################################

cat("Results in", myResultsDir ,"are tabulated!\n")


########################################################################################################

#Now give the latex tables for each measure in allDF
# 
# allDF
# DataSet                        Algo ResultId Time(sec) ApproxRatio    TourCost AWDfromVTX1     MinAWD DelEdgeP
# 1      att48 concaveRingMergeNearestPtV2      AVG     0.018       1.062   35617.959   21415.965  14944.035    93.75
# 2      att48 concaveRingMergeNearestPtV3      AVG     0.019       1.068   35811.198   21470.195  14966.535    93.75
# 3      att48           nearest_insertion      AVG     0.003       1.123   37636.194    20027.21  16191.591   88.542
# 4      att48          farthest_insertion      AVG     0.003       1.053   35307.371   17867.986  14760.577    96.25
# 5      att48          cheapest_insertion      AVG     0.002       1.083   36305.647   18254.169  15441.955   92.604
# 6      att48         arbitrary_insertion      AVG         0       1.057    35439.76   17350.535  14775.081       95
# 7      att48                          nn      AVG         0       1.246    41786.29   20675.131  15733.505   88.125
# 8      att48               repetitive_nn      AVG     0.019        1.17   39236.885   22960.855  15394.334   89.583
# 9      att48                     two_opt      AVG         0       1.109   37176.156   19006.601  14885.482   96.458
# 10     att48                         OPT      OPT      <NA>           1   33523.709   15623.104  14164.066   97.917
# 11  berlin52 concaveRingMergeNearestPtV2      AVG     0.026       1.195    9013.672    4111.265   3370.176   86.538
# 12  berlin52 concaveRingMergeNearestPtV3      AVG     0.024       1.195    9013.672    4111.265   3370.176   86.538
# 13  berlin52           nearest_insertion      AVG     0.004       1.212    9146.656    4616.622   3880.166   88.942
# 14  berlin52          farthest_insertion      AVG     0.004       1.075    8112.208     4166.79   3105.281   95.769
# 15  berlin52          cheapest_insertion      AVG     0.003       1.191    8987.814    4589.605   3755.396   89.712
# 16  berlin52         arbitrary_insertion      AVG     0.001       1.131     8528.92    4317.674    3254.78     92.5
# 17  berlin52                          nn      AVG     0.001       1.243    9374.312    4642.834   3071.349   91.154
# 18  berlin52               repetitive_nn      AVG     0.022       1.085    8182.192    3961.128   3390.523   98.077
# 19  berlin52                     two_opt      AVG         0       1.132    8542.678    4307.289   3216.422   96.442
# 20  berlin52                         OPT      OPT      <NA>           1    7544.366    4131.786   2886.312   98.077



########################################################################################################
#We got 6 metrics: Time(sec) ApproxRatio    TourCost AWDfromVTX1     MinAWD DelEdgeP
#Latex table for each metrics
#Datasets on the rows and algos for columns
#MArk on each row the best and the worst and count them at the end
########################################################################################################


allDFAVG <- allDF[which(allDF$ResultId == "AVG"),]
allDFOPT <- allDF[which(allDF$ResultId == "OPT"),]


########################################################################################################
#Time
########################################################################################################
timeDF <- as.data.frame(dataSets)
colnames(timeDF) <- c("DataSet")

for (a in Algos) {
  #DEBUG
  #a <- Algos[1]
  
  timeCol <- allDFAVG[which(allDFAVG$Algo == a),]$`Time(sec)`
  timeCol <- as.numeric(timeCol)
  timeColDF <- as.data.frame(as.numeric(timeCol))
  #colnames(timeColDF) <- c(a)
  colnames(timeColDF) <- myAlgoNames[which(Algos == a)]
  
  timeDF <- cbind(timeDF, timeColDF)
}

#No need to add OPT col
#Print results into latex file 

results <- xtable::xtable(timeDF, 
                          caption ="Benchmark results for running time in seconds.",
                          label = "runtime",
                          digits=c(0,0,rep(3,ncol(timeDF)-1))
)

xtable::print.xtable(results, 
                     type = "latex", 
                     file = paste(myResultsDir,"time-results.txt", sep = "/"), 
                     caption.placement = "top",
                     include.rownames=FALSE
)



########################################################################################################
#ApproxRatio 
########################################################################################################

arDF <- as.data.frame(dataSets)
colnames(arDF) <- c("DataSet")

for (a in Algos) {
  #DEBUG
  #a <- Algos[1]
  
  arCol <- allDFAVG[which(allDFAVG$Algo == a),]$ApproxRatio
  ar <- as.numeric(arCol)
  arColDF <- as.data.frame(as.numeric(arCol))
  #colnames(timeColDF) <- c(a)
  colnames(arColDF) <- myAlgoNames[which(Algos == a)]
  
  arDF <- cbind(arDF, arColDF)
}

#No need to add OPT col
#Print results into latex file 

results <- xtable::xtable(arDF, 
                          caption ="Approximation ratios against known optimal tour cost.",
                          label = "ar",
                          digits=c(0,0,rep(3,ncol(arDF)-1)),
                          NA.string="NA"
)

xtable::print.xtable(results, 
                     type = "latex", 
                     file = paste(myResultsDir,"ar-results.txt", sep = "/"), 
                     caption.placement = "top",
                     NA.string="NA",
                     include.rownames=FALSE
)




########################################################################################################
#TourCost
########################################################################################################

tourCostDF <- as.data.frame(dataSets)
colnames(tourCostDF) <- c("DataSet")

for (a in Algos) {
  #DEBUG
  #a <- Algos[1]
  
  tourCostCol <- allDFAVG[which(allDFAVG$Algo == a),]$TourCost
  tourCost <- as.numeric(tourCostCol)
  tourCostColDF <- as.data.frame(as.numeric(tourCostCol))
  #colnames(timeColDF) <- c(a)
  colnames(tourCostColDF) <- myAlgoNames[which(Algos == a)]
  
  tourCostDF <- cbind(tourCostDF, tourCostColDF)
}

#add OPT col
optCol <- allDFOPT$TourCost

optCol[which(optCol == -1)] <- NA
optCol <- as.numeric(optCol)
optColDF <- as.data.frame(optCol)
colnames(optColDF) <- c("OPT")
tourCostDF <- cbind(tourCostDF, optColDF)

#Print results into latex file 

results <- xtable::xtable(tourCostDF, 
                          caption ="Approximation tour costs and the optimum tour cost if it is known.",
                          label = "tourCost",
                          digits=c(0,0,rep(1,ncol(tourCostDF)-1)),
                          NA.string="NA"
)

xtable::print.xtable(results, 
                     type = "latex", 
                     file = paste(myResultsDir,"tourCost-results.txt", sep = "/"), 
                     caption.placement = "top",
                     NA.string="NA",
                     include.rownames=FALSE
)






########################################################################################################
#AWDfromVTX1
########################################################################################################

awd1DF <- as.data.frame(dataSets)
colnames(awd1DF) <- c("DataSet")

for (a in Algos) {
  #DEBUG
  #a <- Algos[1]
  
  awd1Col <- allDFAVG[which(allDFAVG$Algo == a),]$AWDfromVTX1
  awd1 <- as.numeric(awd1Col)
  awd1ColDF <- as.data.frame(as.numeric(awd1Col))
  #colnames(timeColDF) <- c(a)
  colnames(awd1ColDF) <- myAlgoNames[which(Algos == a)]
  
  awd1DF <- cbind(awd1DF, awd1ColDF)
}

#add OPT col
optCol <- allDFOPT$AWDfromVTX1


optCol <- as.numeric(optCol)
optColDF <- as.data.frame(optCol)
colnames(optColDF) <- c("OPT")
awd1DF <- cbind(awd1DF, optColDF)

#Print results into latex file 

results <- xtable::xtable(awd1DF, 
                          caption ="AWD (from vertex 1) costs for the approximations and for the optimum tour cost if it is known.",
                          label = "awd1",
                          digits=c(0,0,rep(1,ncol(awd1DF)-1)),
                          NA.string="NA"
)

xtable::print.xtable(results, 
                     type = "latex", 
                     file = paste(myResultsDir,"awd1-results.txt", sep = "/"), 
                     caption.placement = "top",
                     NA.string="NA",
                     include.rownames=FALSE
)


########################################################################################################
#MinAWD 
########################################################################################################

minAwdDF <- as.data.frame(dataSets)
colnames(minAwdDF) <- c("DataSet")

for (a in Algos) {
  #DEBUG
  #a <- Algos[1]
  
  minAwdCol <- allDFAVG[which(allDFAVG$Algo == a),]$MinAWD
  minAwd <- as.numeric(minAwdCol)
  minAwdColDF <- as.data.frame(as.numeric(minAwdCol))
  #colnames(timeColDF) <- c(a)
  colnames(minAwdColDF) <- myAlgoNames[which(Algos == a)]
  
  minAwdDF <- cbind(minAwdDF, minAwdColDF)
}

#add OPT col
optCol <- allDFOPT$MinAWD


optCol <- as.numeric(optCol)
optColDF <- as.data.frame(optCol)
colnames(optColDF) <- c("OPT")
minAwdDF <- cbind(minAwdDF, optColDF)

#Print results into latex file 

results <- xtable::xtable(minAwdDF, 
                          caption ="minAWD costs for the approximations and for the optimum tour cost if it is known.",
                          label = "minAwd",
                          digits=c(0,0,rep(1,ncol(minAwdDF)-1)),
                          NA.string="NA"
)

xtable::print.xtable(results, 
                     type = "latex", 
                     file = paste(myResultsDir,"minAwd-results.txt", sep = "/"), 
                     caption.placement = "top",
                     NA.string="NA",
                     include.rownames=FALSE
)





########################################################################################################
#DelEdgeP
########################################################################################################


delEdgePDF <- as.data.frame(dataSets)
colnames(delEdgePDF) <- c("DataSet")

for (a in Algos) {
  #DEBUG
  #a <- Algos[1]
  
  delEdgePCol <- allDFAVG[which(allDFAVG$Algo == a),]$DelEdgeP
  delEdgeP <- as.numeric(delEdgePCol)
  delEdgePColDF <- as.data.frame(as.numeric(delEdgePCol))
  #colnames(timeColDF) <- c(a)
  colnames(delEdgePColDF) <- myAlgoNames[which(Algos == a)]
  
  delEdgePDF <- cbind(delEdgePDF, delEdgePColDF)
}

#add OPT col
optCol <- allDFOPT$DelEdgeP


optCol <- as.numeric(optCol)
optColDF <- as.data.frame(optCol)
colnames(optColDF) <- c("OPT")
delEdgePDF <- cbind(delEdgePDF, optColDF)

#Print results into latex file 

results <- xtable::xtable(delEdgePDF, 
                          caption ="Percentage of tour edges that come from Delaunay Triangulation edges.",
                          label = "delEdgeP",
                          digits=c(0,0,rep(2,ncol(delEdgePDF)-1)),
                          NA.string="NA"
)

xtable::print.xtable(results, 
                     type = "latex", 
                     file = paste(myResultsDir,"delEdgeP-results.txt", sep = "/"), 
                     caption.placement = "top",
                     NA.string="NA",
                     include.rownames=FALSE
)


#print(xtable(tableResults), include.rownames=FALSE)












########################################################################################################



