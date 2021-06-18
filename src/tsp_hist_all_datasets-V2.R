########################################################################################################
#Clean the environment
remove(list = ls())
########################################################################################################
library(TSP)
library(stringr) 
library(modeest)
library(moments)
library(xtable)
########################################################################################################
myWorkDir <- "/my_hd1/my_dir/my_prg/Camerino/tsp/dataset-hist" #Set the working dir
myDataDir <- "/my_hd1/my_dir/my_prg/Camerino/tsp/dataset-hist/data" #Set the working dir
setwd(myWorkDir)

###################################### MY FUNCTIONS BEGIN #################################################
#########################################################################################################
# Read next line from a connection.
#
# @param con [\code{connection}]\cr
#   Connection object.
# @return [\code{character(1)}] Next line of text from \code{con} or \code{NULL} if no
#   more lines of text are available.
#########################################################################################################
next_line = function(con) {
  str_trim(readLines(con, n = 1, warn = FALSE))
}

# Peek at next line of text from a connection.
#
# @param con [\code{connection}]\cr
#   Connection object.
# @return [\code{character(1)}] Next line of text from \code{con} or \code{NULL} if no
#   more lines of text are available. The line is not removed
#   from the connection so the next call to \code{next_line} will
#   return it.
peek_line = function(con) {
  line = readLines(con, 1, warn = FALSE)
  pushBack(line, con)
  line
}

# Check if one more line exists.
#
# @param con [\code{connection}]\cr
#   Connection object.
# @return [\code{logical(1)}] Boolean value indicating
#   whether another line in connection exists.
next_line_exists = function(con) {
  length(peek_line(con)) > 0
}


#########################################################################################################
# Read in a TSPLIB style Traveling Salesman Problem tour from a file
#
# @param path [\code{character(1)}]\cr
#   Filename of file containing a TSP tour.
# @return [\code{\link[TSP]{TOUR}}]
#  TOUR object from package TSP, containing order of cities, tour length and
#  method name that generated this solution.
#########################################################################################################
my_read_tsplib_tour = function(path) {
  
  # path <- "xqf131.tour"
  # determine file extension
  splitted = strsplit(path, ".", fixed = TRUE)[[1]]
  file_ext = splitted[length(splitted)]
  con = file(path, open = "r")
  on.exit(close(con))
  
  tsp_tour = list()
  # tour is given in tsplib format
  if(file_ext == "tour") {
    dimension = NULL
    
    line = next_line(con)
    while (str_trim(peek_line(con)) != "TOUR_SECTION") {
      field = str_split(line, ":")[[1]]
      field_name = str_trim(field[1])
      field_value = str_trim(field[2])
      # save dimension (neccessary to extract tour )
      if(field_name == "DIMENSION") {
        dimension = as.integer(field_value)
      }
      line = next_line(con)
    }
    
    while (next_line_exists(con)) {
      line = next_line(con)
      
      # what section is it?
      if (line == "TOUR_SECTION") {
        tsp_tour = read_tsplib_tour_section(con)
      }
    }
    
    # tour is given in concorde format
  } else if(file_ext == "sol") {
    line = next_line(con)
    dimension = as.integer(line)
    tsp_tour = numeric()
    while (next_line_exists(con)) {
      line = next_line(con)
      sub_tour = as.integer(unlist(strsplit(line, " ")[[1]]))
      tsp_tour = c(tsp_tour, sub_tour)
    }
    # concorde returns nodes 0 to (n-1), but we need nodes between 1 and n
    tsp_tour = tsp_tour + 1
  } else {
    stop("BAM! Unknown tour file format.")
  }
  return(as.integer(tsp_tour))
}

# read_tsplib_tour_section - extracts the nodes (cities) on the TSP tour
#
# @param con - open connection to file
# @param dimension
#
# @return vector of nodes on the TSP tour.
read_tsplib_tour_section = function(con) {
  nodes = numeric()
  line = next_line(con)
  while(line != "-1") {
    # we have to do some creepy stuff here to achieve correct tour
    temp = as.integer(c(str_split(str_trim(str_split(line, "  ")[[1]]), " "),
                        recursive = TRUE))
    nodes = c(nodes, temp)
    line = next_line(con)
  }
  return(nodes)
}


#https://stackoverflow.com/questions/2547402/how-to-find-the-statistical-mode/8189441#8189441
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#########################################################################################################

###################################### MY FUNCTIONS END #################################################


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

#Select all or certain dataset group
dsetname <- c(tspDataSets, bigTSPDataSets, bigCustomTSPDataSets )

#dsetname <- bigCustomTSPDataSets
lenDF <- length(dsetname)
tspDF <- data.frame(cbind(dsetname, rep(0, lenDF), rep(0, lenDF), rep(0, lenDF), 
                                    rep(0, lenDF), rep(0, lenDF), rep(0, lenDF)))
colnames(tspDF) <- c("dset","N","mean","med","std","skw","krt")

tspNo <- 1
for (d in dsetname) {
#for (d in tspDataSets) {
  
  #DEBUG  
  #d <- "usa13509" 
  #d <-"pr1002"
  #d <- "myRND-400"  
  #Set the tsp file name here without .tsp!!!
  
  #DEBUG theTSP <- dsetname[tspNo]
  
  theTSP <- d  #OK
  
  theTSPFileName <- paste0(theTSP,".tsp")
  
  cat(tspNo,"--- Opening the TSP:", theTSP, "\n")
  theTSPObject <- read_TSPLIB(paste(myDataDir,theTSPFileName, sep = "/"))
  theTSPcoordDF <- as.data.frame(theTSPObject)
  #insert vertex ids
  theTSPcoordDF[,3] <- 1:nrow(theTSPcoordDF)
  colnames(theTSPcoordDF) <- c("x", "y", "vno")
  
  theTSPcoordMTX <- as.matrix(theTSPcoordDF)
  
  #For plot axis
  xmin <- min(theTSPcoordMTX[,1])
  xmax <- max(theTSPcoordMTX[,1])
  ymin <- min(theTSPcoordMTX[,2])
  ymax <- max(theTSPcoordMTX[,2])
  
  theDistMTX <- as.matrix(dist(theTSPcoordMTX[,1:2], diag = TRUE, upper = TRUE))
  x <- as.vector(theDistMTX)
  
  #Mode estimation may take long time
  #xMode <- mlv(x[x != 0], method = "mfv")
  xMean <- mean(x[x != 0])
  xMedian <- median(x[x != 0])
  xSd <- sd(x[x != 0])
  xSkw <- skewness(x[x != 0])
  xKrt <- kurtosis(x[x!=0])
  
  #############################################################################
  #Plot the hist and save it if you want
  #############################################################################
  #png(paste0(theTSP,"-dist-hist.png"), width = myWidth, height = myHeight)
  
  #For EPS output
  cairo_ps(filename = (paste0(theTSP,"-dist-hist.eps")),
           width = 16,
           height = 10,
           pointsize = 12,
           fallback_resolution = 600)
  
  
  h <- hist(x[x != 0],
            main=paste0(theTSP, " - mean=",round(xMean,2) ," - med=",round(xMedian,2),
                        " - std=",round(xSd,2), " - skw=",round(xSkw,2), " - krt=",round(xKrt,2)),
            xlab="Edge Distance",
  )
  
  # h <- hist(x,
  #           main=paste0(theTSP, " - mean=",round(xMean,2) ," - med=",round(xMedian,2), 
  #                       " - mode=",round(xMode,2)," - std=",round(xSd,2)),
  #           xlab="Edge Distance",
  # )
  
  abline(v = xMean, col = "blue", lwd = 2)
  #text(mean(x), max(h$counts)+5, "Mean", col = "blue")
  
  abline(v = xMedian , col = "red", lwd = 2)
  #text(median(x), max(h$counts)-5, "Med", col = "red")
  
  legend("topright", legend=c("Mean", "Median"), col=c("blue", "red"), lwd=2)
  
  # abline(v = xMode, col = "green", lwd = 2)
  # legend("topright", legend=c("Mean", "Median", "Mode"), col=c("blue", "red", "green"), lwd=2)
  
  dev.off()
  
  # plot(density(x[x!=0]), paste0("Density plot ", theTSP))
  #############################################################################
  
  
  
  
  #############################################################################
  #Plot the vtxs and save it if you want
  #############################################################################
  #png(paste0(theTSP,"-vtx.png"), width = myWidth, height = myHeight)
  
  #For EPS output
  
  #for lon lat data rotate 90 counter clockwise
  # (x,y) -> (-y, x)
  
  
  if (d == "usa13509") {
    
    theTSPcoordMTXusa <- theTSPcoordMTX
    
    theTSPcoordMTXusa[,1] <- ymax - theTSPcoordMTX[,2]
    theTSPcoordMTXusa[,2] <- theTSPcoordMTX[,1]
    
    xmin2 <- min(theTSPcoordMTXusa[,1])
    xmax2 <- max(theTSPcoordMTXusa[,1])
    
    ymin2 <- min(theTSPcoordMTXusa[,2])
    ymax2 <- max(theTSPcoordMTXusa[,2])
    
    
    cairo_ps(filename = (paste0(theTSP,"-vtx-plot.eps")),
             width = 16,
             height = 10,
             pointsize = 12,
             fallback_resolution = 600)
    
    plot(theTSPcoordMTXusa[,1:2], col = "black", pch=19, cex=0.05,
         xlab = "x-coordinate", 
         ylab = "y-coordinate", 
         main=paste0(theTSP, " - mean=",round(xMean,2) ," - med=",round(xMedian,2), 
                     " - std=",round(xSd,2), " - skw=",round(xSkw,2), " - krt=",round(xKrt,2)),
         xlim=c(xmin2, xmax2), 
         ylim=c(ymin2, ymax2), asp=1
    )
    
    dev.off()
    
  }else{
    
    cairo_ps(filename = (paste0(theTSP,"-vtx-plot.eps")),
             width = 16,
             height = 10,
             pointsize = 12,
             fallback_resolution = 600)
    
    plot(theTSPcoordMTX[,1:2], col = "black", pch=19, cex=0.05,
         xlab = "x-coordinate", 
         ylab = "y-coordinate", 
         main=paste0(theTSP, " - mean=",round(xMean,2) ," - med=",round(xMedian,2), 
                     " - std=",round(xSd,2), " - skw=",round(xSkw,2), " - krt=",round(xKrt,2)),
         xlim=c(xmin, xmax), 
         ylim=c(ymin, ymax), asp=1
    )
    
    dev.off()
    
  }
  
  
  tspDF[tspNo,2:7] <- c(nrow(theTSPcoordMTX), xMean, xMedian, xSd, xSkw, xKrt)
  
  #############################################################################
  
  tspNo <- tspNo + 1
}


#############################################################################
#Print results into latex file 
for (i in 3:7) {
  tspDF[, i] <- as.numeric(tspDF[, i])
}

results <- xtable::xtable(tspDF, 
                          caption ="Statistics of the edge distances for datasets.",
                          label = "stat-dset",
                          digits=c(0,0,0, rep(2,5)), #precision for each column of the df
                          NA.string="NA"
)

xtable::print.xtable(results, 
                     type = "latex", 
                     file = paste(myDataDir,"edge-dist-stat-big-custom.txt", sep = "/"), 
                     caption.placement = "top",
                     NA.string="NA",
                     include.rownames=FALSE
)
#############################################################################


# 
# 
# #############################################################################
# #Distribution Test for usa13509
# #############################################################################
# 
# #https://stats.stackexchange.com/questions/132652/how-to-determine-which-distribution-fits-my-data-best
# 
# library(fitdistrplus)
# library(logspline)
# 
# d <- "usa13509" 
#  
# #Set the tsp file name here without .tsp!!!
# theTSP <- "usa13509"  #OK
# 
# theTSPFileName <- paste0(theTSP,".tsp")
# 
# cat("Opening the TSP:", theTSP, "\n")
# theTSPObject <- read_TSPLIB(paste(myDataDir,theTSPFileName, sep = "/"))
# theTSPcoordDF <- as.data.frame(theTSPObject)
# #insert vertex ids
# theTSPcoordDF[,3] <- 1:nrow(theTSPcoordDF)
# colnames(theTSPcoordDF) <- c("x", "y", "vno")
# 
# theTSPcoordMTX <- as.matrix(theTSPcoordDF)
# 
# 
# #For plot axis
# xmin <- min(theTSPcoordMTX[,1])
# xmax <- max(theTSPcoordMTX[,1])
# 
# ymin <- min(theTSPcoordMTX[,2])
# ymax <- max(theTSPcoordMTX[,2])
# 
# 
# 
# theDistMTX <- as.matrix(dist(theTSPcoordMTX[,1:2], diag = TRUE, upper = TRUE))
# 
# x <- as.vector(theDistMTX)
# 
# 
# xNZ <- x[x != 0] #Remove zeros
# 
# ##################################################
# # Distr fit is so slow for usa13506
# ###############################################
# # 
# # descdist(xNZ, discrete = FALSE)
# # summary statistics
# # ------
# #   min:  2.777   max:  575461.2 
# # median:  132805.8 
# # mean:  159409.2 
# # estimated sd:  109330.3 
# # estimated skewness:  1.055289 
# # estimated kurtosis:  3.591165 
# 
# # 
# # Took too long for these:
# # fit.weibull <- fitdist(xNZ, "weibull")
# # fit.norm <- fitdist(xNZ, "norm")
# # # Now inspect the fit for the normal:
# # plot(fit.norm)
# # 
# # #And for the Weibull fit:
# # plot(fit.weibull)
# # # 
# # # Both look good but judged by the QQ-Plot, 
# # # the Weibull maybe looks a bit better, especially at the tails. 
# # # Correspondingly, the AIC of the Weibull fit is lower compared to the normal fit:
# # #   
# # fit.weibull$aic
# # fit.norm$aic
# 
# library(MASS)
# 
# res <- fitdistr(xNZ,"weibull",lower=0.001)
# 
# 
# #Too long
# 
# fit.norm <- fitdist(xNZ, "norm")
# fit.norm
# 
# plot(fit.norm)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# descdist(xNZ)
# 
# fn <- fitdistr(xNZ, "normal")
# plot(fn)
# 
# 
# fw <- fitdist(xNZ, "weibull")
# 
# fg <- fitdist(xNZ, "gamma")
# 
# fe <- fitdist(xNZ, "exp")
# 
# par(mfrow = c(2, 2))
# 
# plot.legend <- c("Weibull", "gamma", "expo")
# 
# denscomp(list(fw, fg, fe), legendtext = plot.legend)
# 
# qqcomp(list(fw, fg, fe), legendtext = plot.legend)
# 
# cdfcomp(list(fw, fg, fe), legendtext = plot.legend)
# 
# ppcomp(list(fw, fg, fe), legendtext = plot.legend)
# 
# 
# 
# 
# 









