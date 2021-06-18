#########################################################################################################
#Functions
#########################################################################################################


########################################################################################################
#Given the mtx that formed from reading TSP object and the mtx that contains the coords of the tour vtxs
#this function return the vector that contains vertex ids for the tour as a permutation not as a cycle
#TSP mtx is not cyclical
#Coord mtx for the tour is and not cyclical
########################################################################################################
coord2vtxPerm <- function(theTSPMTX, theCoordMTX) {
  #Lets find the vertex ids on the tour
  nVtx <-  nrow(theCoordMTX)
  vtxIds <- vector()
  for (i in 1:nVtx) {
    idx <- which(theTSPMTX[,1] ==  theCoordMTX[i,1] & theTSPMTX[,2] ==  theCoordMTX[i,2])
    
    if (length(is.na(idx)) == 0) {
      cat("ERROR: The point", i ,"(",theCoordMTX[i,1], theCoordMTX[i,2],") was not found in coords!!!\n")
    }else{
      vtxIds[i] <- idx 
    }
  }
  return(vtxIds)
}
########################################################################################################


########################################################################################################
#Given (lon, lat) or (x, y) in a data frame sort points
#https://raw.githubusercontent.com/skgrange/gissr/master/R/sort_points.R
#Example call:
#sortedBTMtx <- my_sort_points(BTDF, y="lat", x="lon", bsx=BSMtx[1,1], bsy=BSMtx[1,2])
########################################################################################################

my_sort_points <- function(df, x = "V1", y = "V2", clockwise = TRUE) {
  #we need to define delta sector
  #if any two points fall into same sector
  # we need to give the one that is closest to the "last point" visited
  # NA check, if NAs drop them
  if (any(is.na(c(df[, y], df[, x])))) {
    
    # Remove NAs
    df <- df[!(is.na(df[, y]) & is.na(df[, x])), ]
    
    # Raise warning
    warning("Missing coordinates were detected and have been removed.", 
            call. = FALSE)
    
    # Check 
    if (nrow(df) == 0) stop("There are no valid coordinates.", call. = FALSE)
    
  }
  
  #Get centre (-oid) point of points
  x_centre <- mean(df[, x])
  y_centre <- mean(df[, y])
  
  
  
  
  # Calculate deltas
  df$x_delta <- df[, x] - x_centre
  df$y_delta <- df[, y] - y_centre
  
  # Resolve angle, in radians
  df$angle <- atan2(df$y_delta, df$x_delta)
  df$angle_degrees <- df$angle * 180 / pi
  
  # Arrange by angle
  if (clockwise) {
    
    df <- df[order(df$angle, decreasing = TRUE), ]
    
  } else {
    
    df <- df[order(df$angle, decreasing = FALSE), ]
    
  }
  
  #Drop intermediate variables
  df[, c("x_delta", "y_delta", "angle", "angle_degrees")] <- NULL
  
  # Return
  df
  
}





#########################################################################################################
# This function sums the cost of the rings generated in the list
#########################################################################################################

my_sum_cost_rings <- function(theRingList) {
  
  nr <- length(theRingList)
  sumCost <- 0
  for (i in (1:nr)) {
    nElem <- nrow(theRingList[[i]]) - 1
    sumCost <- sumCost + -1 * my_get_cost_perm(theRingList[[i]][1:nElem,3])
  }
  
  return(sumCost)
}
#########################################################################################################





#########################################################################################################
#The fitness function
# evaluates the total distance given the permutation
# Do not repeat the first point at the end like a polygon 
# This function will take care of it!!!!
#########################################################################################################
my_get_cost_perm <- function(thePerm) {
  
  nNode <- length(thePerm)
  
  totDist <- 0
  
  startNode <- thePerm[1]
  
  for (k in 2:nNode) {
    
    stopNode <- thePerm[k]
    d <- sqrt( ((theTSPcoordDF[stopNode,1] -  theTSPcoordDF[startNode,1]) * (theTSPcoordDF[stopNode,1] -  theTSPcoordDF[startNode,1])) + 
                 ((theTSPcoordDF[stopNode,2] -  theTSPcoordDF[startNode,2]) * (theTSPcoordDF[stopNode,2] -  theTSPcoordDF[startNode,2])) )
    totDist <- totDist + d
    startNode <- thePerm[k]
  }
  
  #Close the tour
  startNode <- thePerm[nNode]
  stopNode  <- thePerm[1]
  d <- sqrt( ((theTSPcoordDF[stopNode,1] -  theTSPcoordDF[startNode,1]) * (theTSPcoordDF[stopNode,1] -  theTSPcoordDF[startNode,1])) + 
               ((theTSPcoordDF[stopNode,2] -  theTSPcoordDF[startNode,2]) * (theTSPcoordDF[stopNode,2] -  theTSPcoordDF[startNode,2])) )
  totDist <- totDist + d
  
  #Return the cost of the tour
  return(-totDist)
}

#########################################################################################################

#########################################################################################################
#The awt (avg waiting time) function
# evaluates the avg waiting time for the vertices on the perm
# Of course since we do not assume any speed the awt will be in terms of distance
# Do not repeat the first point at the end like a polygon 
# This function will take care of it!!!!
#########################################################################################################
my_get_awt_perm <- function(thePerm) {
  
  nNode <- length(thePerm)
  
  sumDist <- 0
  pathSoFar <- 0
  
  startNode <- thePerm[1]
  
  for (k in 2:nNode) {
    
    stopNode <- thePerm[k]
    d <- sqrt( ((theTSPcoordDF[stopNode,1] -  theTSPcoordDF[startNode,1]) * (theTSPcoordDF[stopNode,1] -  theTSPcoordDF[startNode,1])) + 
                 ((theTSPcoordDF[stopNode,2] -  theTSPcoordDF[startNode,2]) * (theTSPcoordDF[stopNode,2] -  theTSPcoordDF[startNode,2])) )
    pathSoFar <- pathSoFar + d
    sumDist <- sumDist + pathSoFar
    startNode <- thePerm[k]
  }
  
  
  #Close the tour
  startNode <- thePerm[nNode]
  stopNode  <- thePerm[1]
  d <- sqrt( ((theTSPcoordDF[stopNode,1] -  theTSPcoordDF[startNode,1]) * (theTSPcoordDF[stopNode,1] -  theTSPcoordDF[startNode,1])) + 
               ((theTSPcoordDF[stopNode,2] -  theTSPcoordDF[startNode,2]) * (theTSPcoordDF[stopNode,2] -  theTSPcoordDF[startNode,2])) )
  
  pathSoFar <-pathSoFar + d
  sumDist <- sumDist + pathSoFar
  
  #Return the awt of the tour
  awt <- sumDist / nNode
  
  return(awt)
}

#########################################################################################################



#########################################################################################################
#The awt (acyclical avg waiting time) function
# evaluates the acyclical avg waiting time for the vertices on the perm
# Of course since we do not assume any speed the awt will be in terms of distance
# Do not repeat the first point at the end like a polygon 
# This function will take care of it!!!!
# We care for the time that all the vertices are served except the starting vertex
# return path from the last vertex is not important
#########################################################################################################
my_get_awt_perm_acyclical <- function(thePerm) {
  
  nNode <- length(thePerm)
  
  sumDist <- 0
  pathSoFar <- 0
  
  startNode <- thePerm[1]
  
  for (k in 2:nNode) {
    
    stopNode <- thePerm[k]
    d <- sqrt( ((theTSPcoordDF[stopNode,1] -  theTSPcoordDF[startNode,1]) * (theTSPcoordDF[stopNode,1] -  theTSPcoordDF[startNode,1])) + 
                 ((theTSPcoordDF[stopNode,2] -  theTSPcoordDF[startNode,2]) * (theTSPcoordDF[stopNode,2] -  theTSPcoordDF[startNode,2])) )
    pathSoFar <- pathSoFar + d
    sumDist <- sumDist + pathSoFar
    startNode <- thePerm[k]
  }
  
  #UPDATE: There is no need to "close the tour".
  #Once you serve the last vertex the job is done!!!!
  # #Close the tour
  # startNode <- thePerm[nNode]
  # stopNode  <- thePerm[1]
  # d <- sqrt( ((theTSPcoordDF[stopNode,1] -  theTSPcoordDF[startNode,1]) * (theTSPcoordDF[stopNode,1] -  theTSPcoordDF[startNode,1])) + 
  #              ((theTSPcoordDF[stopNode,2] -  theTSPcoordDF[startNode,2]) * (theTSPcoordDF[stopNode,2] -  theTSPcoordDF[startNode,2])) )
  # 
  # pathSoFar <-pathSoFar + d
  # sumDist <- sumDist + pathSoFar
  
  #Return the awt of the tour
  awt <- sumDist / (nNode - 1)
  
  return(awt)
}

#########################################################################################################

#########################################################################################################
#The awt (avg waiting time) function
# evaluates the avg waiting time (distance!!!) for the vertices on the perm starting from the vertex no 1
# TSPLIB data set optimum tours are from vtx 1
# Of course since we do not assume any speed the awt will be in terms of distance
# Do not repeat the first point at the end like a polygon 
# This function will take care of it!!!!
#########################################################################################################
my_get_awt_from_vtx1 <- function(thePerm) {
  
  nNode <- length(thePerm)
  
  #Lets rotate the tour so that vtxno 1 will be the first vertex of the tour
  vtx1pos <- which(thePerm == 1)
  if (vtx1pos != 1) {
    thePerm <- c(thePerm[vtx1pos:nNode], thePerm[1:(vtx1pos-1)])
  }
  
  #Now we calculate the AWD as usual
  sumDist <- 0
  pathSoFar <- 0
  
  startNode <- thePerm[1]
  
  for (k in 2:nNode) {
    
    stopNode <- thePerm[k]
    d <- sqrt( ((theTSPcoordDF[stopNode,1] -  theTSPcoordDF[startNode,1]) * (theTSPcoordDF[stopNode,1] -  theTSPcoordDF[startNode,1])) + 
                 ((theTSPcoordDF[stopNode,2] -  theTSPcoordDF[startNode,2]) * (theTSPcoordDF[stopNode,2] -  theTSPcoordDF[startNode,2])) )
    pathSoFar <- pathSoFar + d
    sumDist <- sumDist + pathSoFar
    startNode <- thePerm[k]
  }
  
  #Close the tour
  startNode <- thePerm[nNode]
  stopNode  <- thePerm[1]
  d <- sqrt( ((theTSPcoordDF[stopNode,1] -  theTSPcoordDF[startNode,1]) * (theTSPcoordDF[stopNode,1] -  theTSPcoordDF[startNode,1])) + 
               ((theTSPcoordDF[stopNode,2] -  theTSPcoordDF[startNode,2]) * (theTSPcoordDF[stopNode,2] -  theTSPcoordDF[startNode,2])) )
  
  pathSoFar <-pathSoFar + d
  sumDist <- sumDist + pathSoFar
  
  #Return the awt of the tour
  awt <- sumDist / nNode
  
  return(awt)
}
#########################################################################################################





#########################################################################################################
#The awt (avg waiting time) function
# evaluates the avg waiting time (distance!!!) for the vertices on the perm starting from the vertex no 1
# TSPLIB data set optimum tours are from vtx 1
# Of course since we do not assume any speed the awt will be in terms of distance
# Do not repeat the first point at the end like a polygon 
# This function will take care of it!!!!
#########################################################################################################
my_get_awt_from_vtx1_acyclical <- function(thePerm) {
  
  nNode <- length(thePerm)
  
  #Lets rotate the tour so that vtxno 1 will be the first vertex of the tour
  vtx1pos <- which(thePerm == 1)
  if (vtx1pos != 1) {
    thePerm <- c(thePerm[vtx1pos:nNode], thePerm[1:(vtx1pos-1)])
  }
  
  #Now we calculate the AWD as usual
  sumDist <- 0
  pathSoFar <- 0
  
  startNode <- thePerm[1]
  
  for (k in 2:nNode) {
    
    stopNode <- thePerm[k]
    d <- sqrt( ((theTSPcoordDF[stopNode,1] -  theTSPcoordDF[startNode,1]) * (theTSPcoordDF[stopNode,1] -  theTSPcoordDF[startNode,1])) + 
                 ((theTSPcoordDF[stopNode,2] -  theTSPcoordDF[startNode,2]) * (theTSPcoordDF[stopNode,2] -  theTSPcoordDF[startNode,2])) )
    pathSoFar <- pathSoFar + d
    sumDist <- sumDist + pathSoFar
    startNode <- thePerm[k]
  }
  
  # #Close the tour
  # startNode <- thePerm[nNode]
  # stopNode  <- thePerm[1]
  # d <- sqrt( ((theTSPcoordDF[stopNode,1] -  theTSPcoordDF[startNode,1]) * (theTSPcoordDF[stopNode,1] -  theTSPcoordDF[startNode,1])) + 
  #              ((theTSPcoordDF[stopNode,2] -  theTSPcoordDF[startNode,2]) * (theTSPcoordDF[stopNode,2] -  theTSPcoordDF[startNode,2])) )
  # 
  # pathSoFar <-pathSoFar + d
  # sumDist <- sumDist + pathSoFar
  
  #Return the awt of the tour
  
  awt <- sumDist / (nNode - 1)
  
  return(awt)
}
#########################################################################################################









#########################################################################################################
#The awt (avg waiting time) function
# evaluates the min avg waiting time for the vertices on the perm
# Of course since we do not assume any speed the awt will be in terms of distance
# Do not repeat the first point at the end like a polygon 
# This function will take care of it!!!!
#########################################################################################################
my_get_min_awt_perm <- function(thePerm) {
  
  #DEBUG 
  #  thePerm <- c(1,2,3,4,5,6,7)
  
  #Update: consider clockwise and anti-clockwise directions!!!!
  
  
  
  
  #########################################################################################################
  #Clockwise
  #########################################################################################################
  thePermCW <- thePerm
  nNode <- length(thePermCW)
  minAWTCW <- Inf
  minAWTpermCW <- thePermCW
  
  for (startIdx in 1:nNode) {
    
    #create the tour array and calculate the AWT
    permHead <- thePermCW[startIdx:nNode]
    if ((startIdx - 1) > 0 ) permTail <- thePermCW[1:(startIdx-1)] 
    else permTail <- NULL
    
    currentPerm <- c(permHead, permTail)
    #DEBUG     cat(startIdx,"---", currentPerm,"\n")
    
    
    sumDist <- 0
    pathSoFar <- 0
    
    startNode <- currentPerm[1]
    
    for (k in 2:nNode) {
      
      stopNode <- currentPerm[k]
      d <- sqrt( ((theTSPcoordDF[stopNode,1] -  theTSPcoordDF[startNode,1]) * (theTSPcoordDF[stopNode,1] -  theTSPcoordDF[startNode,1])) + 
                   ((theTSPcoordDF[stopNode,2] -  theTSPcoordDF[startNode,2]) * (theTSPcoordDF[stopNode,2] -  theTSPcoordDF[startNode,2])) )
      pathSoFar <- pathSoFar + d
      sumDist <- sumDist + pathSoFar
      startNode <- currentPerm[k]
    }
    
    #Close the tour
    startNode <- currentPerm[nNode]
    stopNode  <- currentPerm[1]
    d <- sqrt( ((theTSPcoordDF[stopNode,1] -  theTSPcoordDF[startNode,1]) * (theTSPcoordDF[stopNode,1] -  theTSPcoordDF[startNode,1])) + 
                 ((theTSPcoordDF[stopNode,2] -  theTSPcoordDF[startNode,2]) * (theTSPcoordDF[stopNode,2] -  theTSPcoordDF[startNode,2])) )
    
    pathSoFar <- pathSoFar + d
    sumDist <- sumDist + pathSoFar
    
    #Return the awt of the tour
    awt <- sumDist / nNode
    
    
    if (awt < minAWTCW) {
      minAWTCW <- awt
      minAWTCWStartVTX <- currentPerm[1]
      minAWTpermCW <- currentPerm
    }
    
  }
  
  resultCW <- list(3)
  resultCW[[1]] <- minAWTCW
  resultCW[[2]] <- minAWTCWStartVTX
  resultCW[[3]] <- minAWTpermCW
  
  
  #########################################################################################################
  #Anti-Clockwise
  #########################################################################################################
  nNode <- length(thePerm)
  thePermACW <- c(thePerm[1] , rev(thePerm[2:nNode]))
  
  minAWTACW <- Inf
  minAWTpermACW <- thePermACW
  
  for (startIdx in 1:nNode) {
    
    #create the tour array and calculate the AWT
    permHead <- thePermACW[startIdx:nNode]
    if ((startIdx - 1) > 0 ) permTail <- thePermACW[1:(startIdx-1)] 
    else permTail <- NULL
    
    currentPerm <- c(permHead, permTail)
    #DEBUG     cat(startIdx,"---", currentPerm,"\n")
    
    
    sumDist <- 0
    pathSoFar <- 0
    
    startNode <- currentPerm[1]
    
    for (k in 2:nNode) {
      
      stopNode <- currentPerm[k]
      d <- sqrt( ((theTSPcoordDF[stopNode,1] -  theTSPcoordDF[startNode,1]) * (theTSPcoordDF[stopNode,1] -  theTSPcoordDF[startNode,1])) + 
                   ((theTSPcoordDF[stopNode,2] -  theTSPcoordDF[startNode,2]) * (theTSPcoordDF[stopNode,2] -  theTSPcoordDF[startNode,2])) )
      pathSoFar <- pathSoFar + d
      sumDist <- sumDist + pathSoFar
      startNode <- currentPerm[k]
    }
    
    #Close the tour
    startNode <- currentPerm[nNode]
    stopNode  <- currentPerm[1]
    d <- sqrt( ((theTSPcoordDF[stopNode,1] -  theTSPcoordDF[startNode,1]) * (theTSPcoordDF[stopNode,1] -  theTSPcoordDF[startNode,1])) + 
                 ((theTSPcoordDF[stopNode,2] -  theTSPcoordDF[startNode,2]) * (theTSPcoordDF[stopNode,2] -  theTSPcoordDF[startNode,2])) )
    
    pathSoFar <- pathSoFar + d
    sumDist <- sumDist + pathSoFar
    
    #Return the awt of the tour
    awt <- sumDist / nNode
    
    
    if (awt < minAWTACW) {
      minAWTACW <- awt
      minAWTACWStartVTX <- currentPerm[1]
      minAWTpermACW <- currentPerm
    }
    
  }
  
  resultACW <- list(3)
  resultACW[[1]] <- minAWTACW
  resultACW[[2]] <- minAWTACWStartVTX
  resultACW[[3]] <- minAWTpermACW
  
  #Pick the smallest of CW=Forward and ACW=Backward results
  result <- list(4)
  if (resultACW[[1]] < resultCW[[1]]){
    result <- resultACW
    result[[4]] <- 0 #Clock wise = Forward/Backward flag
  }else{
    
    result <- resultCW
    result[[4]] <- 1 #Clock wise flag
    if (resultACW[[1]] == resultCW[[1]]) { result[[4]] <- 2}
  }
  
  
  
  
  return(result)
}

#########################################################################################################




#########################################################################################################
#The awt (avg waiting time) function
# evaluates the min avg waiting time for the vertices on the perm
# Of course since we do not assume any speed the awt will be in terms of distance
# Do not repeat the first point at the end like a polygon 
# This function will take care of it!!!!
#########################################################################################################
my_get_min_awt_perm_fast <- function(thePerm) {
  
  #USES theDistMTX
  #DEBUG 
  #  thePerm <- c(1,2,3,4,5,6,7)
  
  #Update: consider clockwise and anti-clockwise directions!!!!
  
  #DEBUG thePerm <- concaveRingMergeNearestPtV2Tour
  
  # myLattice-20x20-400 0.451 sec elapsed
  #Normal awt finding takes 20.689 sec elapsed
  
  #########################################################################################################
  #Clockwise
  #########################################################################################################
  thePermCW <- thePerm
  nNode <- length(thePermCW)
  minAWTCW <- Inf
  minAWTpermCW <- thePermCW
  
  for (startIdx in 1:nNode) {
    
    #create the tour array and calculate the AWT
    permHead <- thePermCW[startIdx:nNode]
    if ((startIdx - 1) > 0 ) permTail <- thePermCW[1:(startIdx-1)] 
    else permTail <- NULL
    
    currentPerm <- c(permHead, permTail)
    #DEBUG     cat(startIdx,"---", currentPerm,"\n")
    
    
    sumDist <- 0
    pathSoFar <- 0
    
    startNode <- currentPerm[1]
    
    for (k in 2:nNode) {
      stopNode <- currentPerm[k]
      pathSoFar <- pathSoFar + theDistMTX[stopNode, startNode]
      sumDist <- sumDist + pathSoFar
      startNode <- currentPerm[k]
    }
    
    #Close the tour
    startNode <- currentPerm[nNode]
    stopNode  <- currentPerm[1]
    pathSoFar <- pathSoFar + theDistMTX[stopNode, startNode]
    sumDist <- sumDist + pathSoFar
    
    #Return the awt of the tour
    awt <- sumDist / nNode
    
    
    if (awt < minAWTCW) {
      minAWTCW <- awt
      minAWTCWStartVTX <- currentPerm[1]
      minAWTpermCW <- currentPerm
    }
    
  }
  
  resultCW <- list(3)
  resultCW[[1]] <- minAWTCW
  resultCW[[2]] <- minAWTCWStartVTX
  resultCW[[3]] <- minAWTpermCW
  
  
  #########################################################################################################
  #Anti-Clockwise
  #########################################################################################################
  nNode <- length(thePerm)
  thePermACW <- c(thePerm[1] , rev(thePerm[2:nNode]))
  
  minAWTACW <- Inf
  minAWTpermACW <- thePermACW
  
  for (startIdx in 1:nNode) {
    
    #create the tour array and calculate the AWT
    permHead <- thePermACW[startIdx:nNode]
    if ((startIdx - 1) > 0 ) permTail <- thePermACW[1:(startIdx-1)] 
    else permTail <- NULL
    
    currentPerm <- c(permHead, permTail)
    #DEBUG     cat(startIdx,"---", currentPerm,"\n")
    
    
    sumDist <- 0
    pathSoFar <- 0
    
    startNode <- currentPerm[1]
    
    for (k in 2:nNode) {
      
      stopNode <- currentPerm[k]
      pathSoFar <- pathSoFar + theDistMTX[stopNode, startNode]
      sumDist <- sumDist + pathSoFar
      startNode <- currentPerm[k]
    }
    
    #Close the tour
    startNode <- currentPerm[nNode]
    stopNode  <- currentPerm[1]
    pathSoFar <- pathSoFar + theDistMTX[stopNode, startNode]
    sumDist <- sumDist + pathSoFar
    
    #Return the awt of the tour
    awt <- sumDist / nNode
    
    
    if (awt < minAWTACW) {
      minAWTACW <- awt
      minAWTACWStartVTX <- currentPerm[1]
      minAWTpermACW <- currentPerm
    }
    
  }
  
  resultACW <- list(3)
  resultACW[[1]] <- minAWTACW
  resultACW[[2]] <- minAWTACWStartVTX
  resultACW[[3]] <- minAWTpermACW
  
  #Pick the smallest of CW=Forward and ACW=Backward results
  result <- list(4)
  if (resultACW[[1]] < resultCW[[1]]){
    result <- resultACW
    result[[4]] <- 0 #Clock wise = Forward/Backward flag
  }else{
    
    result <- resultCW
    result[[4]] <- 1 #Clock wise flag
    if (resultACW[[1]] == resultCW[[1]]) { result[[4]] <- 2}
  }
  
  
  
  
  return(result)
}

#########################################################################################################















#########################################################################################################
#The awt (avg waiting time) function
# evaluates the min avg waiting time for the vertices on the perm
# Of course since we do not assume any speed the awt will be in terms of distance
# Do not repeat the first point at the end like a polygon 
# This function will take care of it!!!!
#########################################################################################################
my_get_min_awt_perm_acyclical <- function(thePerm) {
  
  #DEBUG 
  #  thePerm <- c(1,2,3,4,5,6,7)
  
  #Update: consider clockwise and anti-clockwise directions!!!!
  
  
  
  
  #########################################################################################################
  #Clockwise
  #########################################################################################################
  thePermCW <- thePerm
  nNode <- length(thePermCW)
  minAWTCW <- Inf
  minAWTpermCW <- thePermCW
  
  for (startIdx in 1:nNode) {
    
    #create the tour array and calculate the AWT
    permHead <- thePermCW[startIdx:nNode]
    if ((startIdx - 1) > 0 ) permTail <- thePermCW[1:(startIdx-1)] 
    else permTail <- NULL
    
    currentPerm <- c(permHead, permTail)
    #DEBUG     cat(startIdx,"---", currentPerm,"\n")
    
    
    sumDist <- 0
    pathSoFar <- 0
    
    startNode <- currentPerm[1]
    
    for (k in 2:nNode) {
      
      stopNode <- currentPerm[k]
      d <- sqrt( ((theTSPcoordDF[stopNode,1] -  theTSPcoordDF[startNode,1]) * (theTSPcoordDF[stopNode,1] -  theTSPcoordDF[startNode,1])) + 
                   ((theTSPcoordDF[stopNode,2] -  theTSPcoordDF[startNode,2]) * (theTSPcoordDF[stopNode,2] -  theTSPcoordDF[startNode,2])) )
      pathSoFar <- pathSoFar + d
      sumDist <- sumDist + pathSoFar
      startNode <- currentPerm[k]
    }
    
    # #Close the tour
    # startNode <- currentPerm[nNode]
    # stopNode  <- currentPerm[1]
    # d <- sqrt( ((theTSPcoordDF[stopNode,1] -  theTSPcoordDF[startNode,1]) * (theTSPcoordDF[stopNode,1] -  theTSPcoordDF[startNode,1])) + 
    #              ((theTSPcoordDF[stopNode,2] -  theTSPcoordDF[startNode,2]) * (theTSPcoordDF[stopNode,2] -  theTSPcoordDF[startNode,2])) )
    # 
    # pathSoFar <- pathSoFar + d
    # sumDist <- sumDist + pathSoFar
    
    #Return the awt of the tour
    awt <- sumDist / (nNode - 1)
    
    
    if (awt < minAWTCW) {
      minAWTCW <- awt
      minAWTCWStartVTX <- currentPerm[1]
      minAWTpermCW <- currentPerm
    }
    
  }
  
  resultCW <- list(3)
  resultCW[[1]] <- minAWTCW
  resultCW[[2]] <- minAWTCWStartVTX
  resultCW[[3]] <- minAWTpermCW
  
  
  #########################################################################################################
  #Anti-Clockwise
  #########################################################################################################
  nNode <- length(thePerm)
  thePermACW <- c(thePerm[1] , rev(thePerm[2:nNode]))
  
  minAWTACW <- Inf
  minAWTpermACW <- thePermACW
  
  for (startIdx in 1:nNode) {
    
    #create the tour array and calculate the AWT
    permHead <- thePermACW[startIdx:nNode]
    if ((startIdx - 1) > 0 ) permTail <- thePermACW[1:(startIdx-1)] 
    else permTail <- NULL
    
    currentPerm <- c(permHead, permTail)
    #DEBUG     cat(startIdx,"---", currentPerm,"\n")
    
    
    sumDist <- 0
    pathSoFar <- 0
    
    startNode <- currentPerm[1]
    
    for (k in 2:nNode) {
      
      stopNode <- currentPerm[k]
      d <- sqrt( ((theTSPcoordDF[stopNode,1] -  theTSPcoordDF[startNode,1]) * (theTSPcoordDF[stopNode,1] -  theTSPcoordDF[startNode,1])) + 
                   ((theTSPcoordDF[stopNode,2] -  theTSPcoordDF[startNode,2]) * (theTSPcoordDF[stopNode,2] -  theTSPcoordDF[startNode,2])) )
      pathSoFar <- pathSoFar + d
      sumDist <- sumDist + pathSoFar
      startNode <- currentPerm[k]
    }
    
    # #Close the tour
    # startNode <- currentPerm[nNode]
    # stopNode  <- currentPerm[1]
    # d <- sqrt( ((theTSPcoordDF[stopNode,1] -  theTSPcoordDF[startNode,1]) * (theTSPcoordDF[stopNode,1] -  theTSPcoordDF[startNode,1])) + 
    #              ((theTSPcoordDF[stopNode,2] -  theTSPcoordDF[startNode,2]) * (theTSPcoordDF[stopNode,2] -  theTSPcoordDF[startNode,2])) )
    # 
    # pathSoFar <- pathSoFar + d
    # sumDist <- sumDist + pathSoFar
    
    #Return the awt of the tour
    awt <- sumDist / (nNode - 1)
    
    
    if (awt < minAWTACW) {
      minAWTACW <- awt
      minAWTACWStartVTX <- currentPerm[1]
      minAWTpermACW <- currentPerm
    }
    
  }
  
  resultACW <- list(3)
  resultACW[[1]] <- minAWTACW
  resultACW[[2]] <- minAWTACWStartVTX
  resultACW[[3]] <- minAWTpermACW
  
  #Pick the smallest of CW=Forward and ACW=Backward results
  result <- list(4)
  if (resultACW[[1]] < resultCW[[1]]){
    result <- resultACW
    result[[4]] <- 0 #Clock wise = Forward/Backward flag
  }else{
    
    result <- resultCW
    result[[4]] <- 1 #Clock wise flag
    if (resultACW[[1]] == resultCW[[1]]) { result[[4]] <- 2}
  }
  
  
  
  
  return(result)
}

#########################################################################################################










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
#########################################################################################################


#########################################################################################################
#Delaunay utility functions
#Make sure you have global df for delaunay edges
#########################################################################################################

#########################################################################################################
#checks if the edge is Delaunay edge
#########################################################################################################
is.DelEdge <- function(src, dst, DelEdgeDF){
  
  #undirected !!
  i <- which( (DelEdgeDF$src == src & DelEdgeDF$dst == dst) | (DelEdgeDF$src == dst & DelEdgeDF$dst == src) )
  if (length(i) == 0) flag = FALSE
  else flag = TRUE
  
  return(flag)
}

#########################################################################################################
#Counts Delaunay edges given the tour
#########################################################################################################

cntDelEdge <- function(thePerm, DelEdgeDF){
  
  cnt <- 0
  l <- length(thePerm)
  
  for (e in 1: (l-1)) {
    if (is.DelEdge(thePerm[e], thePerm[e+1], DelEdgeDF)) cnt <- cnt + 1
  }
  
  if (is.DelEdge(thePerm[l], thePerm[1], DelEdgeDF)) cnt <- cnt + 1
  
  return(cnt)
}
################################################################################################################################

#########################################################################################################
#Finds percentage of Delaunay edges given the tour and the DelEdge df
#########################################################################################################
my_get_del_percent_fast <- function(thePerm, DelEdgeDF) {

 nEdges <- length(thePerm)
 nDel <- cntDelEdge(thePerm, DelEdgeDF)
 
 return(100 * nDel / nEdges)
}



#########################################################################################################
# Do not repeat the first point at the end like a polygon 
# This function will take care of it!!!!
#Initialize globals before calling this function
#It counts the what percentage of the edges on the tour comes from the Delaunay triangulation
#########################################################################################################


my_get_del_percent <- function(thePerm) {
  
  
  
  dxy <- deldir(theTSPcoordMTX[,1], theTSPcoordMTX[,2])
  
  #DEBUG
  # plot(theTSPObject, theOptimTour, tour_lty = 1, tour_col = "green")
  # text(theTSPcoordMTX[,1], theTSPcoordMTX[,2], cex = 0.5, col = "red", pos = 2)
  # plot(dxy, wlines = "triang", wpoints = "none", lty=2, lwd=0.2, add=T)
  
  
  #Extract the Delaunay edges, they are undirected!!!
  edgesDelaunayDF <- data.frame(
    src = dxy$delsgs$ind1,
    dst = dxy$delsgs$ind2)
  
  
  
  
  #Find the edge costs, dists, and add them to df
  
  nEdges <- nrow(edgesDelaunayDF)
  
  distVec <- NULL
  
  #Create dist matrix out of Delaunay Edges
  nVertex <- nrow(theTSPcoordDF)
  theDelDistMTX <- matrix(Inf, nrow = nVertex, ncol = nVertex)
  
  for (i in (1:nEdges) ){
    
    #DEBUG
    #i <- 2
    
    srcVtxno <- edgesDelaunayDF$src[i]
    dstVtxno <- edgesDelaunayDF$dst[i]
    srcCoord <- theTSPcoordMTX[which(theTSPcoordMTX[,3] == srcVtxno),]
    dstCoord <- theTSPcoordMTX[which(theTSPcoordMTX[,3] == dstVtxno),]
    d <- sqrt( ((dstCoord[1] - srcCoord[1]) * (dstCoord[1] - srcCoord[1])) + ((dstCoord[2] - srcCoord[2]) * (dstCoord[2] - srcCoord[2])))
    distVec <- c(distVec, d)
    
    #Save data to distmatx for Delaunay dist mtx
    theDelDistMTX[srcVtxno, dstVtxno] <- d
    theDelDistMTX[dstVtxno, srcVtxno] <- d
  }
  
  edgesDelaunayDF <- cbind(edgesDelaunayDF, distVec)
  colnames(edgesDelaunayDF) <- c("src", "dst", "cost")
  
  
  
  tl <- length(thePerm)
  nDelEdge <- 0
  for (i in 1:(tl-1)) {
    
    s <- thePerm[i]
    d <- thePerm[i+1]
    
    if ( (theDelDistMTX[s,d] != Inf) | (theDelDistMTX[d,s] != Inf)) nDelEdge <- nDelEdge +1
  }
  
  #The last edge?
  
  s <- thePerm[tl]
  d <- thePerm[1]
  
  if ( (theDelDistMTX[s,d] != Inf) | (theDelDistMTX[d,s] != Inf)) nDelEdge <- nDelEdge +1
  
  #DEBUG
  #cat("The nn opt tour has", tl, "edges.", nDelEdge, "of them are from Delaunay Edges. Percentage:", 100 * nDelEdge / tl)
  
  return(100 * nDelEdge / tl)
}

#########################################################################################################














