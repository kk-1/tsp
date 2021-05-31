########################################################################################################
#Clean the environment
remove(list = ls())
########################################################################################################
myWorkDir <- "/my_hd1/my_dir/my_prg/Camerino/tsp" #Set the working dir
setwd(myWorkDir)
########################################################################
# Create lattices of TSPLIB data set
########################################################################
library(TSP)
library(sp)
library(sf)



#########################################################################
#Set Global vars
#########################################################################
#N pts in x and y axis
nx <- c(10, 10, 15, 20, 20) #10,  10,  15,  20,  20
ny <- c(10, 20, 20, 20, 25) #10,  20,  20,  20,  25

#Separation between each point in x and y axis
ystep <- 10
xstep <- 10


#Total N pts for uniformly rnd pts
nRND <- c(100, 200, 300, 400, 500) #100, 200, 300, 400, 500

#The width and height of the field for RND pts
nxRND <-  c(1000, 2000, 3000, 4000, 5000) # 1000, 2000, 3000, 4000, 5000
nyRND <-  c(1000, 2000, 3000, 4000, 5000) # 1000, 2000, 3000, 4000, 5000

marginRND <- 20 # The margin that ensures rnd points are that apart from each other for the RND pts
stepRND <- 10 #For rnd points this is the "grid" resolution


#Removal percentage for irregular grid
removePercent <- 0.15 # 0.05,  0.10, 0.15



##########################################################################################################################################
#BEGIN --- IRREGULAR GRIDS
##########################################################################################################################################

nGrid <- length(nx)
for (i in 1:nGrid) {
  
  
  ##########################################################################################################################################
  #########################################################################
  #Create irregular rectangular grid, lattice with rnd pts removed
  #########################################################################
  
  
  #Removal parameter for each axis
  
  removeLenX <- ceiling(nx[i] * removePercent)
  removeLenY <- ceiling(ny[i] * removePercent)
  
  nx2 <- nx[i] + removeLenX
  ny2 <- ny[i] + removeLenY 
  
  
  
  xcoord <- seq(0, ((nx2-1) * xstep), by=xstep)
  
  ycoord <- seq(0, ((ny2-1) * ystep), by=ystep)
  
  tspDF <- as.data.frame(expand.grid(xcoord, ycoord))
  
  
  #remove rnd pts
  tspIDX <- 1:nrow(tspDF)
  removeIDX <- sample(tspIDX, ((removeLenY * nx[i]) + (removeLenX * ny[i]) + (removeLenX * removeLenY)), replace = FALSE, prob = NULL)
  tspDF2 <- tspDF[-removeIDX,]
  cat(nrow(tspDF),"points generated for",nx2,"by",ny2,"lattice.\n")
  cat(length(removeIDX),"rnd points are removed from regular lattice.\n")
  
  colnames(tspDF2) <- c("x", "y")
  cat(nrow(tspDF2),"points selected for",nx2,"by",ny2,"lattice.\n")
  
  tspFileName <- paste0("myRNDLattice-",nx2,"x",ny2,"-",nrow(tspDF2),".tsp")
  
  etsp <- ETSP(tspDF2[,1:2]) 
  write_TSPLIB(etsp, file=tspFileName)
  cat(tspFileName,"is saved!\n\n")
  
  
  #For EPS output
  cairo_ps(filename = paste0("myRNDLattice-",nx2,"x",ny2,"-",nrow(tspDF2),".eps"),
           width = 7,
           height = 7,
           pointsize = 12,
           fallback_resolution = 600)
  
  
  plot(tspDF2, 
       xlab = "x-coordinate", 
       ylab = "y-coordinate", 
       main = paste0("myRNDLattice-",nx2,"x",ny2,"-",nrow(tspDF2),".tsp"), 
       asp=1,
       col.axis="blue", 
       font.axis=4, 
       cex.axis=1.0
  )
  
  dev.off()
  
  
  #DEBUG plot
  
  plot(tspDF2, 
       xlab = "x-coordinate", 
       ylab = "y-coordinate", 
       main = paste0("myRNDLattice-",nx2,"x",ny2,"-",nrow(tspDF2),".tsp"), 
       asp=1,
       col.axis="blue", 
       font.axis=4, 
       cex.axis=1.0
  )
  ##########################################################################################################################################
  
  
  
  
  
  
  
  
  
  ##########################################################################################################################################
  ####################################################################################
  #Create irregular Hexagonal lattice by using sf library
  #All the point sampling functions in R, they all suck
  #Determine your bbox, sample and pray that you well have more than enough pts
  #Then throw away the ones outside of your bbox!!
  #################################################################################
  
  #Removal parameter for each axis
  
  removeLenX <- ceiling(nx[i] * removePercent)
  removeLenY <- ceiling(ny[i] * removePercent)
  
  nx2 <- nx[i] + removeLenX
  ny2 <- ny[i] + removeLenY 
  
  
  myCellSize <- 10
  
  n <- nx2 * ny2
  
  #Some hocus pocus calcs
  a <- myCellSize / (sqrt(3))
  regionx <- nx2 * myCellSize + (myCellSize + 1)
  regiony <- (1.5 * ny2) * a + (ny2/nx2 * myCellSize + 1)
  
  sfc = st_sfc(st_polygon(list(rbind(c(0,0), c(regionx,0), c(regionx,regiony), c(0,regiony), c(0,0)))))
  #sfc <- sfc + c(0, 2*myCellSize)
  #DEBUG plot
  plot(sfc, 
       xlab = "x-coordinate", 
       ylab = "y-coordinate", 
       main = paste0("RNDHexLattice with", n, "pts"), 
       asp=1,
       col.axis="blue", 
       font.axis=4, 
       cex.axis=1.0,
       xlim=c(0, regionx+3*myCellSize), 
       ylim=c(0, regiony+3*myCellSize)
  )
  
  axis(2, at=seq(0, regiony+3*myCellSize, by=myCellSize), col.axis="red", las=2)
  axis(1, at=seq(0, regionx+3*myCellSize, by=myCellSize), col.axis="red", las=2)
  
  
  hexCenters <- st_make_grid(
    sfc,
    #cellsize = myCellSize, 
    n = c(nx2, ny2), 
    what = "centers",
    square = FALSE,
    offset = c(0, 0)
  )
  length(hexCenters)
  hexCenters <- hexCenters + c(myCellSize,myCellSize)
  #DEBUG plot
  plot(hexCenters, add = TRUE)
  
  hexPolys <- st_make_grid(
    sfc,
    #cellsize = myCellSize, 
    n = c(nx2, ny2), 
    what = "polygons",
    square = FALSE,
    offset = c(0,0)
  )
  hexPolys <- hexPolys + c(myCellSize,myCellSize)
  #DEBUG plot
  plot(hexPolys, add = TRUE)
  
  
  
  
  #################################################################################
  #Extract coordinates and write them to tsp file
  #################################################################################
  latticeDF <- as.data.frame(st_coordinates(hexCenters))
  latticeDF <- latticeDF[which(latticeDF$X <= regionx & latticeDF$Y <= regiony),]
  
  
  #remove rnd pts
  latticeIDX <- 1:nrow(latticeDF)
  removeIDX <- sample(latticeIDX, ((removeLenY * nx[i]) + (removeLenX * ny[i]) + (removeLenX * removeLenY)), replace = FALSE, prob = NULL)
  latticeDF2 <- latticeDF[-removeIDX,]
  cat(nrow(latticeDF),"points generated for",nx2,"by",ny2,"lattice.\n")
  cat(length(removeIDX),"rnd points are removed from regular lattice.\n")
  cat(nrow(latticeDF2),"points selected for",nx2,"by",ny2,"lattice.\n")
  latticeMTX <- as.matrix(cbind(as.numeric(latticeDF2$X), as.numeric(latticeDF2$Y)))
  latticeMTX <- round(latticeMTX, digits = 3)
  
  
  
  
  #For EPS output
  cairo_ps(filename = paste0("myRNDHexLattice-",nx2,"x",ny2,"-",nrow(latticeDF2),".eps"), 
           width = 7,
           height = 7,
           pointsize = 12,
           fallback_resolution = 600)
  
  
  
  plot(latticeMTX, 
       xlab = "x-coordinate", 
       ylab = "y-coordinate", 
       main = paste0("myRNDHexLattice-",nx2,"x",ny2,"-",nrow(latticeDF2),".tsp"), 
       asp=1,
       col.axis="blue", 
       font.axis=4, 
       cex.axis=1.0,
       xlim=c(0, regionx+3*myCellSize), 
       ylim=c(0, regiony+3*myCellSize)
       
  )
  
  
  dev.off()
  
  
  
  
  #DEBUG plot
  plot(latticeMTX, 
       xlab = "x-coordinate", 
       ylab = "y-coordinate", 
       main = paste0("myRNDHexLattice-",nx2,"x",ny2,"-",nrow(latticeDF2),".tsp"), 
       asp=1,
       col.axis="blue", 
       font.axis=4, 
       cex.axis=1.0,
       xlim=c(0, regionx+3*myCellSize), 
       ylim=c(0, regiony+3*myCellSize)
       
  )
  sfc <- sfc + c(0, myCellSize)
  plot(sfc, asp=1, add=TRUE)
  nPts <- nrow(latticeMTX)
  cat(nPts,"points generated for",nx2,"by",ny2,"lattice.\n")
  
  
  
  
  
  
  
  #################################################################################
  #Save data to TSPLIB file
  #################################################################################
  tspFileName <- paste0("myRNDHexLattice-",nx2,"x",ny2,"-",nrow(latticeDF2),".tsp")
  etsp <- ETSP(latticeMTX)
  write_TSPLIB(etsp, file=tspFileName)
  cat(tspFileName,"is saved!\n\n")
  ##########################################################################################################################################
  
  
}

##########################################################################################################################################
#END --- IRREGULAR GRIDS
##########################################################################################################################################






























##########################################################################################################################################
#BEGIN --- REGULAR GRIDS
##########################################################################################################################################

nGrid <- length(nx)
for (i in 1:nGrid) {
  
  ##########################################################################################################################################
  #########################################################################
  #Create Rectangular grid, lattice
  #########################################################################
  
  
  xcoord <- seq(0, ((nx[i]-1) * xstep), by=xstep)
  
  ycoord <- seq(0, ((ny[i]-1) * ystep), by=ystep)
  
  tspDF <- as.data.frame(expand.grid(xcoord, ycoord))
  colnames(tspDF) <- c("x", "y")
  cat(nrow(tspDF),"points generated for",nx[i],"by",ny[i],"lattice.\n")
  
  tspFileName <- paste0("myLattice-",nx[i],"x",ny[i],"-",nrow(tspDF),".tsp")
  
  etsp <- ETSP(tspDF[,1:2]) 
  write_TSPLIB(etsp, file=tspFileName)
  cat(tspFileName,"is saved!\n\n")
  
  #For EPS output
  cairo_ps(filename = paste0("myLattice-",nx[i],"x",ny[i],"-",nrow(tspDF),".eps"),
           width = 7,
           height = 7,
           pointsize = 12,
           fallback_resolution = 600)
  
  plot(tspDF,
       xlab = "x-coordinate", 
       ylab = "y-coordinate", 
       main = paste0("myLattice-",nx[i],"x",ny[i],"-",nrow(tspDF),".tsp"), 
       asp=1,
       col.axis="blue", 
       font.axis=4, 
       cex.axis=1.0,
  )
  
  dev.off()
  
  
  #DEBUG plot
  plot(tspDF,
       xlab = "x-coordinate", 
       ylab = "y-coordinate", 
       main = paste0("myLattice-",nx[i],"x",ny[i],"-",nrow(tspDF),".tsp"), 
       asp=1,
       col.axis="blue", 
       font.axis=4, 
       cex.axis=1.0,
  )
  ##########################################################################################################################################
  
  
  
  
  
  
  
  
  ##########################################################################################################################################
  ####################################################################################
  #Create Hexagonal lattice by using sf library
  #All the point sampling functions in R, they all suck
  #Determine your bbox, sample and pray that you well have more than enough pts
  #Then throw away the ones outside of your bbox!!
  #################################################################################
  
  
  myCellSize <- 10
  
  n <- nx[i] * ny[i]
  
  #Some hocus pocus calcs
  a <- myCellSize / (sqrt(3))
  regionx <- nx[i] * myCellSize + (myCellSize + 1)
  regiony <- (1.5 * ny[i]) * a + (ny[i]/nx[i] * myCellSize + 1)
  
  sfc = st_sfc(st_polygon(list(rbind(c(0,0), c(regionx,0), c(regionx,regiony), c(0,regiony), c(0,0)))))
  #sfc <- sfc + c(0, 2*myCellSize)
  
  #DEBUG plot
  plot(sfc, 
       xlab = "x-coordinate", 
       ylab = "y-coordinate", 
       main = paste0("HexLattice with", n, "pts"), 
       asp=1,
       col.axis="blue", 
       font.axis=4, 
       cex.axis=1.0,
       xlim=c(0, regionx+3*myCellSize), 
       ylim=c(0, regiony+3*myCellSize)
  )
  
  axis(2, at=seq(0, regiony+3*myCellSize, by=myCellSize), col.axis="red", las=2)
  axis(1, at=seq(0, regionx+3*myCellSize, by=myCellSize), col.axis="red", las=2)
  
  
  hexCenters <- st_make_grid(
    sfc,
    #cellsize = myCellSize, 
    n = c(nx[i], ny[i]), 
    what = "centers",
    square = FALSE,
    offset = c(0, 0)
  )
  length(hexCenters)
  hexCenters <- hexCenters + c(myCellSize,myCellSize)
  #DEBUG plot
  plot(hexCenters, add = TRUE)
  
  hexPolys <- st_make_grid(
    sfc,
    #cellsize = myCellSize, 
    n = c(nx[i], ny[i]), 
    what = "polygons",
    square = FALSE,
    offset = c(0,0)
  )
  hexPolys <- hexPolys + c(myCellSize,myCellSize)
  #DEBUG plot
  plot(hexPolys, add = TRUE)
  
  
  
  
  #################################################################################
  #Extract coordinates and write them to tsp file
  #################################################################################
  latticeDF <- as.data.frame(st_coordinates(hexCenters))
  latticeDF <- latticeDF[which(latticeDF$X <= regionx & latticeDF$Y <= regiony),]
  latticeMTX <- as.matrix(cbind(as.numeric(latticeDF$X), as.numeric(latticeDF$Y)))
  latticeMTX <- round(latticeMTX, digits = 3)
  
  
  
  
  
  #For EPS output
  cairo_ps(filename = paste0("myHexLattice-",nx[i],"x",ny[i],"-",nrow(latticeMTX),".eps"),
           width = 7,
           height = 7,
           pointsize = 12,
           fallback_resolution = 600)
  
  
  
  plot(latticeMTX, 
       xlab = "x-coordinate", 
       ylab = "y-coordinate", 
       main = paste0("myHexLattice-",nx[i],"x",ny[i],"-",nrow(latticeMTX),".tsp"), 
       asp=1,
       col.axis="blue", 
       font.axis=4, 
       cex.axis=1.0,
       xlim=c(0, regionx+3*myCellSize), 
       ylim=c(0, regiony+3*myCellSize)
       
  )
  
  
  dev.off()
  
  
  #DEBUG plot
  
  plot(latticeMTX, 
       xlab = "x-coordinate", 
       ylab = "y-coordinate", 
       main = paste0("myHexLattice-",nx[i],"x",ny[i],"-",nrow(latticeMTX),".tsp"), 
       asp=1,
       col.axis="blue", 
       font.axis=4, 
       cex.axis=1.0,
       xlim=c(0, regionx+3*myCellSize), 
       ylim=c(0, regiony+3*myCellSize)
       
  )
  
  sfc <- sfc + c(0, myCellSize)
  plot(sfc, asp=1, add=TRUE)
  nPts <- nrow(latticeMTX)
  cat(nPts,"points generated for",nx[i],"by",ny[i],"lattice.\n")
  
  
  
  
  #################################################################################
  #Save data to TSPLIB file
  #################################################################################
  tspFileName <- paste0("myHexLattice-",nx[i],"x",ny[i],"-",nrow(latticeMTX),".tsp")
  etsp <- ETSP(latticeMTX)
  write_TSPLIB(etsp, file=tspFileName)
  cat(tspFileName,"is saved!\n\n")
  ##########################################################################################################################################
  
  
}



nGrid <- length(nRND)
for (i in 1:nGrid) {
  
  
  ##########################################################################################################################################
  #########################################################################
  #Create RND uniform pts in the square
  #########################################################################
  #DEBUG   i <- 4
  
  
  x_coord <- seq(0, nxRND[i], by=stepRND)
  y_coord <- seq(0, nyRND[i], by=stepRND)
  tspDF <- as.data.frame(expand.grid(x_coord, y_coord))
  tspMTX <- as.matrix(tspDF)
  colnames(tspDF) <- c("x","y")
  
  ptsMTX <- matrix(0, nrow = nRND[i], ncol = 2)
  
  #Create rnd points marginRND apart
  for (k in 1:nRND[i]){
    #DEBUG k <- 2
    #Pick a point
    ptsX <- sample(tspMTX[,1], 1)
    ptsY <- sample(tspMTX[,2], 1)
    ptsMTX[k,] <- c(ptsX, ptsY)
    #Erase the neighborhood so that there will be no point close the ones selected
    eraseIDX <- which( ((tspMTX[,1] <= (ptsX + marginRND)) & (tspMTX[,1] > (ptsX - marginRND)) & 
                         (tspMTX[,2] <= (ptsY + marginRND)) & (tspMTX[,2] > (ptsY - marginRND)))  |
                         ((tspMTX[,1] == ptsX) & (tspMTX[,2] == ptsY))
                       )
    
    #DEBUG
    #bf <- nrow(tspMTX)
    #cat(k,"-",bf," pts before erase")
    
    if (length(eraseIDX) > 0) tspMTX <- tspMTX[-eraseIDX, ]
    
    #DEBUG
    #af <-  nrow(tspMTX)
    #cat(" -",af," pts after erase.",(bf -af),"pts deleted!\n")
  }
  
  ptsDF <- as.matrix(ptsMTX)
  colnames(ptsDF) <- c("x","y")
  
  
  #DEBUG plot
  plot(ptsDF)
  
  
  
  nPts <- nrow(ptsMTX)
  
  
  #DEBUG plot
  plot(ptsMTX, 
       xlab = "x-coordinate", 
       ylab = "y-coordinate", 
       main = paste0("myRND-",nPts,".tsp"),
       asp=1,
       col.axis="blue", 
       font.axis=4, 
       cex.axis=1.0
       
  )
  
  
  
  cat(nPts,"points generated for",nxRND[i],"by",nyRND[i],"lattice.\n")
  
  
  #Save the tsp file
  tspFileName <- paste0("myRND-",nPts,".tsp")
  etsp <- ETSP(ptsMTX)
  write_TSPLIB(etsp, file=tspFileName)
  cat(tspFileName,"is saved!\n\n")
  
  
  
  
  #For EPS output
  cairo_ps(filename = paste0("myRND-",nPts,".eps"),
           width = 7,
           height = 7,
           pointsize = 12,
           fallback_resolution = 600)
  
  
  
  plot(ptsMTX, 
       xlab = "x-coordinate", 
       ylab = "y-coordinate", 
       main = paste0("myRND-",nPts,".tsp"),
       asp=1,
       col.axis="blue", 
       font.axis=4, 
       cex.axis=1.0
       
  )
  
  dev.off()
  
  
  
  
  
  #########################################################################
  
  ##########################################################################################################################################
  
}


##########################################################################################################################################
#END --- REGULAR GRIDS
##########################################################################################################################################




