##----------------libraries--------------------
##Install packages if necessary.
# install.packages("plyr")
# install.packages("sp")
# install.packages("maptools")
# install.packages("raster")
# install.packages("grDevices")
# install.packages("rgdal")
# install.packages("rgeos")
# install.packages("gdata")
# install.packages("reshape")
# install.packages("cleangeo")
# install.packages("zoom") #not needed for model, but usefull for examining plots

library(plyr)
library(sp)
library(maptools)
library(raster)
library(grDevices)
library(rgdal)
library(rgeos)
library(gdata)
library(reshape)
library(cleangeo)
library(zoom)

#----SET UP PARAMETERS-------------------------------

#mean dispersal distance
cmd <- 20 #for claustral queens
pmd <- 30 #for parasitic queens
#dispersal standard deviation
csd <- 2 #for claustral queens
psd <- 4 #for parasitic queens
#dispersal standard deviation
cmort <- 0.95 #random mortality for claustral queen dispersal
pmort <- 0.95 #random mortality for parasitic queen dispersal

gf <- 0.2 #growth factor (fraction of original size a colony can grow)

qp <- 0.1 #queen production ratio, proportion of queens produced per worker

cqd <- c(0, 0, 0, 0.2, 0.2, 0.3, 0.2, 0.1, 0, 0, 0, 0)  #Vector of ratios for how many claustral queens disperse each month
pqd <- c(0, 0.5, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0)      #Vector of ratios for how many claustral queens disperse each month

xmax = 50 #how broad is the arena
ymax = 1000 #how tall is the arena

arena <- extent(c(0, xmax, 0, ymax))

gridres = 3 #how many decimal places for rounding polygon coordinates (3 = millimeter resoluton)

gridcells = 4 #how many cells on a side for dividing up overlapping growth polygons

zsurv <- 6 #How many months does a zombied colony survive without a queen?

offenders <- NULL #global variable for dealing with bad topologies

#-----FUNCTIONS----------------------------

disp2 <- function(df, md, sd, xmax) { #dataframe dispersal function: takes a dataframe with x and y coordinates and mean and sd of distribution function
  gpar <-  drop(cbind(shape=(md/sd)^2, rate=md/sd^2)) #get the gamma parameters - function borrowed from SGAT package
  df$dist <- rgamma(nrow(df), gpar[1], gpar[2])
  df$ang <- runif(nrow(df), 0, 2*pi)
  df$x <- df$x + df$dist * cos(df$ang)
  df$y <- df$y + df$dist * sin(df$ang)
  df$x <- df$x %% xmax #get the modulus (remainder) of the x value so it wraps around horizontally if necessary
  drops <- c("dist","ang")
  df <- df[,!(names(df) %in% drops)] #get rid of unwanted columns
  df <- df[c(df$y > 0 & df$y < ymax),] #get rid of colonies that fall beyond y boundaries.
  return(df)
} 

makeSquares <- function(points){ #make a spatial polygons data frame with squares around points 
  #points must be three columns with x and y coordinates and area
  #define corners of the box based on center and area 
  points$x1 <- points$x - points$area^(.5)/2
  points$x2 <- points$x + points$area^(.5)/2
  points$y1 <- points$y - points$area^(.5)/2
  points$y2 <- points$y + points$area^(.5)/2
  #make   polygons from points 
  boxes <- mapply(function(x1,x2,y1,y2,ID) Polygons(list(Polygon(matrix(c(x1,x1,x2,x2,x1,y1,y2,y2,y1,y1), ncol = 2))), ID), points$x1, points$x2, points$y1, points$y2, points$ID)
  boxes <- SpatialPolygons(boxes)
  boxes <- RoundPolygons(boxes)
  return(boxes)
}

af <- function(wrkr){ #expected area based on worker population
  area <- 0.00064*wrkr # 6.4 square cenitimeters per worker from Tschinkel et al 1995.
  return(area) 
}

expWrkr <- function(area) { #gives number of workers expected for a given colony area
  #uses same equation as af - just solved for workers
  wrkr <- area/0.00064
  return(round(wrkr,0)) #round answer to nearest whole number
}

RoundPolygons<-function(shptemp, digitss=gridres) {
  for(i in 1:length(shptemp)) {
    shptemp@polygons[[i]]@Polygons[[1]]@coords<-round(shptemp@polygons[[i]]@Polygons[[1]]@coords,digits=digitss)
  }
  return(shptemp)
}

CatchTopErr <- function(poly, name) { #catch and correct topology errors
  cc<-try(gIsValid(poly), silent=T) #test for topology problems
  if(cc != T) {
    print(paste("Topology error detected in ", name, ". Correcting with cleangeo...", sep = ""))
    testpoly <- try(clgeo_Clean(poly), silent=T)
    if(!is.null(nrow(testpoly))) {poly <- testpoly}
  }
  cc<-try(gIsValid(poly), silent=T) #test for topology problems
  if(cc != T) {
    print(paste("Error persists in ", name, ". Correcting with FindOffenders...", sep = ""))
    offenders <- NULL
    offenders <- FindOffenders(poly)
    if(!is.null(offenders)) {  
      cat("\nRemoving polygons:", offenders, "\n")
      poly <- poly[-offenders,]
      offenders <- NULL
    } 
  }
  cc<-try(gIsValid(poly), silent=T) #test for topology problems
  if(cc != T) {
    print(paste("Error persists in ", name, ". Cannot fix it...", sep = ""))
  }
  return(poly)
}

CleanTop <- function(poly) {
  poly <- gBuffer(poly,  byid = T, width = .01, joinStyle="MITRE", mitreLimit=10) #attempt to clean topology
  poly <- gBuffer(poly,  byid = T, width = -.01, joinStyle="MITRE", mitreLimit=10) 
  return(poly)
}

FindOffenders <- function(sp) { #Sequentially add subpolygons to a polygons object to find and eliminate the ones with topology problems
  plist <- NULL
  offend <- NULL
  incr <- 30
  for(pc in seq(0, length(sp), by = incr)) {
    minp <- pc+1
    maxp <- min(pc+incr, length(sp))
    cat('\r', "\nchecking ", minp, " to ", maxp, sep = ""); flush.console()  
    sp1 <- sp[1:maxp,] #Make a subset of polygons
    if(is.null(offend)) { #remove known offending polygons if there are any
      sp2 <- sp1
    } else {
      sp2 <- sp1[-offend,]
    }
    cc<-try(gIsValid(sp2), silent=T) #test for topology problems
    if(cc != T) {  #if there are problems with the batch of polygons then figure out which particular polygon is the offender
      cat("\n  ") #new line
      for(p in minp:maxp) {
        cat('\r', "checking ", p, sep = ""); flush.console()  
        sp1 <- sp[1:p,] #Make a subset of polygons
        if(is.null(offend)) { #remove known offending polygons if there are any
          sp2 <- sp1
        } else {
          sp2 <- sp1[-offend,]
        }
        cc<-try(gIsValid(sp2), silent=T) #test for topology problems
        if(cc != T) {  #if there are problems with the batch of polygons then figure out which particular polygon is the offender
          cat("\nPolygon ", p, " is an offender.\n", sep = "")
          offend <- append(offend, p)
          cc<-try(gIsValid(sp[-offend,]), silent=T) #test for topology problems
          if(cc == T) {  #if there are no problems then all offenders have been found
            return(offend)
          }
        }
      }
    }
  }
  return(offend)
}

FindOffendersU2 <- function(sp) { #Sequentially add subpolygons to a polygons object to find and eliminate the ones with topology problems
  plist <- NULL
  offend <- NULL
  incr <- 30
  for(pc in seq(0, length(sp), by = incr)) {
    minp <- pc+1
    maxp <- min(pc+incr, length(sp))
    cat('\r', "\nchecking ", minp, " to ", maxp, sep = ""); flush.console()  
    sp1 <- sp[1:maxp,] #Make a subset of polygons
    if(is.null(offend)) { #remove known offending polygons if there are any
      sp2 <- sp1
    } else {
      sp2 <- sp1[-offend,]
    }
    usp <- tryCatch(gUnaryUnion(sp2), error = function(e) e = NULL)
    if(is.null(usp)) {  #if there are problems with the batch of polygons then figure out which particular polygon is the offender
      cat("\n  ") #new line
      for(p in minp:maxp) {
        cat('\r', "checking ", p, sep = ""); flush.console()  
        sp1 <- sp[1:p,] #Make a subset of polygons
        if(is.null(offend)) { #remove known offending polygons if there are any
          sp2 <- sp1
        } else {
          sp2 <- sp1[-offend,]
        }
        usp <- tryCatch(gUnaryUnion(sp2), error = function(e) e = NULL)
        if(is.null(usp)) {  #if there are problems with the batch of polygons then we now know which polygon is the offender
          cat("\nPolygon ", p, " is an offender.", sep = "")
          offend <- append(offend, p)
          usp <- tryCatch(gUnaryUnion(sp[-offend,]), error = function(e) e = NULL)
          if(!is.null(usp)) {  #if there are no problems then all offenders have been found
            return(offend)
          }
        }
      }
    }
  }
  return(offend)
}

rmTiny <- function(poly) {
  d <- disaggregate(poly)   #...separate them
  d2 <- d[which(gArea(d, byid = T) > 0.01),] #see if any are smaller than 1 square cm
  if(length(d2) > 0) {  #If there are some greater than 1 square cm...
    poly2 <- gUnaryUnion(d2) #...then combine them back into one polygon object
  } else { #if there are no polygons greater than 1cm 
    poly2 <- NULL #skip this loop instance and do not assign the shapes to any colonies
  }
  return(poly2)
}

expnd <- function(spdf){ #expand all colonies according to their growth parameters
  
  cc<-try(gIsValid(spdf), silent=F) #test for topology problems
  if(cc != T) {
    cat("\nPolygon topology problem detected in spdf....Fixing with cleangeo")
    spdf <-  clgeo_Clean(spdf)
  }
  #spdf <- RoundPolygons(spdf)              
  
  #Make sure ID column matches polygon IDs
  spdf$ID <- sapply(slot(spdf, "polygons"), function(x) slot(x, "ID")) #get the IDs
  
  rng <- .5 #fraction of orginal shape extent that becomes the working area when expanding territories 
  #smaller = faster
  #larger = less chance of inaccuracy
  buff <- 0.1
  
  ymax2 <- max(spdf$y) + 10
  
  #Make two big polygons on either side of the area. They get used later to identify landlocked colonies.
  p1 <- matrix(data = c(0,-5,-5,0,0,-5,-5,ymax2,ymax2,-5), nrow=5) #establish coordinaates for the four corners
  p2 <- matrix(data = c(xmax,xmax+5,xmax+5,xmax,xmax,-5,-5,ymax2,ymax2,-5), nrow=5) #establish coordinaates for the four corners
  p1 <- Polygon(p1) #make polygons
  p2 <- Polygon(p2)
  p11 <- Polygons(list(p1), ID = row.names(spdf[1,1])) #make Polygons Object
  p12 <- Polygons(list(p2), ID = row.names(spdf[2,1]))
  sides <- SpatialPolygons(list(p11, p12)) #make Spatial Polygons Object
  sides <- SpatialPolygonsDataFrame(sides, data=(spdf@data[1:2,])) #Assign some dummy data to make a sptial polygons data frame
  row.names(sides) <- c("a", "b") #assign some random row names
  
  #Bind sides to the main spatial polygons data frame
  
  all <-rbind(spdf, sides) 
  
  cc<-try(gIsValid(all), silent=F) #test for unary union problems
  if(cc != T) {
    cat("\nPolygon union problem with all detected....Fixing with cleangeo")
    #all <- RoundPolygons(all)
    allarea <- sapply(all@polygons[[1]]@Polygons, function(x) slot(x, "area"))
    polys <-  all@polygons[[1]]@Polygons[which(allarea > 0.0001)]
    all <- Polygons(polys, ID = "all")
    all <- SpatialPolygons(list(all))
    all <- clgeo_Clean(all)
  }
  
  all2 <- tryCatch(gUnaryUnion(all),
                   error = function(er) {cat("\ngUnaryUnion error when creating all2");
                                         offenders <<- FindOffendersU2(all); 
                                         cat("\nremoving polygon(s):", offenders, "\n");
                                         gUnaryUnion(all[-offenders,]) } )  
  
  if(!is.null(offenders)) {  
    spdf <- spdf[-offenders,]
    offenders <- NULL
  }
  
  #Make a SPDF that excludes any colonies that cannot grow
  excld <- gContainsProperly(all2, spdf, byid = T) #T/F list for landlocked colonies (colonies that cannot grow)
  excld <- data.frame(ID = row.names(excld), tf = excld[,1])
  excld <- excld[excld$tf,] #Just the lines with TRUE (the landlocked colonies)
  if(nrow(excld) == 0) {
    free <- spdf #Do this if all the colonies can grow
  } else {
    free <- spdf[which(!spdf$ID %in% excld$ID),] #These are the colonies that can grow (not landlocked)
  }
  
  #Determine if the population is large enough to warrant growth
  free$area <- gArea(free, byid = T) #Get accurate area measure loaded in
  expSize <- lapply(free$pop, FUN=af) #find the expected size of the colony
  togrow <- which(expSize > free$area)
  if(length(togrow)==0) { #Exit the function if there is no growth needed
    return(spdf)
  } 
  grow <- free[togrow,] #Isolate the colonies that can and need to grow
  
  #Generate buffers around each growing colony
  buff1 <- gBuffer(grow,  byid = T, width = (grow$area^0.5)*gf, joinStyle="MITRE", mitreLimit=10) #define the buffer
  #buff1 <- RoundPolygons(buff1)
  buff1 <- CatchTopErr(buff1, "buff1")
  buff2 <- gDifference(buff1, all2, byid = T) #remove occupied space from the buffer1
  buff2 <- buff2[which(gArea(buff2, byid = T) > 0.001),] #get rid of tiny polygons (less than 10 sq cm)
  #buff2 <- RoundPolygons(buff2)
  
  buff2 <- CatchTopErr(buff2, "buff2")
  
  bIDs <- sapply(slot(buff2, "polygons"), function(x) gsub(' 1','',slot(x, "ID"))) #Get the polygon IDs and remove the useless text
  buff2 <- spChFIDs(buff2, bIDs) #assign the modified IDs to the polygons
  inters <- gIntersects(buff2, spgeom2 = NULL, byid = T) #Find all intersecting polygons
  upperTriangle(inters, diag=T) <- F  #Eliminate the upper triangle of the matrix of the T/F overlap matrix (it's redundant)
  mint <- melt(inters) #melt the matrix to get row and column pairs for each matrix value
  int <- mint[which(mint$value),] #Toss out all of the FALSE values (no overlap): The remaining lines have overlaps that need to be dealt with
  
  #Do this if there are no overlapping areas of growth
  if(nrow(int) == 0) {
    #Get the buffers into a Spatial Polygons Data Frame
    #bdata <- grow@data[which(grow$ID %in% bIDs),] #probably good enough, but may be vulnerable to ID mixups
    bIDs <- as.data.frame(bIDs) #Make the list of buffer IDs a data frame
    names(bIDs) <- c("ID") #Change the column name to match the ID column in grow@data
    bdata <- join(bIDs, grow@data, by = "ID", match = "all")  #create data for SPDF
    bdata <- bdata[,c(2,3,4,5,6,7,8,9,1,10,11,12,13)] #reorder the columns
    bdata$ID <- as.character(bdata$ID)
    row.names(bdata) <- bdata$ID   
    simpleBuffs <- SpatialPolygonsDataFrame(buff2, data = bdata) #Make a SPDF ready to bind to other SPDFs 
    simpleBuffs <- spChFIDs(simpleBuffs, paste("b", bIDs$ID, sep = "")) #Change IDs again to make them unique
    
    #Merge the buffers with their parent polygons and reformat the SPDF
    grown <- rbind(grow, simpleBuffs)
    grown <- unionSpatialPolygons(grown, IDs = grown$ID)
    gIDs <- sapply(slot(grown, "polygons"), function(x) slot(x, "ID")) #Get the polygon IDs
    gIDs <- as.data.frame(gIDs)
    names(gIDs) <- c("ID") #Change the column name to match the ID column in grow@data
    gdata <- join(gIDs, grow@data, by = "ID", match = "all")  #create data for SPDF
    gdata <- gdata[,c(2,3,4,5,6,7,8,9,1,10,11,12,13)] #reorder the columns
    gdata$ID <- as.character(gdata$ID)
    row.names(gdata) <- gdata$ID 
    grown <- SpatialPolygonsDataFrame(grown, data = gdata) #Make a SPDF ready to bind to other SPDFs 
    
    #Replace the old parent polygons with the augmented ones.
    if(sum(spdf$ID %in% grown$ID) == length(spdf)) { #do this if all of the colonies have grown
      spdf2 <- grown
    } else { #otherwise bind the colonies that have grown to the ones that have not
      spdf2 <- spdf[which(!spdf$ID %in% grown$ID),] #remove the polygons to be replaced
      spdf2 <- rbind(spdf2, grown) #add in the altered colonies
    }
    spdf2$area <- sapply(slot(spdf2, "polygons"), slot, "area") #update the areas
    spdf2$ID <- as.numeric(spdf2$ID) #Make IDs numeric before return
    spdf2 <- RoundPolygons(spdf2)
    return(spdf2)
  }
  
  #Continue here if we have to deal with overlapping growth areas
  buff3 <- buff2 #temp storage -- working polygons object
  bID <- bIDs #temp storage -- working polygons object
  parentIDs <- NULL #vector for parent polygon IDs -- will become a data frame later
  fragmentIDs <- NULL #vector for fragment IDs -- will become a data frame later 
  
  gridall <- NULL  #set gridall to zero so it will be evident if it remains empty
  
  #loop through each overlapping growth zone
  for(i in 1:nrow(int)) {
    #for(i in 1:25) {
    #Make a polygon of the overlapp
    P1ind <- int$X1[i]    #ID of parent Polygon 1
    P2ind <- int$X2[i]    #ID of parent Polygon 2
    
    buff3_1 <- buff3[which(bID==P1ind),]
    buff3_2 <- buff3[which(bID==P2ind),]
    
    if(length(buff3_2) * length(buff3_1) == 0) {next}
    
    if(length(buff3_1@polygons[[1]]@Polygons) > 1) { 
      buff3_1 <- rmTiny(buff3_1)
    }
    if(length(buff3_2@polygons[[1]]@Polygons) > 1) { 
      buff3_2 <- rmTiny(buff3_2)
    }
    
    if(is.null(buff3_1)) {next}
    if(is.null(buff3_2)) {next}
    
    Pint <-  gIntersection(buff3_1, buff3_2, drop_lower_td = T) #This is the overlap that needs to be dealt with
    if(is.null(Pint)) {next} #abort the loop instance is Pint is null (the area was already consumed in a previous loop instance)
    if(gArea(Pint) < 0.0001) {next} #If the area is less than 1 square cm, then forget about it. Skip loop instance and do not allocate the space. This will avert some errors and save time
    #Pint <- RoundPolygons(Pint)
    tryval <- try(gIsValid(Pint), silent==T)
    test <- CatchTopErr(Pint, "Pint")
    if(tryval != T) {next} #skip growth if the geometry is not valid
    #remove subpolygons that are too small to matter
    if(length(Pint@polygons[[1]]@Polygons) > 1) {  #If you have multiple sub polygons...
      Pint <- rmTiny(Pint)
    }
    if(is.null(Pint)) {next}
    
    #plot(Pint, add = T, col = "red")
    
    #Remove the overlaping area from the Buffer polygons object. This avoids problems where three or more buffers overlap
    buff3 <- CatchTopErr(buff3, "buff3")
    buff3 <- gDifference(buff3, Pint, byid = T) #overlap removed from working polygons object - This avoids problems where three or more buffers overlap
    bID <- sapply(slot(buff3, "polygons"), function(x) gsub(' 1','',slot(x, "ID"))) #edit the Polygon IDs to remove useless text
    buff3 <- spChFIDs(buff3, bID) #restore IDs that match parent IDs
    
    #Divide the overlapping bit into a grid of cells
    ext <- extent(Pint)
    grid <- raster(ext, nrows = gridcells, ncols = gridcells)
    grid <- rasterToPolygons(grid)
    #grid <- RoundPolygons(grid)
    grid <- gIntersection(grid, Pint, byid = T, drop_lower_td = T) #crop grid to relevant area
    #grid <- RoundPolygons(grid)
    bigenuf <- which(gArea(grid, byid = T) > 0.0001) #identify polygons grid cells big enough to consider
    if(length(bigenuf) == 0) { next } #If no polygons were big enough go to the next loop instance
    grid <- grid[bigenuf,] #Only keep fragments that are bigger than 1 cm squared
    
    #Give the grid cells IDs and add them to the general collection of grid cells in gridall
    maxID <- ifelse(is.null(fragmentIDs), 0, as.numeric(max(fragmentIDs)))  #Find the maximum ID number from past grid cell generation
    newIDs <- (maxID+1) : (maxID+length(grid)) #Forge some new ID numbers 
    newIDs <- sprintf("%08d", newIDs) #Make the ID a character with lots of leading zeros
    grid <- spChFIDs(grid, newIDs) #assign new and unique IDs to grid
    if(is.null(gridall)) {
      gridall <- grid #define gridall for the first time
    } else {
      gridall <- rbind(gridall, grid) #on subsequent iterations append new data to gridall
    }
    
    #Determine where to assign the grid cells 
    centroids <- as.data.frame(coordinates(grid)) #Get the centroid coordinates of each grid cell
    colnames(centroids) <- c("x", "y") #rename the columns
    coordinates(centroids) <- c("x", "y") #make a Spatial Points Object from the centroids
    
    #Generate a dataframe with the grid cell ids and the distances from each grid cell centroid to the centroid of the competing, growing parent colonies
    griddata <- data.frame(ID = sapply(slot(grid, "polygons"), function(i) slot(i, "ID")), #ID numbers of the grid polygons
                           dist1 = spDistsN1(centroids, coordinates(grow[which(grow$ID==P1ind),]), longlat = FALSE), #distances from parent polygon 1
                           dist2 = spDistsN1(centroids, coordinates(grow[which(grow$ID==P2ind),]), longlat = FALSE), #distances from parent polygon 2
                           stringsAsFactors = F)
    
    gridtemp <- griddata #temp storage for griddata
    
    for(j in 1:nrow(gridtemp)) {  #go through all rows in the gridtemp dataframe and assign cells to parent colonies
      if(j %% 2 == 1) {  #alternate between the first parent and second parent
        P1 <- which(gridtemp$dist1 == min(gridtemp$dist1))[1] #find the grid cell with the minimum distance for polygon 1, the extra bracket means only one value will get selected
        Px <- P1ind
      } else {
        P1 <- which(gridtemp$dist2 == min(gridtemp$dist2))[1] #find the grid cell with the minimum distance for polygon 1
        Px <- P2ind
      }
      gridtemp$dist1[P1] <- 99999 #make it so this cell does not have the minimum distance anymore 
      gridtemp$dist2[P1] <- 99999 #make it so this cell does not have the minimum distance anymore
      parentIDs <- append(parentIDs, Px) #add that grid cell ID to the list
      fragmentIDs <- append(fragmentIDs, gridtemp$ID[P1])
    } #end of loop to assign each fragment
    
    #parentIDs #vector for parent polygon IDs -- will become a data frame later
    #fragmentIDs #vector for fragment IDs -- will become a data frame later 
    #gridtemp
  } #end of loop that goes through each overlap
  
  if(is.null(gridall)) {  #Make sure gridall was filled with something
    return(spdf) #nothing actually grew
  }
  
  #Now do a big unification of the assigned fragments and their parent polygons. 
  #first get the fragments into a spatial polygons data frame with the same IDs as their parents
  
  #Prepare the non-overlapping buffers (undisputed expansion area) for merging with parent polygons
  #Buff 3 should only have non-overlapping (non-disputed) expansion areas left. all the overlaps have been removed
  #So just prep Buff3 to merge with grow
  bID <- as.data.frame(bID) #Make the list of buffer IDs a data frame
  names(bID) <- c("ID") #Change the column name to match the ID column in grow@data
  b3data <- join(bID, grow@data, by = "ID", match = "all")  #create data for SPDF
  b3data <- b3data[,c(2,3,4,5,6,7,8,9,1,10,11,12,13)] #reorder the columns
  b3data$ID <- as.character(b3data$ID)
  row.names(b3data) <- b3data$ID   
  buff4 <- SpatialPolygonsDataFrame(buff3, data = b3data) #Make a SPDF ready to bind to other SPDFs 
  buff4 <- spChFIDs(buff4, paste("b", bID$ID, sep = "")) #Change Polygon IDs again to make them unique
  #buff4 <- RoundPolygons(buff4) #causes errors sometimes (WHY??)
  buff4 <- CatchTopErr(buff4, "buff4")
  
  #Prepare the fragments of overlapping buffers for merging with parent polygons
  fragind <- data.frame(ID = parentIDs, frag = fragmentIDs)
  fIDs <- sapply(slot(gridall, "polygons"), function(x) slot(x, "ID")) #Get all of the polygon IDs of the grid cells. They should be in neumerical order
  ord <- order(fragmentIDs) #Get an order permutation from the fragments data frame
  fragind <- fragind[ord,] #Sort the fragments data frame to correspond to the order of the polygons in gridall (should just put them in numerical order) - The parent column is now an index to the ID column in the main Spatial Polygons dataframe.  
  fragdat <- join(fragind, grow@data, by = "ID", match = "all") 
  fragdat <- fragdat[,c(3,4,5,6,7,8,9,10,1,11,12,13,14)] #reorder the columns and drop frag column
  row.names(fragdat) <- fIDs
  #gridall <- RoundPolygons(gridall)
  gridspdf <- SpatialPolygonsDataFrame(gridall, data = fragdat)
  gridspdf$ID <- as.character(gridspdf$ID)
  gridspdf <- CatchTopErr(gridspdf, alist(gridspdf))
  if(try(gIsValid(gridspdf), silent=F) != T) {
    offenders <- NULL
    gridspdf <- FindOffenders(gridspdf) 
    if(!is.null(offenders)) {  
      cat("\nRemoving polygons:", offenders, "\n")
      gridspdf <- gridspdf[-offenders,]
      offenders <- NULL
    }
  }
  
  #Merge the fragments and undisputed buffers with the parent colonies
  grown <- rbind(grow, gridspdf, buff4) #Bind the three SPDFs together (grow=original colonies that need to grow, gridspdf=the divided expansion area, undisputed=undisputed expansion area)
  grown <- CatchTopErr(grown, alist(grown))
  grown <- unionSpatialPolygons(grown, IDs = grown$ID, avoidUnaryUnion=T)
  grown <- CatchTopErr(grown, alist(grown))
  gIDs <- sapply(slot(grown, "polygons"), function(x) slot(x, "ID")) #Get the polygon IDs
  gIDs <- as.data.frame(gIDs)
  names(gIDs) <- c("ID") #Change the column name to match the ID column in grow@data
  gdata <- join(gIDs, grow@data, by = "ID", match = "all")  #create data for SPDF
  gdata <- gdata[,c(2,3,4,5,6,7,8,9,1,10,11,12,13)] #reorder the columns
  gdata$ID <- as.character(gdata$ID)
  row.names(gdata) <- gdata$ID 
  grown <- SpatialPolygonsDataFrame(grown, data = gdata) #Make a SPDF 
  
  #Remove fragments and holes
  for(i in 1:nrow(grown)) {
    if(length(grown@polygons[[i]]@Polygons)>1) { #see if there is more than one polygon (if not no problem)
      areas <- sapply(slot(grown@polygons[[i]], "Polygons"), slot, "area") #Make a vector of all the area fragments
      mxar <- which(areas == max(areas)) #index of largest area
      newpoly <-  grown@polygons[[i]]@Polygons[[mxar]] #Isolate the largest fragment
      newpoly <- Polygons(list(newpoly), ID = grown$ID[i]) #make the largest fragment a Polygons (plural) object
      grown@polygons[[i]] <- newpoly #replace the fragmented polygon with the single largest fragment (smaller fragments go away)
    }
  }
  
  #Replace the original colonies with the grown ones
  if(sum(spdf$ID %in% grown$ID) == length(spdf)) { #do this if all of the colonies have grown
    spdf2 <- grown
  } else { #otherwise bind the colonies that have grown to the ones that have not
    spdf2 <- spdf[!(spdf$ID %in% grown$ID),] #Get the polygons that were not grown
    spdf2 <- rbind(spdf2, grown) #bind the colonies that have grown to the ones that have not
  }
  
  #Kill colonies that are completely engulfed by others
  trapped <- gContainsProperly(spdf2, byid = T) 
  trapped <- melt(trapped) #melt the matrix to get row and column pairs for each matrix value
  trapped <- trapped[which(trapped$value),] #Toss out all of the FALSE values: The remaining lines show trapped colonies
  
  #test <- spdf2[(spdf2$ID %in% trapped$X1),]
  
  spdf3 <- spdf2[!(spdf2$ID %in% trapped$X1),]
  spdf3 <- RoundPolygons(spdf3)
  
  #update area field in data frame
  spdf3$area <- sapply(slot(spdf3, "polygons"), slot, "area") 
  spdf3$ID <- as.numeric(spdf3$ID)
  return(spdf3)
  
} #end of function

pop6 <- function(pop, m, queen, type, age) { #population growth function
  r <- 1.26/12 #growth rate
  mmax <- 1 #month when worker population is greatest)
  #Colony size is calculated based on 3 equations:
  # eq 1: P = 8.41*exp(0.972m) or P = 165000/(1+(83*(exp(-r*(m+1))))) based on current population sizes
  #   P is the unaltered population size before seasonal fluctuations are added
  # eq2: f = cos((2*pi*(m-mmax))/12)*55000 
  #   Calcualtes seasonal oscillaiton as a function of month
  # eq 3: g = 1388407*pop^(-1.1678667) 
  #   Weighting factor for oscillation (small populations are not really affected much)
  
  #to derive P, we must use different functions for colonies at different stages
  #For a new colony, a certain number of workers (less than 100) are generated
  #For a population above 2500, use a logistic growth equation
  #For population below 2500 but above 100, use simple exponential growth
  #For a population without a queen, the growth is negative
  
  if(queen & age == 1 & pop < 2500 & type == "C") { #Claustral queens starting out, test queen present, age, and population size
    return(40) #workers produced during first month
  }
  
  if(queen & pop >= 2500) {
    #for larger populaitons, use a logistic funciton
    #pop must be >= 1965 and <= 165000.
    #Might need a provision for pop exceeding 165000, can that happen?
    if(pop >= 165000) {
      P = 165000 #population cannot exceed this maximum
    } else {
      #Use the time function to get new populaitons size:
      #First determine which time step corresponds to the current population
      t <- (log((165000/pop-1)*1/83))/-r #determine the age (in months) that corresponds with the current populaiton size
      #then calculate P for the next time step
      P = 165000/(1+(83*(exp(-r*(t+1)))))
    }
    #P = r*pop*(1 - pop/165000)
    f = cos((2*pi*(m-mmax))/12)*55000 #Seasonal oscillaiton added to population size
    #g =(1+(150*exp(-1.26*(a/12)))) #Weighting factor for oscillation based on age (should be based on size)
    g = 1388407*pop^(-1.1678667) #Weighting factor for oscillation based on population size (power function approximation)
    return(P + f/g)
  }
  
  if(queen & pop < 2500) { #test queen present, and population size
    #First determine which time step corresponds to the current population
    t = log(pop/8.41)*1/0.972
    #then calculate P for the next time step
    P = 8.41*exp(0.972*(t+1)) #simple exponental growth for small populations
    #fit to values in review by Booth and Dhami
    f = cos((2*pi*(m-mmax))/12)*55000 #Seasonal oscillaiton added to population size
    #g =(1+(150*exp(-1.26*(a/12)))) #Weighting factor for oscillation based on age (should be based on size)
    g = 1388407*pop^(-1.1678667) #Weighting factor for oscillation based on population size (power function approximation)
    return(P + f/g)
  }
  
  if(!queen) { #TRUE if queen is dead.
    #negative population growth when no queen is present
    newpop <- pop - pop/(zsurv + age + 1) #if the queen is dead then age should be negative 
    return(newpop) #need something better here!!!
  }
}

agecol <- function(age, queen) {
  if(queen) {
    newage <- age + 1 #If the queen is alive, just advance the age
  } else {
    if(age >= 0) { #IF the queen is dead, start aging the zombie colony from -1 
      newage <- -1
    } else {
      newage <- age - 1
    }
  }
}

makeQueens <- function(col, qp) { #make a data frame for new queens
  #New data frame for the queens
  if(!any(col$pop >= 30000)) { #If not any colony has at least 30,000 workers then no queens are produced
    queens <- NULL
  } else {
    colq <- col[col$pop >= 30000 & col$queen != 0,] #subset down to just the colonies with at least 30,000 workers and a live queen
    
    wmass <- (0.086*colq$pop^0.178)*colq$pop #total worker mass
    #avg worker mass increases with population size => 0.086*colq$pop^0.178; from Tschinkel 1993 Fig 23.
    qmass <- wmass * (ifelse(colq$pop >= 50000, 0.35, 0.16))  #16% of mass devoted to queens in colonies of size class 1. 34-36% in larger size classes (Tschinkle 1993)
    qmass <- qmass*0.75
    qmassc <- qmass*colq$ptype #mass of claustral queens
    qmassp <- qmass*(1-colq$ptype) #mass of parasitic queens
    
    qnumc <- round((qmassc/7.2)*(1-cmort))  #number of claustral queens assuming average mass of 7.2mg and given the random mortality rate
    qnump <- round((qmassp/4.7)*(1-pmort))  #number of parasitic queens assuming average mass of 4.7mg and given the random mortality rate 
    
    queens <- data.frame(x = colq$x, y = colq$y,   #Starting coordinates (same as parent)
                         area = .1,    #initial area in sq meters
                         pop = 0,     #inital population size
                         queen = 1, #1 for queen alive, 0 for queen dead
                         #qcnt = 0, #queen production, always zero to start out
                         type = rep("C",nrow(colq)), #make all claustral at first
                         ptype = colq$ptype, #Which proportion is claustral 
                         age = 0, #age in MONTHS!!
                         ID = 0, #incremenal ID no. fill this in later
                         gen = colq$gen + 1, #generation (always 1 + parent generation)
                         birth = tick$t, #when colqony was founded
                         lin = colq$lin, #orignial lineage ID
                         fill = colq$fill, #same as parent colony
                         stringsAsFactors = F)
    
    #Determine baseline queen production per colony.
    #Make a vector qprod (queen production) for number of queens produced by each colony.
    #Queen production differs for each month (May be too simpe)
    pqns <- queens[rep(seq_len(nrow(queens)), qnump),]
    if(nrow(pqns) > 0) {pqns$type <- "P"}
    cqns <- queens[rep(seq_len(nrow(queens)), qnumc),]
    queens <- rbind(pqns, cqns)
  }
  return(queens)
}

queenMort <- function(col) { #queen survival function
  #combines two functions to derive queen survivial
  #1) y = 1-exp(-0.002*(col$pop+1250)) #asymptotic function for monthly survival based on colony size--Queens in smaller colonies are more likely to die (just made this up - needs better parameterization) 
  #2) y = (2/10^30)*col$age^16 + 1 #a very flat power function that induces queen death at about six years.
  qsp <- (1-exp(-0.002*(col$pop+1250))) * (-(2/10^30)*col$age^16 + 1)
  rnd <- runif(n=nrow(col))  #generate some random numbers
  col$queen <- as.numeric(qsp > rnd) * col$queen  #true/false conditional for survival probability > a random number; multiply by col$queen so that dead queens stay dead.
  return(col)
}

#------------MAIN PROGRAM-----------------------------

#------------ESTABLISH INITIAL COLONIES----------------

startcol <- 50 # how many colonies to startwith
startext <- (extent(c(0,xmax, 0, 18))) #extent of the starting zone
startsize <- 3 # approximate length/width of starting grid cell
startgrid <- raster(startext, nrows = round(startext@ymax/startsize), ncols = round(xmax/startsize))
startgrid <- rasterToPolygons(startgrid)
startgrid <- RoundPolygons(startgrid)
startgrid <- startgrid[sample(1:length(startgrid), startcol),]
gridarea <- startgrid@polygons[[1]]@area
maxpop <- gridarea/0.00064  #maximum population give area constraint

# Make initial lineages

lineages <- data.frame(x = rep(1,5),       #x coordinates
                       y = rep(1,5),       #y coordinates
                       area = 1,    #initial area in sq meters
                       pop = 1000, #initial population size
                       queen = 1, #1 for queen alive, 0 for queen dead
                       type = "C",
                       ptype = c(0.98, 0.95, 0.90, 0.75, 0.50), #What proportion is claustral 
                       age = 12, #initial age in MONTHS!!
                       ID = c(1:5), #incremenal ID no.
                       gen = 1, #generation (always 1 + parent generation)
                       birth = 0, #when colony was founded
                       lin = c(1:5), #orignial lineage ID
                       fill = (c("red", "orange", "yellow", "green", "blue")), #colors
                       stringsAsFactors=F) #avoid pesky factor designations

#Expand representation from each lineage
lineages <- lineages[rep(1:nrow(lineages),length.out=startcol),]
lineages <- lineages[sample(startcol),] #randomize order
lineages$ID <- c(1:50) #assign a sequence of  IDs
row.names(lineages) <- lineages$ID #match IDs with row names

#Assign lineage data to the grid cells
startgrid <- spChFIDs(startgrid, as.character(lineages$ID))
startgrid <- SpatialPolygonsDataFrame(startgrid, data = lineages)

#establish the actual x,y coordinates and iniital populations
coords <- as.data.frame(coordinates(startgrid))
startgrid$x <- coords$V1
startgrid$y <- coords$V2
startgrid$age <- sample(1:36, size=50, replace=T) #assign ages randomly from 1 to 36 months
#Calculate starting population size based on age - Using logistic equation from Tschinkle et al 1993.
startgrid$pop <- 165000/(1+83*exp(-1.26*(startgrid$age/12))) + (cos(6.28*startgrid$age/12)*55000)/(1+150*exp(-1.26*(startgrid$age/12)))

#match up initial areas
startgrid$area <- gArea(startgrid, byid = T)

#----------------Main loop--------------------------

cper <- c(1, 0.98, 0.95, 0.90, 0.75, 0.50) #array of claustral queen proportions

for(g in cper) {
  col <- startgrid #define colony spatial polygons data frame
  col$ptype <- g
  oldcol1 <- NULL
  tick <- list(t = 1, m = 1) #establish a list of loop parameters #t = time step (month), m = month of year
  j = 1
  while(j <= 300) {
    
    maxID <- max(col$ID) + 1
    print(Sys.time())
    cat("START step =", tick$t, " mo =", tick$m, " colonies =", nrow(col), "\n", sep = " ")
    
    #Caculate queen mortality 
    col <- queenMort(col)
    col$fill[col$queen == 0] <- "gray50" #recolor colonies with dead queens
    #Determine expected new population size
    col$pop <- mapply(FUN = pop6, pop = col$pop, m = tick$m, queen = col$queen, type = col$type, age = col$age)
    
    #REMOVE DEAD COLONIES (no queen and population below 100)
    cat("dead colonies:", col$ID[which(col$queen == 0 & col$pop < 100)], "\n", sep = " ")
    col <- col[c(which(col$queen != 0 | col$pop > 100)),] #keep only viable colonies 
    
    
    #save old status of simulation in case there is an error with expand
    
    if(tick$m == 1) {
      oldcol2 <- oldcol1
      oldcol1 <- list(ocol = col, otick = tick, oj = j)
    }
    
    #Determine colony area based on current area and projected growth
    cat("EXPAND\n")
    col <- tryCatch(expnd(col),
                    error = function(restart) {cat("major error with Expnd function. Backing up and starting over.");
                                               tick <<- oldcol2$otick;
                                               j <<- oldcol2$oj;
                                               oldcol2$ocol } )          
    cat("EXPAND DONE\n")
    col <- CleanTop(col) #Attempt to clean the topology
    col$ID <- as.numeric(col$ID) #Make sure ID is numeric
    maxID <- max(col$ID) + 1 #reestablish the max ID number
    col <- RoundPolygons(col)
    col <- CatchTopErr(col, "col")
    #trim data to arena and adjust area values in the data frame
    col <- crop(col, arena)
    col$area <- gArea(col, byid = T) #update area field in data frame
    areas <- sapply(slot(col, "polygons"), slot, "area") #get areas from spatial polygons datafram.
    #expand colony populations according to the area each colony has accrued.
    expPop <- sapply(areas, expWrkr) #get expected populations for given area
    col$pop <- mapply(FUN = min, col$pop, expPop) #final population is the minimum of the potential growth and the area-limited growth    
    #plot colonies
    try(dev.off)
    plot(col, add = F, col = col$fill, ylim = c(0,max(col$y)+pmd))
    
    #Calculate queen production
    if(tick$m == 1) {  #calculate queen production based on january populaiton sizes
      queens <- makeQueens(col, qp) #determine queen production
      pqueens <- cqueens <- ptot <- ctot <- data.frame() #set initial valuse to 0
      if(!is.null(queens)) {
        pqueens <- queens[queens$type == "P",] #subset of parasitic queens
        ptot <- nrow(pqueens) #total number of parasitic queens produced
        cqueens <- queens[queens$type == "C",] #subset of claustral queens
        ctot <- nrow(cqueens) #total number of claustral queens produced
      } 
    }
    
    #DISPERSE PARASITIC QUEENS
    if(nrow(pqueens) > 0) {
      samp <- sample(1:nrow(pqueens), size=ptot*pqd[tick$m]) #generate a random vector of row numbers which will constitute the parasitic queens that disperse during this time step
      if(length(samp) > 0) { #disperse parsitic queens, but only during indicated months and if there are zombie colonies present
        cat("Dispersing", length(samp), "parasitic queens.", "\n", sep = " ")
        parasitic <- pqueens[samp,] #subset of parasitic queens to disperse during this time step
        pqueens <- pqueens[-samp,] #subset of parasitic queens remaining after dispersal
        if(any(col$queen == 0)) { #only disperse parasitic queens if there are any zombie colonies
          parasitic <- disp2(parasitic, pmd, psd, xmax) #run the dispersal function to get new coordinates
          parasitic <- parasitic[sample(nrow(parasitic)),] #randomize row order
          coordinates(parasitic) <- c("x", "y")
          ############points(parasitic$x, parasitic$y, pch = 2, cex = 0.5, col = parasitic$fill) #plot the points
          zombies <- col[col$queen == 0,] #subset of colonies without living queens
          alive <- over(parasitic, zombies) #Dataframe index of parasitic queens that landed in zombie colonies
          if(any(!is.na(alive$ID))) { #test to see if any parasitic queens landed in zombie colonies 
            alive <- alive[!is.na(alive$ID),]
            alive <- alive[!duplicated(alive$ID),]
            for(k in 1:nrow(alive)){
              colrow <- which(col$ID==alive$ID[k]) #which row in col dataframe has its data replaced
              parrow <- which(row.names(parasitic)==row.names(alive)[k]) #which row from the parasitic dataframe provides the replacement data
              col@data[colrow, 5:13] <- parasitic@data[parrow, 3:11] #replace all data except x coordinate, y coordinate, area, and population
              col@data$ID[colrow] <- maxID #give a new ID number to the parasitized colony
              maxID = maxID + 1 #increase the max ID number for next use
            }
            newparID <- seq(from = maxID - 1, by = -1, length.out = k)
          }    
        }  
      } #end of parasitic queen dispersal routine
    } #end of if statement
    
    #DISPERSE CLAUSTRAL QUEENS
    if(nrow(cqueens) > 0)  {
      samp <- sample(1:nrow(cqueens), size=ctot*cqd[tick$m]) #generate a random vector of row numbers which will constitute the claustral queens that disperse during this time step
      if(length(samp) > 0) { #disperse claustral queens, but only when there are queens available
        cat("Dispersing", length(samp), "claustral queens.", "\n", sep = " ")
        claustral <- cqueens[samp,] #subset of claustral queens to disperse during this time step
        cqueens <- cqueens[-samp,] #subset of claustral queens remaining after dispersal
        claustral <- disp2(claustral, cmd, csd, xmax) #run the dispersal function to get new coordinates
        if(nrow(claustral) > 0) { #check to see that there is still at least one queen (queens that fall outside the allowed y axes are eliminated)
          claustral <- claustral[sample(nrow(claustral)),] #randomize row order
          claustral$x1 <- claustral$x #repeated coordinates columns for making spatial points class
          claustral$y1 <- claustral$y
          coordinates(claustral) <- c("x1", "y1")
          points(claustral$x, claustral$y, cex = 0.5, col = claustral$fill)
          claustral <- claustral[is.na(over(claustral, col)$x),] #only claustral queens that land in empty space (where x is NA) survive  
          if(nrow(claustral) > 0) {
            claustral$ID <- seq(from=maxID, by=1, length.out= nrow(claustral)) #establish unique ID numbers (necessary for spatial object designation)
            maxID <- max(claustral$ID) + 1 #increase max ID appropriately
            
            #Define initial areas (squares) for colonies
            sqrs <- makeSquares(claustral@data) 
            row.names(claustral@data) <- claustral$ID
            claustral <- SpatialPolygonsDataFrame(sqrs, claustral@data)
            #Colonies in overlapping squares compete and a winner is chosen randomly
            #first randomize order of spatial data
            claustral <- claustral[sample(nrow(claustral)),] #may not need to do this. It has been randomized once already
            if(sum(is.na(over(claustral, col)$x)) > 0) { #only proceed if there are any non-overlapping colonies
              claustral <- claustral[is.na(over(claustral, col)$x),] #remove new colonies that overlap with old ones
              overs <- gIntersects(claustral, spgeom2 = NULL, byid = T) #Find all intersecting polygons
              upperTriangle(overs, diag=T) <- F  #Eliminate the upper triangle of the matrix of the T/F overlap matrix (it's redundant)
              melted <- melt(overs) #melt the matrix to get row and column pairs for each matrix value
              overs <- melted[which(melted$value),] #Toss out all of the FALSE values (no overlap): The remaining lines have overlaps that need to be dealt with
              if(nrow(overs)>0) { #if there are some overlaps, get rid of one of the founding colonies
                claustral <- claustral[!(claustral$ID %in% overs$X1),]
              }
              col <- spRbind(col, claustral)
            } #end check
          } #end check on whether any claustral queens survived
        } #end check to see if there are any queens
      } #end of claustral queen dispersal.
    } 
    
    #print report
    cat("END step =", tick$t, "mo =", tick$m, "colonies =", nrow(col), "\n\n", sep = " ")
    
    #advance counters
    col$age <- mapply(FUN = agecol, age = col$age, queen = col$queen)
    tick$t <- tick$t + 1
    tick$m <- ifelse(tick$m==12, 1, tick$m+1)
    
    j = j + 1
  }
  sname <- paste("colonies",g, Sys.time())
  sname <- gsub(" ", "_",sname, fixed = T)
  save(col, file=sname)
}