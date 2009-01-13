# ldheatmap - Plots measures of pairwise linkage disequilibria for SNPs
# Copyright (C) 2004  J.Shin, S. Blay, N. Lewin-Koh, J.Graham, B.McNeney

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

########################################################################

"LDheatmap"<-
   function (gdat, genetic.distances=NULL, 
             distances="physical", LDmeasure="r", title="Pairwise LD",
             add.map=TRUE, add.key=TRUE, geneMapLocation=0.15, 
             geneMapLabelX=0.5, geneMapLabelY=0.3, 
             SNP.name=NULL, color=grey.colors(20),
             newpage=TRUE, name="ldheatmap", vp.name=NULL,
             pop=FALSE)
{

  makeImageRect <- function(nrow, ncol, cols, name) {
    xx <- (1:ncol)/ncol   
    yy <- (1:nrow)/nrow
    right <- rep(xx, nrow)
    top <- rep(yy, each=ncol)
    rectGrob(x=right, y=top, 
           width=1/ncol, height=1/nrow, 
           just=c("right", "top"), 
           gp=gpar(col=NA, fill=cols),
           name=name)
  }

  #_______________________Color Bar _________________________##
  # Draw the Color Key
  LDheatmap.Legend.add <- function(color, vp=heatmapVP){
    ImageRect<- makeImageRect(2,length(color), col=c(rep(NA,length(color)),color[length(color):1]),
                              "colorKey")
    keyVP <- viewport(x=1.1, y=-.10, height=.10, width=.5, just=c("right","bottom"))
    #Adding the label 'Color key'
    title<-textGrob("Color Key", x=0.5, y=1.25, name="title", gp=gpar(cex=0.8))

    #Adding labels to the color key
    labels<-textGrob(paste(0.2*0:5), x=0.2*0:5,y=0.25, gp=gpar(cex=0.6), name="labels")

    #Drawing ticks at the bottom axis of the color key
    ticks<-segmentsGrob(x0=c(0:5)*0.2 , y0=rep(0.4,6), x1=c(0:5)*0.2 , y1=rep(0.5,6),name="ticks")

    #Drawing a box around the color key
    #box <- rectGrob(x=0, y=0.5, height=0.5, just=c("left", "bottom"), name="box")
    box <- linesGrob(x=c(0,0,1,1,0), y=c(0.5,1,1,0.5,0.5), name="box")

    key <- gTree(children=gList(ImageRect, title, labels, ticks, box), name = "Key",
                 vp=vpTree(vp,vpList(keyVP)))
    key
  }

  #_______________________Genetic Map Editing _________________________##
  # adds a genetic map to the heatmap plot along the diagonal
  LDheatmap.Map.add <- function(nsnps, add.map, genetic.distances, 
                       geneMapLocation=0.15,
                       geneMapLabelX=0.65, geneMapLabelY=0.3,
                       distances="physical", vp=heatmapVP, 
                       SNP.name=NULL, ind=0){
    snp <- ((1:nsnps-1) + 0.5) / nsnps
    if(add.map){
    min.dist <- min(genetic.distances) 
    max.dist <- max(genetic.distances)
    total.dist <- max.dist - min.dist
  
    # Drawing the diagonal line 
    seq.x <- c(0.5*geneMapLocation  + 1/(nsnps*2),
             1+0.5*geneMapLocation - 1/(nsnps*2))
    seq.y <- c(-0.5*geneMapLocation + 1/(nsnps*2),
             1-0.5*geneMapLocation - 1/(nsnps*2))
    diagonal<-linesGrob(seq.x, seq.y, gp=gpar(lty=1), name="diagonal")

    ## Adding line segments to the plot: (point1 <-> point2) 
    ## point1: relative position of a SNP on the scaled line
    ## point2: position of that SNP on the LD image  
    regionx <- seq.x[1] +
             ((genetic.distances-min.dist)/total.dist)*(seq.x[2]-seq.x[1])
    regiony <- seq.y[1] +
             ((genetic.distances-min.dist)/total.dist)*(seq.y[2]-seq.y[1]) 
    segments<-segmentsGrob(snp, snp, regionx, regiony, name="segments")

    ## Adding the text indicating Physical length of the region under study
    if (distances=="physical")
      mapLabel <- paste("Physical Length:", round((total.dist/1000),1),
                      "kb", sep="")
    else 
      mapLabel <- paste("Genetic Map Length:", round(total.dist,1),"cM",sep="")
      title <- textGrob(mapLabel, geneMapLabelX, geneMapLabelY,
              gp=gpar(cex=0.9), just="left", name="title")

    ## Labelling some SNPs 
    geneMap <- gTree(children=gList(diagonal, segments, title),name="geneMap", vp=vp)
    if (!is.null(SNP.name) && (any(ind!=0))){
      symbols <- pointsGrob(snp[ind],snp[ind],pch="*",
                gp=gpar(cex=1.25, bg="blue", col="blue"), name="symbols")
      SNPnames <- textGrob(paste(" ", SNP.name), just="left", rot=-45,
              regionx[ind], regiony[ind], gp=gpar(cex=0.6, col="blue"), name="SNPnames")
      geneMap <- gTree(children=gList(diagonal, segments, title, symbols, SNPnames),name="geneMap", vp=vp)
    }} # if(add.map) end

    else if (!add.map && !is.null(SNP.name) && (any(ind!=0))){
      geneMap <- textGrob(paste(" ", SNP.name), just="left", rot=-45,
                          snp[ind], snp[ind], gp=gpar(cex=0.6, col="blue"),
                          name="SNPnames", vp=vp)
    }

    else geneMap <- NULL

    geneMap
  }

  #--------------------------------------------------------------------#
  ## If genetic.distances is missing, calculate an equispaced default:
  if(is.null(genetic.distances)) {
     if (inherits(gdat,"data.frame"))
        genetic.distances=1000*(1:ncol(gdat))
     else if(inherits(gdat,"matrix"))
        genetic.distances=1000*(1:length(gdat[1,]))
     else    # gdat is of class LDheatmap
        genetic.distances = gdat$genetic.distances
  }
  #--------------------------------------------------------------------#
  ## Calculate or extract LDmatrix

    if(inherits(gdat,"data.frame")){
      for(i in 1:ncol(gdat)) {
        if(!genetics::is.genotype(gdat[,i]))
          stop("column ",i," is not a genotype object\n")
      }

      ## Exclude SNPs with less than 2 alleles:
      gvars <- unlist(sapply(gdat, function(x) genetics::nallele(x) == 2))
      genetic.distances <- genetic.distances[gvars]
      gdat <- gdat[gvars]

      ## Sort data in ascending order of SNPs map position:
      if(!is.vector(genetic.distances))
        {stop("Distance should be in the form of a vector")}
      o<-order(genetic.distances)
      genetic.distances<-genetic.distances[o]
      gdat<-gdat[,o]
      myLD <- genetics::LD(gdat)
      if(LDmeasure=="r")
        LDmatrix <- myLD[[LDmeasure]]^2   
      else if (LDmeasure=="D'")
        LDmatrix <- abs(myLD[[LDmeasure]])  
      else 
        stop("Invalid LD measurement, choose r or D'.")      
    }
    else if(inherits(gdat,"LDheatmap")){
      LDmatrix <- gdat$LDmatrix
      distances <- gdat$distances
    }
    else if(inherits(gdat,"matrix")){
      if(nrow(gdat) != ncol(gdat))
        stop("The matrix of linkage disequilibrium measurements must be a square matrix")
      LDmatrix <- gdat
      LDmatrix[lower.tri(LDmatrix, diag=TRUE)] <- NA
    }
    else if(!missing(gdat))  
      stop(paste("No method for an object of class",class(gdat)))
    else
      stop("Need to supply LD matrix or genotypes")

  #--------------------------------------------------------------------#
  ## Draw the heat map
  heatmapVP <- viewport(width = unit(.8, "snpc"), height = unit(.8, "snpc"),
                       name=vp.name)

  if(color[1]=="blueToRed") color = rainbow(20, start=4/6, end=0, s=.7)[20:1]
  if(newpage)grid.newpage()
  mybreak <- 0:length(color)/length(color)
  colcut <- as.character(cut(1-LDmatrix,mybreak,labels=as.character(color), 
include.lowest=TRUE))
  if(is.numeric(color)) colcut <- as.integer(colcut)
  ImageRect<-makeImageRect(dim(LDmatrix)[1],dim(LDmatrix)[2],colcut, name="heatmap")
  title <- textGrob(title, 0.5, 1.05, gp=gpar(cex=1.0), name="title")
  heatMap <- gTree(children=gList(ImageRect, title), name="heatMap", vp=heatmapVP)
  #--------------------------------------------------------------------#
  ## Draw a diagonal line indicating the physical or genetic map positions of the SNPs
  nsnps <- ncol(LDmatrix)
  step <- 1/(nsnps-1)
  ind <- match(SNP.name, names(LDmatrix[,1]), nomatch=0)
  geneMap <- LDheatmap.Map.add (nsnps, genetic.distances=genetic.distances,
                     geneMapLocation=geneMapLocation,add.map, 
                     geneMapLabelX=geneMapLabelX,
                     geneMapLabelY=geneMapLabelY,
                     distances=distances, vp=heatmapVP, 
                     SNP.name=SNP.name, ind=ind)
  #--------------------------------------------------------------------#
  if(add.key) Key <- LDheatmap.Legend.add(color, vp=heatmapVP)  # draw the Color Key
  else Key <- NULL 
  #heatmapVP <- viewport(width = unit(.8, "snpc"), height = unit(.8, "snpc"),
  #                      name="heatmapVP")
  LDheatmapGrob<-gTree(children=gList(heatMap, geneMap, Key),
                       childrenvp=heatmapVP, 
                       name=name, cl="ldheatmap")
  grid.draw(LDheatmapGrob)
  if(pop){
    downViewport(heatmapVP$name)
    popViewport()} #pop the heat map viewport

  ldheatmap <- list(LDmatrix=LDmatrix, LDheatmapGrob=LDheatmapGrob, heatmapVP=heatmapVP,                  
                    genetic.distances=genetic.distances, distances=distances)
  class(ldheatmap) <- "LDheatmap"
  invisible(ldheatmap)
} # function LDheatmap ends


preDrawDetails.ldheatmap <- function(x) {
  fontsize <- convertX(unit(1/20,"grobwidth", rectGrob()), "points")
  pushViewport(viewport(gp=gpar(fontsize=fontsize)))
}

postDrawDetails.ldheatmap <- function(x) {
  popViewport()
}

#_______________________HeatMap Coordinates_________________________##

# Adds a symbol to the i,jth cell of the heatmap. 
# The default is to add a symbol to the diagonal (j=i). 
# i and j can be vectors.
LDheatmap.marks <- function(LDheatmap, i, j=NULL, pch=20, gp=gpar(...), ...){
    nSNP <- dim(LDheatmap$LDmatrix)[1]
    if(is.null(j)) j<-i
    ind <- i>j
    if(any(ind)){
      ind <- as.numeric(ind)
      ind[ind>0] <- i[ind>0]
      i[ind>0] <- j[ind>0]
      j[ind>0] <- ind[ind>0]
    }
    heatmap.vp <- LDheatmap$heatmapVP$name
    #If heatmap.vp is on the grid display list, i.e., it is included in the 
    #returned value of current.vpTree(), a[1] <- 1 else a[1] <- NA
    a <- grep(paste("[", heatmap.vp, "]", sep=""), as.character(current.vpTree()), 
fixed=TRUE)
    if(!is.na(a[1]))   seekViewport(heatmap.vp)
    else               pushViewport(LDheatmap$heatmapVP)
    Symbols <- pointsGrob((i-0.5)*1/nSNP, (j-0.5)*1/nSNP, pch=pch, gp=gp, name="symbols")
    symbols <- gTree(children=gList(Symbols), name="Symbols", cl="symbols")
    grid.draw(symbols)
    if(!is.na(a[1]))  upViewport(0)  #back to the root viewport
    else              popViewport() 
    invisible(list(x=(i-0.5)*1/nSNP,y=(j-0.5)*1/nSNP))
}

preDrawDetails.symbols <- function(x) {
  fontsize <- convertX(unit(1/20,"grobwidth", rectGrob()), "points")
  pushViewport(viewport(gp=gpar(fontsize=fontsize)))
}

postDrawDetails.symbols <- function(x) {
  popViewport()
}


LDheatmap.highlight <- function(LDheatmap, i, j, fill="NA", col="black", lwd=1, lty=1){

  # Highlights the perimeter of selected cells in the heatmap as a block
  backbone <- function(i,j,nSNP){
     x <- c(i-1,i-1,j-1)/nSNP
     y <- c(i,j,j)/nSNP
     cbind(x,y)
  }

  jiggle <- function(i,j,nSNP){
    c1 <- j-i
    nvert <- (2*c1)-1
    x <-c(j-1,rep((j-2):(j-c1),each=2))
    y <- c(rep((j-1):(j-(c1-1)),each=2),j-c1)
    cbind(x,y)/nSNP 
  }
                       
  nSNP <- dim(LDheatmap$LDmatrix)[1]
  if(length(i)>1 | length (j) > 1) stop("i and j must be scalar indices")
  if((i<1 | i>nSNP) |(j<1 | j>nSNP) )
    stop(paste("index out of bounds, i and j must be in (1,",nSNP,")",sep=""))
  if(i==j) stop("i cannot be equal to j")
  if(i>j){
     h<-i
     i <- j
     j <- h
  }
  pgon <- data.frame(rbind(backbone(i,j,nSNP), jiggle(i,j,nSNP)))
  ## Square or almost square interior Blocks
  names(pgon) <- c("x","y")
  heatmap.vp <- LDheatmap$heatmapVP$name
  #If heatmap.vp is on the grid display list, i.e., it is included in the 
  #returned value of current.vpTree(), a[1]=1 else a[1]=NA:
  a <- grep(paste("[", heatmap.vp, "]", sep=""), as.character(current.vpTree()), fixed=TRUE)
  if(!is.na(a[1]))   seekViewport(heatmap.vp)
  else               pushViewport(LDheatmap$heatmapVP)
  highlight <- polygonGrob(x=pgon$x, y=pgon$y, 
     gp=gpar(col=col, fill=fill, lwd=lwd, lty=lty), name="highlight")
  grid.draw(highlight)
  if(!is.na(a[1]))  upViewport(0)  #back to the root viewport
  else              popViewport() 
  invisible(pgon)
}

