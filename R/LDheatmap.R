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
             add.map=TRUE, geneMapLocation=0.15, 
             geneMapLabelX=0.5, geneMapLabelY=0.3, 
             SNP.name=NULL, color=heat.colors(20),
             newpage=TRUE, name="ldheatmap")
{

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
      LDmatrix[lower.tri(LDmatrix)] <- NA
    }
    else if(!missing(gdat))  
      stop(paste("No method for an object of class",class(gdat)))
    else
      stop("Need to supply LD matrix or genotypes")

  #--------------------------------------------------------------------#
  ## Draw the heat map 

  if(color[1]=="blueToRed") color = rainbow(20, start=.7, end=0, s=.7)[20:1]

  if(newpage)grid.newpage()
  mybreak <- 0:length(color)/length(color)
  colcut <- as.character(cut(1-LDmatrix,mybreak,labels=as.character(color)))
  if(is.numeric(color)) colcut <- as.integer(colcut)
  fontsize <- convertX(unit(1/30,"grobwidth", rectGrob()), "points")
  ImageRect<-makeImageRect(dim(LDmatrix)[1],dim(LDmatrix)[2],colcut, name="heatmap")
  title <- textGrob(title, 0.5, 1.05, gp=gpar(cex=1.2), name="title")
  heatMap <- gTree(children=gList(ImageRect, title), name="heatMap")
  #--------------------------------------------------------------------#
  ## Draw a diagonal line indicating the physical or genetic map positions of the SNPs
  nsnps <- ncol(LDmatrix)
  step <- 1/(nsnps-1)
  if(add.map) {
      ind <- match(SNP.name, names(LDmatrix[,1]))
      geneMap <- LDheatmap.Map.add (nsnps, genetic.distances=genetic.distances,
                     geneMapLocation=geneMapLocation,
                     geneMapLabelX=geneMapLabelX,
                     geneMapLabelY=geneMapLabelY,
                     distances=distances,
                     SNP.name=SNP.name, ind=ind)
  }  else geneMap <- NULL
  #--------------------------------------------------------------------#
  Key <- LDheatmap.Legend.add(fontsize, color)  # draw the Color Key
  heatmapVP <- viewport(width = unit(.8, "snpc"), height = unit(.8, "snpc"), 
                        gp=gpar(fontsize=fontsize), name="heatmapVP")
  LDheatmapGrob<-gTree(children=gList(heatMap, geneMap, Key), vp=heatmapVP, name=name)
  grid.draw(LDheatmapGrob)
  ldheatmap <- list(LDmatrix=LDmatrix, LDheatmapGrob=LDheatmapGrob, heatmapVP=heatmapVP,                       
                    genetic.distances=genetic.distances, distances=distances)
  class(ldheatmap) <- "LDheatmap"
  invisible(ldheatmap)
} # function LDheatmap ends


#_______________________HeatMap Coordinates_________________________##

# Adds a symbol to the i,jth cell of the heatmap. 
# The defaultis to add a symbol to the diagonal (j=i). 
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
    pushViewport(LDheatmap$heatmapVP)
    Symbols <- pointsGrob((i-0.5)*1/nSNP, (j-0.5)*1/nSNP, pch=pch, gp=gp, name="Symbols")
    grid.draw(Symbols)
    popViewport()
    invisible(list(x=(i-0.5)*1/nSNP,y=(j-0.5)*1/nSNP))
}

#_______________________Block Highlighting_________________________##
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

LDheatmap.highlight <- function(LDheatmap, i, j, fill="NA", col="black", lwd=1, lty=1){                       
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
  pushViewport(LDheatmap$heatmapVP)
  highlight <- polygonGrob(x=pgon$x, y=pgon$y, gp=gpar(col=col, fill=fill, lwd=lwd, lty=lty),                name="highlight")
  grid.draw(highlight)
  popViewport()
  invisible(pgon)
}

#_______________________Color Bar _________________________##
# Draw the Color Key
LDheatmap.Legend.add <- function(fontsize, color){
    ImageRect<- makeImageRect(2,length(color), col=c(rep(NA,length(color)),color[length(color):1]), "colorKey")
    #Adding the label 'Color key'
    title<-textGrob("Color Key", x=0.5, y=1.25, name="title")

    #Adding labels to the color key
    labels<-textGrob(paste(0.2*0:5), x=0.2*0:5,y=0.25, gp=gpar(cex=0.9), name="labels")

    #Drawing ticks at the bottom axis of the color key
    ticks<-segmentsGrob(x0=c(0:5)*0.2 , y0=rep(0.4,6), x1=c(0:5)*0.2 , y1=rep(0.5,6),name="ticks")

    #Drawing a box around the color key
    #box <- rectGrob(x=0, y=0.5, height=0.5, just=c("left", "bottom"), name="box")
    box <- linesGrob(x=c(0,0,1,1,0), y=c(0.5,1,1,0.5,0.5), name="box")

    keyVP <- viewport(x=1.1, y=-.10, height=.10, width=.5, just=c("right","bottom"), 
    gp=gpar(fontsize=fontsize))
    key <- gTree(children=gList(ImageRect, title, labels, ticks, box), name = "Key", vp=keyVP)
    key
}
#_______________________Genetic Map Editing _________________________##
# adds a genetic map to the heatmap plot along the diagonal
LDheatmap.Map.add <- function(nsnps, genetic.distances,
                       geneMapLocation=0.15,
                       geneMapLabelX=0.65, geneMapLabelY=0.3,
                       distances="physical",
                       SNP.name=NULL, ind=NULL){
  snp <- ((1:nsnps-1) + 0.5) / nsnps
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
              gp=gpar(cex=1.1), just="left", name="title")

  ## Labelling some SNPs 
  geneMap <- gTree(children=gList(diagonal, segments, title),name="geneMap")
  if (!is.null(SNP.name) && !is.null(ind)){
    symbols <- pointsGrob(snp[ind],snp[ind],pch="*",
                gp=gpar(cex=2, bg="blue", col="blue"), name="symbols")
    SNPnames <- textGrob(paste(" ", SNP.name), just="left", rot=-45,
              regionx[ind], regiony[ind], gp=gpar(cex=0.8, col="blue"), name="SNPnames")
    geneMap <- gTree(children=gList(diagonal, segments, title, symbols, SNPnames),name="geneMap")
  } # if end
 
 geneMap
}


myRainbow.colors <- function (n)
{
    if ((n <- as.integer(n[1])) > 0) {
        i <- n%/%5
        j <- n - 2*i
                k <- as.integer(j)%/%2
                l <- j-k
        c(rainbow(2*i, start = 0, end = 1/6),
                if (k > 0) hsv(h = seq(from = 23/60, to = 11/60, length=k)[k:1]),
        if (l > 0) hsv(h = seq(from = 41/60, to = 29/60, length = l)[l:1]))
    }
    else character(0)
}


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
