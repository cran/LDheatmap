"LDheatmap" <- function (gdat, map.distance=1000*(1:ncol(gdat)), 
                        distances="physical", LDmeas="r", 
                        title=NULL, add.map=TRUE, x.image.show = 0.2,
                        y.image.show = 0.2, line.position=0.2, 
                        x.length.position=0, y.length.position=0, 
                        SNP.name=NULL, color=heat.colors(20))
{
    for(i in 1:ncol(gdat)) {
        if(!genetics::is.genotype(gdat[,i])) stop("column ",i," is not a genotype object\n")
    }

   if(color=="blueToRed")
	color = rainbow(20, start=.7, end=0, s=.7)[20:1]

    # Sort data in ascending order of SNPs map position:
    if(is.vector(map.distance)){
        o<-order(map.distance)
        map.distance.names<-names(gdat)
        map.distance<-map.distance[o]
        names(map.distance)<-map.distance.names
        gdat<-gdat[,o]}
    else {stop("Distance should be in the form of a vector")}
    
    myLD <- genetics::LD(gdat)
    tem <- c(1:5)
    ind <- names(myLD) == LDmeas
    if(LDmeas=="r")
        LD.measurement <- myLD[[tem[ind]]]^2
    
    else if (LDmeas=="D'")
        LD.measurement <- abs(myLD[[tem[ind]]])
    
    else
        stop("Invalid LD measurement, choose r or D'.")
    
    nsnps <- ncol(gdat)
    step <- 1/(nsnps-1)
    min.dist <- min(map.distance) 
    max.dist <- max(map.distance)
    total.dist <- max.dist - min.dist
    
    seq.xstart <- (0 - step*0.5)
    seq.xend <- 1.0 + step*0.5 
    seq.ystart <- 0 - step*0.5
    seq.yend <- 1.0 + step*0.5

    my.image.show <- c(x.image.show , y.image.show) 
    imagex.0 <- seq.xstart
    imagey.0 <- seq.ystart - my.image.show[2]

    imagex.1 <- seq.xend + my.image.show[1]
    imagey.1 <- seq.yend
    
    mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
    oma <- c(2,1,1,2)
    def.par <- par(no.readonly = TRUE) # save default, for resetting
    m <- matrix(c(1,1,1,2), 2,2)
    h <- (3 + mar.orig[2]) * par("csi") * 2.54 
    layout(m, heights=c(1, lcm(0.5*h))) 
    if (is.null(title)) { title = "Pairwise LD" }
    
    mybreak <- 0:length(color)/length(color)
    image(1-LD.measurement, xlim=c(imagex.0, imagex.1),
          ylim=c(imagey.0,imagey.1), axes=FALSE, main=title, cex.main=1,
          col=color, breaks=mybreak)
      
    if(add.map) {
        SNP.id <- 1:nsnps
        snp <- (SNP.id-1)*step
        
        my.lp<- line.position
        seq.x <- c(0.5*(my.lp), 1+0.5*(my.lp))
        seq.y <- c(-0.5*(my.lp), 1-0.5*(my.lp))
        lines(seq.x, seq.y, lty=1)

        # Adding line segments to the plot: (point1 <-> point2) 
        # point1 = relative position of a SNP on the scaled line
        # point2 = position of that SNP on the LD image  
        regionx <- seq.x[1] + ((map.distance-min.dist)/total.dist)*(seq.x[2]-seq.x[1])
        regiony <- seq.y[1] + ((map.distance-min.dist)/total.dist)*(seq.y[2]-seq.y[1]) 

	# D1, D2 are the lengths of the short line segments
	# to the middle of the bottom of the imaginery square
	# and the the middle of the right side of the imaginery square
	# just below the color squares of the LD plot.
#        D1 <- (regionx-snp)^2 + (regiony-snp+0.5*step)^2
#        D2 <- (regionx-snp-0.5*step)^2 + (regiony-snp)^2

#        ind <- (D1 < D2) 

#        if (all(ind==1))
#            segments(snp[ind] , (snp[ind]-0.5*step), regionx[ind], regiony[ind])
#        else if (all(ind==0))
#            segments((snp[!ind]+0.5*step),snp[!ind], regionx[!ind], regiony[!ind])
#        else{
            segments((snp), (snp), regionx, regiony)
#        }
 
        # Adding the text indicating Physical length of the region under study
        my.length.position <- c(x.length.position, y.length.position)
        xpos <- 0.65*(imagex.1+imagex.0) + my.length.position[1]
        ypos <- 0.3*(imagey.1+imagey.0) + my.length.position[2]

	if (distances=="physical")
        	text(xpos,ypos,paste("Physical Length:", round((total.dist/1000),1),"kb",sep=""), adj=0, cex=0.9)
	else
        	text(xpos,ypos,paste("Genetic Map Length:", round(total.dist,1),"cM",sep=""), adj=0, cex=0.9)
    }
    # Labelling some SNPs 
    if (!is.null(SNP.name)){
        for(i in 1:length(SNP.name)){
            ind <- names(gdat) == SNP.name[i]
            snpx <- snp[ind]
            snpy <- snp[ind]
            points(snpx,snpy, pch="*", cex=2, bg="blue", col="blue")
            text(regionx[ind], regiony[ind], paste(" ", SNP.name[i]), 
                                   cex=0.8, adj=0, srt=-45, col="blue")
        } # for end
    } # if end

    # Drawing the Color Key:
    a <- matrix(1:length(color), ncol=1)
    par(mar=c(2,0,1,2))
    image(a, axes=FALSE, col=color[length(color):1])
    box(bty="o")
    mylable <- as.character(seq(from=0, to=1, length=11))
    axis(side=1, at=seq(from=0-1/11/2, to=1+1/11/2, length=11), labels=mylable)
    mtext("Color Key", side=3, line=0.5, cex=0.75)
    par(def.par)

} # function ends
