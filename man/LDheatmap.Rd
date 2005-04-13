\name{LDheatmap}
\alias{LDheatmap}
\title{ This function produces a pairwise LD plot }
\description{ \dQuote{LDheatmap} is used to produce a graphical display,
  as a heat map, of linkage disequilibrium measures for SNPs. The heatmap is
     a false color image in the upper-left diagonal of a square plot.
     Optionally, a line parallel to the diagonal of the image indicating
     the physical or genetic map positions of the SNPs may be added, along
     with text reporting the total length of the genomic region considered.

}
\usage{
LDheatmap(gdat, map.distance = 1000*(1:ncol(gdat)), distances="physical",
LDmeas = "r", title = NULL, add.map = TRUE, x.image.show = 0.2, 
y.image.show = 0.2, line.position = 0.2, x.length.position = 0, 
y.length.position = 0, SNP.name = NULL)
}
\arguments{
  \item{gdat}{ SNP data: a data frame of genotype objects. }
  \item{map.distance}{ Vector of physical or genetic map distances, if known.
          Physical distances should be in bases, genetic distances
          should be in centiMorgans (cM). Default value is equispaced
          markers, 1000 bases apart. }
  \item{distances}{ A character string indicating whether distances are 
 	  "physical" or "genetic" distances. Default is "physical". If
          distances="physical", the text describing the total length of
          the region will be "Physical Length:XXkb" where XX is the length
          of the region in kilobases. If distances="genetic", the text
          will be "Genetic Map Length:YYcM" where YY is the length of the
          region in centiMorgans.}
  \item{LDmeas}{ A character value of LD measurement - either
    allelic correlation r\eqn{\mbox{\textasciicircum}}{^}2 or Lewontin's |D'|;
    default = \dQuote{r} for r\eqn{\mbox{textasciicircum}}{^}2; type
    \dQuote{D'} for |D'|. }
  \item{title}{ Character string for the title of the
    plot. Default = NULL will print out \dQuote{Pairwise LD}. }
  \item{add.map}{ If TRUE, a diagonal line indicating the physical
    or genetic map positions of the SNPs will be added to the plot. }
  \item{x.image.show}{ A numeric value which will stretch the image
    window horizontally by the given amount; default = 0.2
    will expand the plot area horizontally by 0.2 units to the right. }
  \item{y.image.show}{ A numeric value which will stretch the image
    window vertically by the given amount; default = 0.2 will
    expand the plot area vertically by 0.2 units downwards. }
  \item{line.position}{ Modifies the position of the line
    parallel to the diagonal of the matrix; the larger the value, the
    farther it lies from the matrix diagonal.}
  \item{x.length.position}{ Modifies x-coordinate of the text
    labelling the total length of the genomic region being considered; e.g.,
    x.length.position=0.1 will move the text to the right by
    0.1 units from the original position in the upper right
    quadrant of the plot. }
  \item{y.length.position}{ Modifies y-coordinate of the text
    labelling the total length of the genomic region being considered; e.g.,
    y.length.position=0.1 will move the text upwards by 0.1 units
    from the original position in the upper right quadrant of the plot. }
  \item{SNP.name}{ A vector of character string(s) of SNP name(s) to
    be labelled in the plot. }
}
\details{ The distance vector of SNP positions, if provided, is
  assumed to be in the same order as that of the SNPs in the data. }
\section{Warnings}{ It is not possible to draw multiple LDheatmap plots
  on a single graph area. }
\author{ Ji-hyung Shin <jshinb@sfu.ca>, Sigal Blay <sblay@sfu.ca>, Brad McNeney <mcneney@stat.sfu.ca>, Jinko Graham <jgraham@cs.sfu.ca> }
\seealso{ \code{\link[genetics]{LD}}, \code{\link[genetics]{genotype}}, \code{\link[graphics]{image}} }
\examples{ demo(LDheatmap) }
\keyword{ hplot }
