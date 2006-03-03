\name{LDheatmap.marks}
\alias{LDheatmap.marks}
\title{ Plots a symbol in the centers of cells of the heat map image}
\description{  The function \code{LDheatmap.marks} is used to plot
  a symbol in the centers of cells representing the 
  pairwise linkage disequilibria of specified pairs of SNPs. }
\usage{
LDheatmap.marks(LDheatmap, i, j = NULL, pch = 20, gp=gpar(...), ...)
}
\arguments{
  \item{LDheatmap}{ An object of class \code{"LDheatmap"} 
     (returned by the function \code{\link{LDheatmap}}). }
  \item{i}{ A vector of indices of the first set of SNPs. }
  \item{j}{ A vector of indices of the second set of SNPs. }
  \item{pch}{ Either an integer value or a single character specifying
     the symbol to be plotted. See \code{\link[graphics]{points}} 
     for possible values and their corresponding symbols. }
  \item{gp}{ Graphical parameters; See \code{\link[grid]{gpar}}. }
  \item{...}{ Graphical parameter settings to be passed onto \code{\link[grid]{gpar}} 
	function.} 
}
\value{
  \item{ x }{ The vector of x coordinate(s) of the plotted symbol(s). }
  \item{ y }{ The vector of y coordinate(s) of the plotted symbol(s). }
}
\details{ 
   The lengths of the vectors \code{i} and \code{j} must be the same and
   greater than or equal to 1.
   If the lengths are greater than 1, the function plots the specified
   symbol in the centers of the (i\eqn{\mbox{\textasciicircum}}{^}k, 
   j\eqn{\mbox{\textasciicircum}}{^}k)-th cells (for k=1,...K; K =
   length of the vectors \code{i} and \code{j}), where 
   i\eqn{\mbox{\textasciicircum}}{^}k and 
   j\eqn{\mbox{\textasciicircum}}{^}k are
   the k-th elements of vectors \code{i} and \code{j}, respectively.
   For example, if \code{i=c(1,2)} and \code{j=c(3,5)}, \code{LDheatmap}
   plots a symbol in the centers of the cells representing pairwise
   linkage disequilibria between the first and third SNPs and between the
   second and fifth SNPs in the genome of interest.  Note that the order
   of the sets of indices does not matter; for example,
   \code{LDheatmap.marks(LDheatmap, i=c(1,2), j=c(3,5))} is equivalent
   to \code{LDheatmap.marks(LDheatmap, i=c(3,5), j=c(1,2))}.
}
\author{ Nicholas Lewin-Koh <nikko@hailmail.net> }
\examples{
data(LDheatmapData)
tt <- LDheatmap(HapMap.dat, genetic.distances=distance)
LDheatmap.marks(tt, 15, 3, cex=1.6, col="blue")
}
\keyword{ aplot }
