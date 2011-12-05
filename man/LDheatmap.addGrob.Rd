\name{LDheatmap.addGrob}
\alias{LDheatmap.addGrob}
\title{
Add a graphical object to an LDheatmap plot
}
\description{
Add a graphical object to an LDheatmap plot such that the x-axis corresponds
to the physical map on the heatmap.
}
\usage{
LDheatmap.addGrob(LDheatmap, grob, height = 0.2)
}
\arguments{
  \item{LDheatmap}{
An object of class LDheatmap.
}
  \item{grob}{
A graphical object of class\code{\link[grid]{grob}}.
}
  \item{height}{
The height of the viewport in which the grob will be placed.
}
}
\value{
An object of class LDheatmap given as an argument, with the \code{grob}
\code{LDheatmapGrob} modified to inclue the new child grob.
}
\author{
Sigal Blay <sblay@sfu.ca>
}
\seealso{
\code{\link{LDheatmap}}
}
\examples{
% Add an empty rectangle frame
data(CEUData)
smallmap <- LDheatmap(CEUSNP, genetic.distances = CEUDist,
  color = grey.colors(20))
ll<-LDheatmap.addGrob(smallmap, rectGrob())
}
\keyword{ hplot }
