\name{LDheatmapData}
\alias{LDheatmapData}
\alias{HapMap.dat}
\alias{distance}
\docType{data}
\title{ Example data set for LDheatmap }
\description{
  HapMap.dat:  genotypes on 15 SNPs for 60 people

  distance:    Physical map positions of the 15 SNPs in HapMap.dat
}
\usage{data(LDheatmapData)}
\format{
  HapMap.dat: A dataframe of SNP genotypes.
  Each row represents an individual.
  Each column represents a SNP.

  distance: a vector of integers, representing SNP physical map
  locations on the chromosome.
}
\details{ Data on SNPs with minor allele frequency greater
  than 5\% from a 9kb region of chromosome 7 (base positions 126273659
  through 126282556 from release 7 of the International HapMap Project). 
  Genotypes from 30 parent-offspring trios (both
  parents, one offspring) were obtained.
  The 30 trios are taken from the so-called CEPH families, a set of
  multi-generational Caucasian families from Utah. From this set of 90
  people, 60 parents were extracted.
}
\source{ International HapMap Project \url{www.hapmap.org}
}
\examples{ data(LDheatmapData) }
\keyword{datasets}


