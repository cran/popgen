\name{LDmat}
\alias{LDmat}
\title{Calculate LD measures for haplotype data.}
\description{Linkage disequilibrium measures (D' and
  $r^2$) for a set of genotypes or haplotypes.
}
\usage{
LDmat(mat, typ = c("genotype", "haplotype"), pos = NULL, plotmat = TRUE)
}
\arguments{
  \item{mat}{An nxp matrix of genotypes or haplotypes. In the case of
    genotypes the matrix should contain a row for each individual. The
    genotypes should be coded 0, 1 and 2 and use -1 for missing. In the case of
    haplotypes the matrix should contain a row for each haplotype. The
    alleles at each site should be coded 0 and 1 and use -1 for missing. }
  \item{typ}{Flag determining whether the data matrix contains genotypes
  or haplotypes.}
  \item{pos}{Vector with the positions of the markers in the range [0, 1].}
  \item{plotmat}{Flag specifying whether to plot the LD measures.}
}
\value{Returns a pxp matrix whose upper triangle contains the values of
  D' between pairs of markers and the lower triangle contains the values
  of $r^2$ between the pairs of markers. If plotmat == TRUE then an image is
  plotted in which the upper left segment is D' and the
  bottom right half is $r^2$.
}
\references{B. S. Weir (1996) Genetic Data Analysis II. Sinauer}
\author{Jonathan Marchini}

\keyword{utilities}
