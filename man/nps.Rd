\name{nps}
\alias{nps}
\title{Non-parametric clustering of SNP genotype data.}
\description{Uses classical hierarchical clustering methods and the Gap
  statistic to identify clusters of genotype data.
}
\usage{
nps(X, dmetric = "manhattan", method = "mcquitty", gap.n = 100)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{X}{Data array of dimension c(n, 2, L) where n is the number of
    people and L is the number of loci. Entry [i, j, k] contains the jth
    allele of the ith person at the kth locus. The function only handles
  SNP data and the alleles should be coded 0 and 1.}
  \item{dmetric}{The distance measure to be used. This must be one of
          `"euclidean"', `"maximum"', `"manhattan"', `"canberra"' or
          `"binary"'. Any unambiguous substring can be given. For
	  genotype data we recommend the use of the `"manhattan"' metric.}
  \item{method}{The agglomeration method to be used. This should be (an
          unambiguous abbreviation of) one of `"ward"', `"single"',
          `"complete"', `"average"', `"mcquitty"', `"median"' or
          `"centroid"'.}
  \item{gap.n}{The number of simulations used to compute the Gap statistic.}
}
\details{The nps function can be used to cluster genotype data without
  making any parametric assumptions about the properties of the clusters
  themselves. In order to do this we use the existing clustering
  functions available in R and implement the recently proposed Gap
  statistic [1] to estimate the number of clusters in a given dataset. This method compares the average within-cluster dissimilarity to its expected value under an appropriate null reference distribution for each partition in the hierarchy. In practice the expectation is approximated using samples from the null reference distribution. The null reference distribution used assumes that all the datapoints lie in just one cluster. 
}
\value{A list with components
  \item{num}{The number of clusters selected by using the Gap statistic.}
  \item{clust}{The output of the function hclust applied to the real dataset.}
  \item{gap}{A vector of Gap statistics for each value of $k$.}
  \item{sk}{The adjusted standard error of each Gap statistic.}
}
\references{R. Tibshirani and G. Walther and T. Hastie (2001) Estimating
  the number of clusters in a dataset via the gap statistic. JRSS (B)
}
\author{Jonathan Marchini}

\seealso{\code{\link{nps.plot}}, \code{\link{ps}}}

\examples{

 ## NB. The value of gap.n is set unrealistically low in this example so
 ##     that the examples run efficiently at compile time. We suggest
 ##     gap.n be set to 100 (the default) for reliable results.

 X <- simMD(100, 3, 100, p = NULL, c.vec1 = c(0.1, 0.2, 0.3), c.vec2 = 1, ac = 2, beta = 1) 

 res <- nps(X - 1, gap.n = 2)

 nps.plot(res, k.max = 2)


}
\keyword{utilities}
