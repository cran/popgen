\name{nps.plot}
\alias{nps.plot}
\title{Gap statistic plot for nps function}
\description{Plot the Gap statistic for different numbers of clusters
  using the output of the SNP genotype clustering function nps.
}
\usage{
nps.plot(res, k.max)
}
\arguments{
  \item{res}{Object produced by the function nps}
  \item{k.max}{Maximum number of clusters to plot the Gap statistic for.}
}
\details{Plot the Gap statistic for different numbers of clusters
  using the output of the SNP genotype clustering function nps. A red
  line is used to indicate the selected number of clusters.
}
\references{R. Tibshirani and G. Walther and T. Hastie (2001) Estimating
  the number of clusters in a dataset via the gap statistic. JRSS (B)
}
\author{Jonathan Marchini}

\seealso{\code{\link{nps}}}

\keyword{utilities}
