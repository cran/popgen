\name{treesim}
\alias{treesim}
\title{Standard neutral coalescent simulator}
\description{
Generates a sample from the standard neutral coalescent model.
}
\usage{
treesim(n, Sn, theta, rho.vec, mutations = FALSE, sample = FALSE)
}
\arguments{
  \item{n}{The number of alleles to generate.}
  \item{Sn}{The number of sites between which recombination can occur.}
  \item{theta}{The mutation rate (= 2Nu where u is the mutation rate per
  gene per generation)}
  \item{rho.vec}{Cumulative recombination rate over the Sn sites.}
  \item{mutations}{Logical flag indicating whether the relative
    positions of the mutations should be returned.}
  \item{sample}{Logical flag indicating whether the sample of alleles
    should be returned.}
}
\value{A list with components
  \item{nss}{The number of segregating sites in the sample.}
  \item{tree.time}{The total time on the tree (measured in units of 2N generations)}
  \item{total.time}{The total time on the tree (measured in units of 2N generations)}
  \item{mutations}{A vector containing the positions of the mutations in
    the range [0, 1].}
  \item{sample}{A (n x nss) matrix containing the sample of alleles.}
}

\references{Hudson, R, R. (1990) Gene genealogies and the coalescent
  process. In Oxford Surveys in Evolutionary Biology, Vol. 7, pages
  1--43, OUP, Oxford.}
\author{Jonathan Marchini}

\examples{

 ## No recombination

 n <- 20
 Sn <- 1000
 theta <- 5
 rho.vec <- rep(0, Sn)
 a <- treesim(n, Sn, theta, rho.vec, mutations = TRUE, sample = TRUE)
 mat <- LDmat(a$sample, typ = "haplotype", pos = a$mutations)

 ## 2 hotspots of recombination

 n <- 20
 Sn <- 1000
 theta <- 20
 rho.vec <- c(rep(0, 200), rep(2000, 500), rep(4000, 300))
 a <- treesim(n, Sn, theta, rho.vec, mutations = TRUE, sample = TRUE)
 mat <- LDmat(a$sample, typ = "haplotype", pos = a$mutations)


}
\keyword{utilities}
