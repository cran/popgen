\name{simMD}
\alias{simMD}
\title{Simulate multi-population unlinked genotype data from a
  Multinomial-Dirichlet model.}
\description{Simulate multi-population unlinked genotype data from a
  Multinomial-Dirichlet model.
}
\usage{
simMD(N, P, L, p = NULL, c.vec1, c.vec2 = 1, ac = 2, beta = 1)
}
\arguments{
  \item{N}{The number of people per population.}
  \item{P}{The number of populations.}
  \item{L}{The number of unlinked loci.}
  \item{p}{A (ac x L) with column l containing the ac global allele
    frequencies for locus l. If NULL the frequencies are generated at
    random from a Dirichlet distribution with parameter beta (default is
    beta = 1 which is a uniform distribution.}
  \item{c.vec1}{A vector of length P which contains the level 1 variance
    parameter for each subpopulation.}
  \item{c.vec2}{An (optional) vector which contains the level 2 variance
    parameters for each subpopulation of the level 1 subpopulation.}
  \item{ac}{The number of alleles at each locus.}
  \item{beta}{The parameter of Dirichlet distribution used to simulate
    the global allele frequencies.}
}
\details{The data is simulated from a Multinomial-Dirichlet structure with a
  specific parameterisation suitable for the situation.
  
  \eqn{x_{ipal}} = Allele of the \eqn{a}th chromosome from the \eqn{i}th person in the \eqn{p}th population at locus \eqn{l}.
   
  \eqn{\alpha_{pl_j}} = locus \eqn{l} type \eqn{j} allele frequency of the
  population \eqn{p}.
  
  \eqn{\pi_{l_j}} = locus \eqn{l} type \eqn{j} allele frequency of the `global' population.
  
  \eqn{c_{p}} = variance parameter of population \eqn{p}.

  \eqn{ac} = number of alleles at each locus.
  
  \eqn{x_{ipal} \sim \textrm{Multinom}(\alpha_{pl_1},
    \ldots, \alpha_{pl_{ac}})}
  
  \eqn{\{\alpha_{pl_{1}}, \ldots, \alpha_{pl_{ac}}\} \sim
    \textrm{Dirichlet}(\pi_{l_1}\frac{(1 - c_{p})}{c_{p}},
    \ldots, pi_{l_{ac}}\frac{(1 - c_{p})}{c_{p}})}

  If the global allele frequencies are not specified then 
  \eqn{\{\pi_{l_1}, \ldots, \pi_{l_{ac}}\}  \sim
    \textrm{Dirichlet}(\beta, \ldots, \beta)}
  
}
\value{If the c.vec2 parameter is left as its default value then the
  data is simulated from the standard multinomial-direchlet model (see above)
  and the result is stored in an array with dimensions (N * P, 2, L) containing the data. Entry
  ((p - 1) * P + i, a, l) contains the simulated allele from the ith
  person in the pth population of the ath chromosome at the lth locus.

  If the c.vec2 vector is specified then each of the P populations has
  length(c.vec2) subpopulations with the c.vec2 parameters containing
  the variance parameters of the subpopulations. 
  
}
\references{Nicholson et al. (2002), Assessing population
  differentiation and isolation from single-nucleotide polymorphism
  data. JRSS(B), 64, 695--715

  Marchini and Cardon (2002) Discussion of Nicholson et
  al. (2002). JRSS(B), 64, 740--741

  Nichols and Balding (1995) A method for quantifying differentiation
  between populations at multi-allelic loci and its implications for
  investigating identity and paternity. Genetica, 96, 3--12.
  }
\author{Jonathan Marchini}

\examples{

  X <- simMD(60, 3, 100, p = NULL, c.vec1 = c(0.1, 0.2, 0.3), c.vec2 =
1, ac = 2, beta = 1) 

}
\keyword{utilities}
