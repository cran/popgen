\name{popdiv}
\alias{popdiv}
\title{Population diversity function}
\description{
  This function assesses between population variability using
  unlinked multi-locus multi-allelic genotype data from the populations of interest.
}
\usage{
popdiv(L, P, NUMA, N, X, burnin = 10000, iter = 10000, m = 300, csd = 0.01, outc = FALSE)
}
\arguments{
  \item{L}{Number of loci}
  \item{P}{Number of populations}
  \item{NUMA}{Vector of the number of alleles at each loci}
  \item{N}{LxP matrix containing the number of genotypes in each
    population at each locus.}
  \item{X}{Matrix with L rows, one for each locus. Each row contains a
    list of allele counts from each population. For example, if there
    are 2 populations and the locus has 3 alleles then the first 3
    entries in the row specify the allele counts in population 1 and the
    next 3 entries specify the allele counts in population 2. The rest
    of the row is set to 0.}
  \item{burnin}{Number of burn-in iterations in MCMC run}
  \item{iter}{Number of sampling iterations in MCMC run}
  \item{m}{Scale parameter of Dirichlet distribution used to update
    global allele frequencies (p).}
  \item{csd}{Standard deviation of Normal distribution used to update
    the population variance parameters (c).}
  \item{outc}{Flag that specifies whether the samples of c are returned (Logical)}
}
\details{The model fitted has a Multinomial-Dirichlet structure with a
  specific parameterisation suitable for the situation.
  
  \eqn{x_{il_{j}}} = number of type \eqn{j} alleles observed at locus \eqn{l} in population \eqn{i}.
  
  \eqn{n_{il_{j}}} = number of type \eqn{j} alleles observed at locus \eqn{l} in population \eqn{i}.
  
  \eqn{\alpha_{il_j}} = locus \eqn{l} type \eqn{j} allele frequency of the
  population i.
  
  \eqn{\pi_{l_j}} = locus \eqn{l} type \eqn{j} allele frequency of the `global' population.
  
  \eqn{c_{i}} = variance parameter of population i.

  \eqn{J_l} = number of alleles at locus \eqn{l}.
  
  \eqn{\{x_{il_{1}}, \ldots, x_{il_{J_l}}\} \sim \textrm{Multinom}(\alpha_{il_1},
    \ldots, \alpha_{il_{J_l}})}
  
  \eqn{\{\alpha_{il_{1}}, \ldots, \alpha_{il_{J_l}}\} \sim
    \textrm{Dirichlet}(\pi_{l_1}\frac{(1 - c_{i})}{c_{i}},
    \ldots, pi_{l_{J_l}}\frac{(1 - c_{i})}{c_{i}})}

  The model is specified in a Bayesian framework by placing uniform
  priors on the \eqn{c} and \eqn{\pi} parameters.

  A Metropolis-Hastings algorithm is used to sample from the posterior
  distribution of the parameters \eqn{c} and \eqn{\pi}.

  The proposal distribution for the \eqn{c} parameters is

  \eqn{c^{\textrm{new}} \sim \textrm{N}(c^{\textrm{old}},
    \textrm{csd}^2)}

  The proposal distribution for the \eqn{\pi} parameter is

  \eqn{\pi_l^{\textrm{new}} \sim \textrm{Dirichlet}(m\pi_{l1}, \ldots, m\pi_{lJ_l})}
    
  The conjugate nature of the model allows the \eqn{\alpha} parameters to be
  integrated out of the likelihood so that these parameters are not
  sampled.

  The \eqn{c} parameters are analogous to Fst estimates which are
  traditional measures of population diversity. See [1] for more details.
  
  

}
\value{List with components
  \item{C.sample}{An (iter x P) matrix containing the samples of the c parameters.}
  \item{muc}{The posterior mean of the c parameters.}
  \item{sdc}{The posterior standard deviation of the c parameters.}
  \item{mup}{The posterior mean of the \eqn{\pi}  parameters.}
  \item{sdp}{The posterior standard deviation of the \eqn{\pi}  parameters.}
}
\references{
  [1] Nicholson et al. (2002), Assessing population
  differentiation and isolation from single-nucleotide polymorphism
  data. JRSS(B), 64, 695--715

  [2] Marchini and Cardon (2002) Discussion of Nicholson et
  al. (2002). JRSS(B), 64, 740--741

  [3] Nichols and Balding (1995) A method for quantifying differentiation
  between populations at multi-allelic loci and its implications for
  investigating identity and paternity. Genetica, 96, 3--12.
  
}
\author{Jonathan Marchini}



\examples{

## Note : In the examples below the number of burn-in and sampling
## iterations is set to 100 so that these examples run in a reasonable
## time. In practice we recommend setting these parameters much higher so
## that convergence of the Markov chain to the stationary distribution is
## more likely. 

 ###############
 ## EXAMPLE 1 ##
 ###############

 ## SNP dataset from Nicholson et al. (2002). 66 loci, 4 populations 

 ex1 <- read.table(system.file("example_data1", package = "popgen")) 
 ex1 <- matrix(unlist(ex1), 66, 13, byrow = FALSE)
 NUMA <- ex1[,1]
 N <- ex1[, 2:5]
 X <- ex1[, 6:13]
 L <- 66
 P <- 4

 res <- popdiv(L, P, NUMA, N, X, burnin = 100, iter = 100, m = 10, csd = 0.01, outc = TRUE)

 plot(c(0, 100), c(0, max(res$C.sample)), type = "n")
 lines(1:100, res$C.sample[, 1])
 lines(1:100, res$C.sample[, 2], col = 2)
 lines(1:100, res$C.sample[, 3], col = 3)
 lines(1:100, res$C.sample[, 4], col = 4)

 ###############
 ## EXAMPLE 2 ##
 ###############

 ## Simulated dataset from the Multinomial-Dirichlet model.
 ## 100 loci, 3 populations.
 ## Variance parameters set to 0.01, 0.02, 0.04.
 ## Each loci has a variable number of alleles (between 3 and 8).
 ## Dataset stored in the following format :-
 ##  - each individuals genotypes are stored in 2 rows of the file.
 ##  - the first entry in each row specifies the individuals population.
 ##  - each row has (L + 1) entries where L is the number of loci. 
 ## To examine the datafile use edit(ex2)

 
 ex2 <- read.table(system.file("example_data2", package = "popgen")) 
 X <- matrix(unlist(ex2), nrow(ex2), ncol(ex2), byrow = FALSE)

 X1 <- popdiv.convert(X)

 res1 <- popdiv(L = nrow(X1$N), P = ncol(X1$N), NUMA = X1$NUMA, N = X1$N, X = X1$X, burnin = 100, iter = 100, m = 10, csd = 0.01, outc = TRUE)

 plot(c(0, 100), c(0, max(res1$C.sample)), type = "n")
 lines(1:100, res1$C.sample[, 1])
 lines(1:100, res1$C.sample[, 2], col = 2)
 lines(1:100, res1$C.sample[, 3], col = 3)

  
}
\keyword{utilities}
