\name{ps}
\alias{ps}
\title{Parametric model-based inference of population structure.}
\description{Parametric model-based inference of population structure for unlinked
  SNP datasets.
}
\usage{
ps(X, K, burn.in = 10^4, num.sample = 10^4, thin = 1, alpha.start = 1.0, alpha.sd = 0.05, alpha.step = 10, alpha.max = 20, z.init = NULL, sample.type = rep(1, 5), na.lab, c.init = NULL, pi.init = NULL, c.prior.mu = 0.1, c.prior.sd = 0.1, c.sd = 0.01, fix.alpha = FALSE, fix.pi  = FALSE, popdif.flag = TRUE, beta = 1, fix.z = FALSE)
}
\arguments{
  \item{X}{Data array of dimension c(n, 2, L) where n is the number of
    people and L is the number of loci. Entry [i, j, k] contains the jth
    allele of the ith person at the kth locus.}
  \item{K}{The number of populations.}
  \item{burn.in}{The length of the burn-in.}
  \item{num.sample}{The number of sampling iterations after the burn-in.}
  \item{thin}{The thinning frequency (default = 1)}
  \item{alpha.start}{The initial value of the admixture parameter alpha
    (default = 1)}
  \item{alpha.sd}{The standard deviation of the proposal distribution
    used to update alpha (default = 0.05)}
  \item{alpha.step}{The frequency at which alpha is updated (default = 10)} 
  \item{alpha.max}{The upper limit for alpha (default = 20)}
  \item{z.init}{Array with the same dimension as X which contains the
    initial values of the population labels for each allele at each
    locus in each person (the default is NULL which implies that Z is
    set randomly)}
  \item{sample.type}{Vector of 0/1 flags that specify which parameter samples
    are recorded. The order of the flags is (P, Q, Z, pi, c). The
    default is to sample all the parameters i.e. sample.type = rep(1, 5).}
  \item{na.lab}{Label for missing data.}
  \item{c.init}{Vector of length K containing the initial values of the
    c parameters. If NULL then set randomly.}
  \item{pi.init}{Vector of length L containing the initial values of the
    pi parameters. If NULL then set randomly.}
  \item{c.prior.mu}{Mean of the prior for c (default = 0.01)}
  \item{c.prior.sd}{Standard deviation of the prior for pi (default = 0.05)}
  \item{c.sd}{The standard deviation of the proposal distribution
    used to update c (default = 0.01)}
  \item{fix.alpha}{Flag that specifies that alpha be fixed at its
    initial value.}
  \item{fix.pi}{Flag that specifies that pi be fixed at its
    initial value.}
  \item{popdif.flag}{Flag to turn on the model of population
    differentiation used in the function popdif.}
  \item{beta}{Parameter that scales the dirichlet proposal for pi
    (default = 1)}
  \item{fix.z}{Flag that specifies that Z be fixed at its initial value
    i.e. the structure is fixed and the other parameters infered.}
}
\details{The ps function is an implementation of a model-based
  clustering algorithm for genotype data developed in [1]. The implementation here includes an extension which uses a more realistic prior structure for each subpopulations allele
  frequencies (as developed in [2], [3], and [4]). At present the
  function should only be used with unlinked SNP data.

  The model is specified in a Bayesian framework and samples from the
  posterior distribution are obtained from an MCMC algorithm. Details of the algorithm can be found in the references below.

  The ps function clusters a given dataset for a given number of
  populations K and it is often desirable to obtain inference on the
  number K. In [1] a novel approximation to the evidence was used to
  approximate the posterior distribution on K by combining the results
  of several runs. This can be achieved by combining the results of several
  runs using different values of K. The details are given in [1] and an
  example is given below in the examples section.
  
}
\value{List with components
  \item{Alpha}{Sample of the alpha parameter.}
  \item{Z}{Mean of the (n, 2, L) array Z i.e. the mean population labels for each allele.}
  \item{P}{Mean of the (K, 2, L) array P i.e. the mean allele frequencies at each
    locus in each population.}
  \item{Q}{Mean of the (n, K) matrix Q i.e. the mean proportion of each
    persons genome from each population.}
  \item{mu}{Mean value of log-likelihood}
  \item{v}{Variance of log-likelihood}
  \item{log.PX.K}{Estimated log-probability of the data given K}
  \item{PX}{Sample of the log-likelihood}
  \item{pi}{Sample of the pi parameters.}
  \item{c}{Sample of the c parameters.}
}
\references{
  [1] Pritchard, Stephens and Donnelly (2000) Inference of population
  structure using multilocus genotype data. Genetics 155:945-959
  
  [2] Nichols and Balding (1995) A method for quantifying differentiation
  between populations at multi-allelic loci and its implications for
  investigating identity and paternity. Genetica, 96, 3--12.
  Marchini and Cardon (2002) Discussion of Nicholson et
  al. (2002). JRSS(B), 64, 740--741

  [3] Nicholson et al. (2002), Assessing population
  differentiation and isolation from single-nucleotide polymorphism
  data. JRSS(B), 64, 695--715

  [4] Marchini and Cardon (2002) Discussion of Nicholson et
  al. (2002). JRSS(B), 64, 740--741


}
\author{Jonathan Marchini}
\examples{
  X <- simMD(60, 3, 100, p = NULL, c.vec1 = c(0.1, 0.2, 0.3), c.vec2 =
1, ac = 2, beta = 1)

  ## infer the population structure
  ## NB in this example the burn-in and no. of iterations have been set
  ##    way too low so that the examples run in a reasonable time at time
  ##    of installation. To get reliable answers set these parameters
  ##    much higher and check convergence of the chain using multiple runs.

  res2 <- ps(X,
           2,
           num.sample = 100,
           burn.in = 100,
           na.lab = (-1),
           c.init = rep(0.1, 2),
           c.prior.mu = 0.01,
           c.prior.sd = 0.05,
           c.sd = 0.05,
           alpha.start = 1.0,
           fix.alpha = FALSE,
           fix.pi  = FALSE,
           popdif.flag = TRUE,
           beta = 1,
           z.init = NULL,
           fix.z = FALSE)

  res3 <- ps(X,
           3,
           num.sample = 100,
           burn.in = 100,
           na.lab = (-1),
           c.init = rep(0.1, 3),
           c.prior.mu = 0.01,
           c.prior.sd = 0.05,
           c.sd = 0.05,
           alpha.start = 1.0,
           fix.alpha = FALSE,
           fix.pi  = FALSE,
           popdif.flag = TRUE,
           beta = 1,
           z.init = NULL,
           fix.z = FALSE)

  res4 <- ps(X,
           4,
           num.sample = 100,
           burn.in = 100,
           na.lab = (-1),
           c.init = rep(0.1, 4),
           c.prior.mu = 0.01,
           c.prior.sd = 0.05,
           c.sd = 0.05,
           alpha.start = 1.0,
           fix.alpha = FALSE,
           fix.pi  = FALSE,
           popdif.flag = TRUE,
           beta = 1,
           z.init = NULL,
           fix.z = FALSE)

  k.dist <- c(res2$log.PX.K, res3$log.PX.K, res4$log.PX.K)
  k.dist <- k.dist - max(k.dist)
  k.dist <- exp(k.dist) / sum(exp(k.dist))
  k.dist

  ## print the mean values of the c parameters
  mean(res3$c[, 1])
  mean(res3$c[, 2])
  mean(res3$c[, 3])

  ## print out the results
  par(mfcol = c(3, 3))
  hist(res3$c[, 1], n = 100, xlim = c(0, 0.4))
  hist(res3$c[, 2], n = 100, xlim = c(0, 0.4))
  hist(res3$c[, 3], n = 100, xlim = c(0, 0.4))

  plot(res3$c[, 1], typ = "l")
  plot(res3$c[, 2], typ = "l")
  plot(res3$c[, 3], typ = "l")

  image(res3$Q)

}

\keyword{utilities}

