\name{popdiv.convert}
\alias{popdiv.convert}
\title{Convert a dataset for use in function popdif.}
\description{Convert a matrix (in which the multi-locus genotypes of each
  individual are stored in 2 rows) to a format used by the function popdif.
}
\usage{
popdiv.convert(x, miss = 0)
}
\arguments{
  \item{x}{(2 * N, L + 1) matrix where N is the number of individuals
    and L is the number of loci. Each individual's genotypes are
    contained in 2 rows. The first entry of each specifies the
    individuals population. Each row has (L + 1) entries where L is the number of loci.  }
  \item{miss}{The code for missing values.}
}
\value{List with components
  \item{NUMA}{Vector of the number of alleles at each loci}
  \item{N}{LxP matrix containing the number of genotypes in each
    population at each locus.}
  \item{X}{Matrix with L rows, one for each locus. Each row contains a
    list of allele counts from each population. For example, if there
    are 2 populations and the locus has 3 alleles then the first 3
    entries in the row specify the allele counts in population 1 and the
    next 3 entries specify the allele counts in population 2. The rest
    of the row is set to 0.}
  }
\author{Jonathan Marchini}


\seealso{\code{\link{popdiv}}}

\examples{

 ex <- read.table(paste(.libPaths(), "popgen", "example_data2", sep = .Platform$file.sep))
 X <- matrix(unlist(ex), nrow(ex), ncol(ex), byrow = FALSE)

 X1 <- popdiv.convert(X)

}
\keyword{utilities}
