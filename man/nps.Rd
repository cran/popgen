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

Specifically, suppose we have clustered the data into $k$ clusters $C_1, \ldots, C_k$ with $n_r=\arrowvert C_r\arrowvert$. Then let
\begin{equation}
   D_r = \sum_{i, j \in C_r} d_{ij}
\end{equation}
 be the total within-cluster dissimilarity of a given cluster and define the within-cluster dissimilarity of the partition as 
 \begin{equation}
   W_k = \sum_{r=1}^k \frac{1}{2n_r} D_r
 \end{equation}
The Gap statistic is then defined to be
  \begin{equation}
    \textrm{Gap}_n (k) = \bf{E}_n^* \big( \log(W_k) \big) - \log (W_k)
  \end{equation}
where $\bf{E}_n^*$ denotes expectation under a sample of size $n$ from the reference distribution. The estimate of $\hat{k}$ will be the value maximising $\textrm{Gap}_n (k)$ after the sampling distribution is taken into account.

\subsection{Computation of the Gap statistic} \label{sec:Gap}
\begin{enumerate}
\item Apply a hierarchical clustering method to the dataset and for each partition (level of the hierarchy) calculate $W_k, \;k=1,\ldots, K$.
\item Generate $R$ datasets from a null model, cluster each one and calculate $W_{kr}^*, \;r=1,\ldots,R, \;k=1,\ldots,K$. Let
  \begin{eqnarray*}
    m_k & = & (1/R) \sum_{r=1}^R \log(W_{kr}^*\\
    sd_k & = & \sqrt{(1/R) \sum_{r=1}^R (\log(W_{kr}^* - m_k)^2}
  \end{eqnarray*}
\item The (estimated)Gap statistic is then defined to 
  \begin{equation}
    \textrm{Gap} (k)= m_k - \log(W_k)
  \end{equation}
and define 
  \begin{equation}
    s_k = sd_k \sqrt{1 + 1/R}
  \end{equation}
The term $\sqrt{1 + 1/R}$ is introduced to take account of the simulation error in $\bf{E}_n^* \big( \log(W_k) \big)$. Estimate the number of clusters to be
  \begin{equation}
    \hat{k} = \textrm{smallest } k \textrm{ such that Gap}(k) \geq \textrm{ Gap}(k+1) - s_{k+1}
  \end{equation}
\end{enumerate}
A nice property of this method is that the statistic is defined for $k =
1$ clusters whereas some previosly proposed statistics are not. Compared
to other methods of determining the number of clusters Tibshirani et al.
(2001) find that this method yields the highest performance using simulated datasets. 

\subsection{Applying the Gap statistic to SNP genotype data}

For a given clustering method (specified through the arguments dmetric
and method) we compute the Gap statistic by simulating datasets with one
cluster using the sample genotype frequencies at each locus.

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
