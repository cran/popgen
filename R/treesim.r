
treesim <- function(n, Sn, theta, rho.vec, mutations = FALSE, sample = FALSE) {

  es <- theta * sum(1 / (1:(n - 1)))
  vs <- es + theta^2 * sum(1 / ((1:(n - 1))^2))
  s.max <- floor(es + 50 * sqrt(vs))
 
  a <- .C("treesim",
          as.integer(n),
          as.integer(Sn),
          as.double(theta),
          as.double(rho.vec),
          num.sites = integer(1),
          tree.time = double(1),
          total.time = double(1),
          as.integer(mutations),
          mutations = double(s.max),
          as.integer(sample),
          sample = integer(n * s.max),
	PACKAGE = "popgen"
          )
  
  if(sample) {
    sample <- matrix(a$sample[1:(n * a$num.sites)], n, a$num.sites, byrow = TRUE)
  }
  else
    sample <- NULL

  if(mutations) {
    mutations <- a$mutations[1:a$num.sites]
  }
  else {
    mutations <- NULL
  }

  return(list(nss = a$num.sites, tree.time = a$tree.time, total.time = a$total.time, mutations = mutations, sample = sample))
}
  
