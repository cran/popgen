
simMD <- function(N, P, L, p = NULL, c.vec1, c.vec2 = 1, ac = 2, beta = 1) {
  
  rdirichlet <- function(n, a)
    ## pick n random deviates from the Dirichlet function with shape
    ## parameters a
    {
      l <- length(a);
      x <- matrix(rgamma(l * n, a), ncol = l, byrow = TRUE);
      sm <- x %*% rep(1, l);
      x / as.vector(sm);
    }
  
  rmultinom <- function(p) length(p) - sum(runif(1) <  cumsum(p)) + 1 
  
  if(is.null(p)) {
    p <- matrix(0, ac, L)
    for(l in 1:L)
      p[, l] <- rdirichlet(1, rep(beta, ac))
  }
  
  X <- array(0, dim = c(N * P * length(c.vec2), 2, L))
  
  for(nn in 1:P) {
    for(l in 1:L) {
      temp1 <- rdirichlet(1, p[, l] * (1 - c.vec1[nn]) / c.vec1[nn])

      for(j in 1:length(c.vec2)) {
        if(length(c.vec2 == 1)) {
          temp2 = temp1
        } else {
          temp2 = rdirichlet(1, temp1 * (1 - c.vec2[j]) / c.vec2[j])
        }
        
        for(i in 1:N)
          for(a in 1:2) {
            X[i + (nn - 1) * N * length(c.vec2) + (j - 1) * N, a, l] <- rmultinom(temp2)
          } 
      }
    }
  }
  
  return(X)

}


ps <- function(X,
               K,
               burn.in = 10^4,
               num.sample = 10^4,
               thin = 1,
               alpha.start = 1.0,
               alpha.sd = 0.05,
               alpha.step = 10,
               alpha.max = 20,
               z.init = NULL,
               sample.type = rep(1, 5),
               na.lab,
               c.init = NULL,
               pi.init = NULL,
               c.prior.mu = 0.1,
               c.prior.sd = 0.1,
               c.sd = 0.01,
               fix.alpha = FALSE,
               fix.pi  = FALSE,
               popdif.flag = TRUE,
               beta = 1,
               fix.z = FALSE){
    
    ## set up constants
    n <- dim(X)[1]
    L <- dim(X)[3]
    print.step <- round(0.1 * (num.sample * thin), 0)
    
    ## set missing values equal to (-99)
    adj <- 1 - min(X[X != na.lab])
    X <- X + adj
    X[X == (na.lab + adj)] <- (999)

    ## set number of classes vector and re-label X
    J <- apply(X, 3, FUN = function(x){ length(unique(x))})
    Jmiss <- apply(X, 3, FUN = function(x){ sum(x == 999) > 0})
    J <- J - Jmiss
    J <- pmax(J, 2) ## minimum number of alleles is 2
    maxJ <- max(J)
    
    ## re-label X
    labs <- matrix(0, maxJ, L)
    for(l in 1:L) {
        temp <- X[, , l]
        temp <- as.vector(temp[temp != 999])
        tmp <- sort(unique(temp))
        labs[1:length(tmp), l] <- tmp
    }
    for(l in 1:L) {
        tmp <- X[, , l]
        for(j in 1:J[l]) {
          
            tmp[tmp == labs[j, l]] <- j
        }
        X[, , l] <- tmp
    }

    ## initialize Z
    init.z <- 1
    if(is.null(z.init)) {
        init.z <- 0
        z.init <- X
    }

    ## initialize pi and c
    pi.init.set <- function(X, maxJ, J) {
        
        n <- dim(X)[1]
        L <- dim(X)[3]
        
        pi <- matrix(0, L, maxJ)
        
        for(l in 1:L) 
            for(j in 1:J[l])
                pi[l, j] <- runif(1, min = 0.1, max = 0.9)
        
        return(pi)
    }
    pi.init.mat <- pi.init
    if(is.null(pi.init)) pi.init.mat <- pi.init.set(X, maxJ, J)
    
    c.init.mat <- c.init
    if(is.null(c.init)) c.init.mat <- rep(0.1, K)

    temp1 <- 1
    temp2 <- 1
    if(popdif.flag) { 
        if(sample.type[4] == 1) temp1 <- num.sample * maxJ * L
        if(sample.type[5] == 1) temp2 <- num.sample * K
    }
    
    res <- .C("structure_JM",
              as.integer(aperm(X, c(3, 2, 1))),
              as.integer(J),
              as.integer(maxJ),
              as.integer(n),
              as.integer(L),
              as.integer(K),
              as.integer(burn.in),
              as.integer(num.sample),
              as.integer(thin),
              as.integer(print.step),
              as.integer(sample.type),
              as.double(alpha.start),
              as.double(alpha.sd),
              as.integer(alpha.step),
              as.double(alpha.max),
              as.integer(fix.alpha),
              as.double(beta),
              as.integer(popdif.flag),
              as.integer(init.z),
              as.integer(aperm(z.init, c(3, 2, 1))),
              as.integer(fix.z),
              as.double(t(pi.init.mat)),
              as.integer(fix.pi),
              as.double(c.init.mat),
              as.double(c.prior.mu),
              as.double(c.prior.sd),
              as.double(c.sd),
              Alpha.sample = double(num.sample),
              Q.mean = double(n * K),
              P.mean = double(K * maxJ * L),
              Z.mean = integer(n * 2 * L),
              pi.sample = double(temp1),
              c.sample = double(temp2),
              PX.sample = double(num.sample),
	PACKAGE = "popgen")
    
    ## output arrays
    Alpha.sample <- res$Alpha.sample
    Z.mean <- aperm(array(res$Z.mean, dim = c(L, 2, n)), c(3, 2, 1))
    P.mean <- aperm(array(res$P.mean, dim = c(L, maxJ, K)), c(3, 2, 1))
    Q.mean <- aperm(array(res$Q.mean, dim = c(K, n)), c(2, 1))

    if(popdif.flag) { 
        pi.sample <- aperm(array(res$pi.sample, dim = c(maxJ, L, num.sample), c(3, 2, 1)))
        c.sample <- aperm(array(res$c.sample, dim = c(K, num.sample), c(2, 1)))
    }
    else {
        pi.sample <- NULL
        c.sample <- NULL
    }

    mu <- mean(res$PX.sample)
    v <- var(res$PX.sample)
    log.px.k <- mu - v / 2
    
    cat(paste("Estimated Ln Prob of Data   =", log.px.k, "\n", sep = "")) 
    cat(paste("Mean value of ln likelihood =", mu, "\n", sep = "")) 
    cat(paste("Variance of ln likelihood   =", v, "\n", sep = ""))
    cat(paste("Mean value of alpha         =", mean(Alpha.sample), "\n", sep = ""))
     
    return(list(Alpha = Alpha.sample,
                Z = Z.mean / num.sample,
                P = P.mean / num.sample,
                Q = Q.mean / num.sample,
                mu = mu,
                v = v,
                log.PX.K = log.px.k,
                PX = res$PX.sample,
                pi = pi.sample,
                c = c.sample)
           )
}
    




nps <- function(X, dmetric = "manhattan", method = "mcquitty", gap.n = 100) {

  require(cluster)
  
  data <- apply(X, c(1, 3), sum)
  
  d <- dist(data, method = dmetric)
  hh <- hclust(d, method = method)
  
  f.W <- function(clust, d) {
      
      n <- nrow(clust$merge) + 1
      W <- rep(0, n)
      D <- as.matrix(d)
      
      for(i in 1:n) {
          cl <- cutree(clust, i)
          for(j in 1:i) {
              d1 <- D[cl == j, cl == j]
              W[i] <- W[i] + sum(d1) / (2 * sum(cl == j))
          }
      }
      
      return(W)
  }
  
  gap <- log(f.W(hh, d))
  
  f.calc.freq <- function(data) {
      
      freq <- matrix(0, 3, ncol(data))
      for(i in 1:ncol(data)) {
          for(j in 0:2) {
              freq[j + 1, i] <- sum(data[, i] == j, na.rm = TRUE)
          }
          freq[, i] <- freq[, i] / sum(freq[, i])
      }
      
      return(freq)
  }
  
  gap.ref <- matrix(0, gap.n, length(gap))
  freq <- f.calc.freq(data)

  f.gen.data.c <- function(freq, n) {

     
      a <- .C("gen_data",
              as.double(freq),
              as.integer(n),
              as.integer(ncol(freq)),
              dat = integer(n * ncol(freq)),
	PACKAGE = "popgen")
      
      new.data <- matrix(a$dat, n, ncol(freq), byrow = TRUE)
      
      return(new.data)
  }
  
  for(i in 1:gap.n) {
      print(i)
      ref.data <- f.gen.data.c(freq, nrow(data))
      d1 <- dist(ref.data, method = dmetric)    
      hh1 <- hclust(d1, method = method)
      gap.ref[i, ] <- log(f.W(hh1, d1))
  }
  
  l <- colSums(gap.ref) / gap.n
  sdk <- sqrt(apply(gap.ref, 2, FUN = var))
  sk <- sdk * sqrt(1 + 1 / gap.n)
  
  gap <- l - gap
  
  G <- rep(0, length(gap))
  for(i in 1:(length(gap) - 1)) {
    G[i] <- (gap[i] > (gap[i + 1] - sdk[i + 1]))
  }
  
  a <- (G == 1) * 1:length(gap)
  a <- a[a != 0]
  num <- min(a, na.rm = TRUE)
  
  return(list(num = num, clust = hh, gap = gap, sk = sk))
}

nps.plot <- function(res, k.max) {
    
    plot(1:k.max, res$gap[1:k.max], xlab = "Number of Clusters", ylab = "Gap")
    for(i in 1:k.max) {
        arrows(i, res$gap[i] + res$sk[i], i ,res$gap[i] - res$sk[i], length = 0.03, angle = 90, code = 3) 
    }
    
    lines(c(res$num, res$num), range(res$gap, na.rm = TRUE), col = 2)
}





















