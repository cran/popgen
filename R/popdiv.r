





popdiv <- function(L, P, NUMA, N, X, burnin = 10000, iter = 10000, m = 300, csd = 0.01, outc = FALSE) {
    
    
    if(outc == 1) lc = P * (iter)
    max.alleles = max(NUMA)
    
    a <- .C("popdiv",
            as.integer(L),
            as.integer(P),
            as.integer(NUMA),
            as.double(t(N)),
            as.double(t(X)),
            as.integer(burnin),
            as.integer(iter),
            as.integer(max.alleles),
            as.double(m),
            as.double(csd),
            as.integer(outc),
            cc = double(lc),
            muc = double(P),
            sdc = double(P),
            mup = double(L * max.alleles),
            sdp = double(L * max.alleles),
            prate = double(1),
            crate = double(1),
	PACKAGE = "popgen")

    C.sample <- matrix(a$cc, iter, P, byrow = TRUE)
    muc = a$muc
    sdc = a$sdc
    mup = matrix(a$mup, L, max.alleles, byrow = TRUE)
    sdp = matrix(a$sdp, L, max.alleles, byrow = TRUE)
    prate = a$prate / ((iter + burnin) * L)
    crate = a$crate / ((iter + burnin) * P)

    
    return(list(C.sample = C.sample, muc = muc, sdc = sdc, mup = mup, sdp = sdp, prate = prate, crate = crate))
    
}

    

popdiv.convert <- function(x, miss = 0) {

    if(!is.matrix(x)) stop(message = "x is not a matrix")
        
    P <- max(x[, 1])
    L <- ncol(x) - 1
    n <- nrow(x)
    
    NUMA <- matrix(0, 1, L)
    
    for(l in 1:L) {
        a <- unique(x[, l + 1])
        a <- a[a != miss]
        NUMA[l] <- length(a)
    }

    max.alleles <- max(NUMA)

    N <- matrix(0, L, P)
    X <- matrix(0, L, P * max.alleles)

    for(l in 1:L) {
       
        a <- unique(x[, l + 1])
        a <- a[a != miss]
        a <- sort(a)
        
        for(i in 1:n) {
            p <- x[i, 1]
            j <- x[i, l + 1]
            j <- sum((j == a) * 1:length(a))

            if(j != 0) {
                N[l, p] =  N[l, p] + 1
                X[l, (p - 1) * NUMA[l] + j] = X[l, (p - 1) * NUMA[l] + j] + 1
            }
        }
    }

    return(list(NUMA = NUMA, N = N, X = X))
}
            
