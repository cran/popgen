
LDmat = function (mat, typ = c("genotype", "haplotype"), pos = NULL, plotmat = TRUE)
{

    s <- ncol(mat)
    res <- matrix(1.0, s, s)
    
    typ <- match.arg(typ)

    
    flag = 0
    if(typ == "genotype")
      flag = 1
    
    ans <- .C("ld_JM",
              as.integer(nrow(mat)),
              as.integer(ncol(mat)),
              as.integer(flag),
              as.integer(t(mat)),
              m = as.double(res),
              PACKAGE = "popgen")$m

    dim(ans) = c(s, s)
    ans = t(ans)
    
    if(plotmat) {
      if(is.null(pos)) pos = seq(0, 1, len = ncol(mat))
        cols <- c("black", grey(0.4), "blue", "green", "orange", "red")
        par(mar = c(5, 5, 4, 2) + 0.1)
        plot(c(0, 1.1) * max(pos), c(0, 1) * max(pos), typ = "n", xlab = "", ylab = "", axes = FALSE, cex.lab = 2)
        if(is.null(pos)) pos <- seq(0, 1, len = nrow(res))
        image(pos, pos, ans, col = cols, add = TRUE)
        axis(1)
        axis(2)
        text(0.50 * max(pos), -max(pos) / 6, labels = expression(r^2), cex = 2, xpd = NA)
        text(-max(pos)/ 6, 0.50 * max(pos), labels = "D'", cex = 2, xpd = NA)

        image(c(1.060, 1.085) * max(pos), seq(0, 1, len=101) * max(pos),
              matrix(seq(0, 1, len = 100), 1), add = TRUE, col = cols)
        for(n in seq(0, 1, 0.2)) {
            expr <- substitute(expression(n), list(n=n))
            mode(expr) <- "expression"
            text(1.09 * max(pos), n * max(pos), expr, adj=c(0, 0.5))
        }
        lines((min(pos) - 1):(max(pos) + 1), (min(pos) - 1):(max(pos) + 1))

        par(mar = c(5, 4, 4, 2) + 0.1)
        
      }
    
    return(ans)
  }
