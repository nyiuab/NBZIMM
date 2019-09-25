

plot.fixed <- function(res, threshold = 0.05, main = " ", col.pts = "black", gap = 0, 
                    show.all.vars = FALSE, show.pvalues = TRUE, 
                    cex.main = 0.9, xlim = NULL, cex.axis = 0.8, cex.var = 0.8, cex.pts = 1, 
                    pch.pts = 20, type = "p", lwd = 1, lty = 1, line = 0,  
                    add = FALSE)
{
  res <- as.matrix(res)
  if (is.null(rownames(res))) {
    rownames(res) <- paste("v", 1:nrow(res), sep = "")
    warning("Rownames are not given !")
  }
  coefs <- res[, 1]
  sds <- res[, 2]
  pvalues <- res[, 3]
  
  n <- length(coefs)
  coef.l <- coefs - 2*sds
  coef.h <- coefs + 2*sds
  coefs <- coefs[n:1]
  coef.l <- coef.l[n:1] 
  coef.h <- coef.h[n:1]
  pvalues <- pvalues[n:1]
  
  p <- signif(pvalues, 2)
  varnames <- names(coefs)
  if (!show.all.vars) varnames <- ifelse(p <= threshold, varnames, "")
  if (length(col.pts) == 1) col.pts[2] <- "gray"

  z <- c(1:n)
  if (n > 1)
    for(i in 2:n) z[i] <- ifelse(p[i] <= threshold, z[i-1] + 1 + gap, z[i-1] + 1)
  if (is.null(xlim)) xlim = c(min(coef.l), max(coef.h))
  if (!add) plot(coefs, z, main = main, cex.main = cex.main, xlab = "", xlim = xlim,
                 ylab = "", yaxt = "n", type = "n", frame.plot = FALSE, cex.axis = cex.axis)
  axis(side = 3, cex.axis = cex.axis)
  lines(c(0, 0), c(-10, max(z) + 10), lty = 2, lwd = 1)
  col <- rep(col.pts[1], n)
  for (i in 1:n) {
    col[i] <- ifelse(p[i] <= threshold, col.pts[1], col.pts[2])
    lines(c(coef.l[i], coef.h[i]), c(z[i], z[i]), col = col[i], lwd = lwd, lty = lty)
    cex <- ifelse(p[i] <= threshold, cex.var, cex.var - 0.1)
    if (!add) mtext(paste(varnames[i]), side = 2, at = z[i], cex = cex, las = 1, col = col[i], line = line)
    cex <- ifelse(p[i] <= threshold, cex.var - 0.1, cex.var - 0.2)
    if(show.pvalues & p[1] > 0 & varnames[i] != "")
      mtext(p[i], side = 4, at = z[i], cex = cex, las = 1, col = col[i], line = line/2)
  }
  points(coefs, z, pch = pch.pts, cex = cex.pts, col = col, type = type, lwd = lwd, lty = lty)

}



