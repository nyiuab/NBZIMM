

nonzero <- function(y, total, min.p=0, sort=TRUE, plot=FALSE)
{
  y <- as(y, "matrix")
  if (is.null(colnames(y))) colnames(y) <- paste("y", 1:ncol(y), sep = "")
  if(missing(total)) total <- rowSums(y)
  total <- unlist(total)
  nonzero.p <- apply(y, 2, function(z) {length(z[z != 0])/length(z)} )
  if (sort) nonzero.p <- sort(nonzero.p, decreasing=T)
  sub <- names(nonzero.p[nonzero.p > min.p])
  if (length(sub) == 0) stop("min.p is too large") 
  else y0 <- y[, sub, drop=FALSE]
  total.mean <- round(mean(total), 2)
  total.sd <- round(sd(total), 2)
  if (plot){
    par(mfrow = c(1, 2), mar = c(5, 4, 4, 4))
    hist(total, nclass=1000, xlab="Total reads", main="")
    zero <- sort(1 - nonzero.p, decreasing=F)
    plot(x=1:length(zero), y=zero, xlab="Taxa", ylab="Zero Proportion")
  }
  list(nonzero.p=nonzero.p, total=total, total_mean_sd=c(total.mean, total.sd), y.filter=y0)
}

