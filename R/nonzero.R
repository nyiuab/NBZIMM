
nonzero <- function(y, total, plot = FALSE)
{
  y <- as.matrix(y)
  if(missing(total)) total <- rowSums(y)
  total <- unlist(total)
  nonzero.p <- apply(y, 2, function(z) {length(z[z != 0])/length(z)} )
  nonzero.p <- sort(nonzero.p, decreasing = T)
  total.mean <- round(mean(total), 2)
  total.sd <- round(sd(total), 2)
  if (plot){
    par(mfrow = c(1, 2), mar = c(5, 4, 4, 4))
    hist(total, nclass=1000, xlab="Total reads", main="")
    plot(x=1:length(nonzero.p), y=1-nonzero.p, xlab="Taxa", ylab="Zero Proportion")
  }
  list(nonzero.p = nonzero.p, total = total, total_mean_sd = c(total.mean, total.sd))
}

