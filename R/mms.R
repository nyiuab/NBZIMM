
mms <- function(y, fixed, random, data, method = c("nb", "lme", "zinb", "zig"),
                correlation, zi.random = FALSE, niter = 30, epsilon = 1e-05,  
                min.p = 0, sort = FALSE, verbose = TRUE)
{
  if (!requireNamespace("nlme")) install.packages("nlme")
  library(nlme)
  start.time <- Sys.time()
  
  call <- match.call()
  method <- method[1]
  y <- as.matrix(y)
  if (is.null(rownames(y))) rownames(y) <- paste("v", 1:nrow(y), sep = "")
  if (is.null(colnames(y))) colnames(y) <- paste("y", 1:ncol(y), sep = "")
  nonzero.p <- apply(y, 2, function(z) {length(z[z != 0])/length(z)} )
  if (sort){
    nonzero.p <- sort(nonzero.p, decreasing = T)
    y <- y[, names(nonzero.p), drop = FALSE]
  }
  if (ncol(y) > 1 & min.p > 0){
    sub <- names(nonzero.p[nonzero.p > min.p])
    if (length(sub) == 0) stop("min.p is too large: no response") 
    else y <- y[, sub, drop = FALSE]
  }
  if (missing(data)) data <- NULL
  if (missing(correlation)) correlation <- NULL
  
  if (verbose) cat("Analyzing", ncol(y), "responses: \n")
  fm <- y.one ~ .
  fm[[3]] <- fixed[[2]]
  fit <- vector(mode="list", length=NCOL(y))
  names(fit) <- colnames(y)
  for (j in 1:ncol(y)){
    y.one <- y[, j]
    data1 <- as.data.frame(cbind(y.one, data))
    tryCatch( {
      if (method == "nb")
        fit[[j]] <- glmm.nb(fixed = fm, random = random, data = data1, correlation = correlation,
                            niter = niter, epsilon = epsilon, verbose = FALSE) 
      if (method == "lme")
        fit[[j]] <- lme(fixed = fm, random = random, data = data1, correlation = correlation) 
      if (method == "zinb")
        fit[[j]] <- glmm.zinb(fixed = fm, random = random, data = data1, correlation = correlation, zi.random = zi.random, 
                            niter = niter, epsilon = epsilon, verbose = FALSE) 
      if (method == "zig")
        fit[[j]] <- lme.zig(fixed = fm, random = random, data = data1, correlation = correlation, zi.random = zi.random, 
                              niter = niter, epsilon = epsilon, verbose = FALSE) 
      if (verbose) cat(j, "")
      
    }, error = function(e) {cat("\n", "y", j, " error: ", conditionMessage(e), sep="")} )
  } 
  fit <- fit[!sapply(fit, is.null)]
  
  responses <- names(fit)
  dist <- colnames(coef(fit[[1]]))
  variables <- list(dist=dist)
  if (method %in% c("zinb", "zig")) variables$zero <- colnames(fit[[1]]$xz)
  
  res <- list(fit = fit, responses = responses, variables = variables, call = call)
  class(res) <- c("mms", method)
  
  stop.time <- Sys.time()
  minutes <- round(difftime(stop.time, start.time, units = "min"), 3)
  if (verbose) 
    cat("\n Computational time:", minutes, "minutes \n")
  
  res
}

