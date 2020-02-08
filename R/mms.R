
mms <- function(y, fixed, random, data, method=c("nb", "lme", "zinb", "zig"),
                correlation, zi_fixed= ~1, zi_random=NULL, 
                niter=30, epsilon=1e-05, min.p=0, verbose=TRUE)
{
  if (!requireNamespace("nlme")) install.packages("nlme")
  library(nlme)
  start.time <- Sys.time()
  if (missing(data)) stop("'data' should be specified")
  
  call <- match.call()
  method <- method[1]
  y <- nonzero(y=y, min.p=min.p, sort=FALSE)$y.filter
  if (missing(correlation)) correlation <- NULL
  
  if (verbose) cat("Analyzing", NCOL(y), "responses: \n")
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
        fit[[j]] <- glmm.zinb(fixed = fm, random = random, data = data1, correlation = correlation, 
                              zi_fixed = zi_fixed, zi_random = zi_random, 
                              niter = niter, epsilon = epsilon, verbose = FALSE) 
      if (method == "zig")
        fit[[j]] <- lme.zig(fixed = fm, random = random, data = data1, correlation = correlation, 
                            zi_fixed = zi_fixed, zi_random = zi_random, 
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


mms.GLMMadaptive <- function(y, fixed, random, data, family,
                            zi_fixed=NULL, zi_random=NULL, penalized=FALSE,
                            min.p=0, verbose=TRUE)
{
  if (!requireNamespace("GLMMadaptive")) install.packages("GLMMadaptive")
  library(GLMMadaptive)
  start.time <- Sys.time()
  if (missing(data)) stop("'data' should be specified")
  
  call <- match.call()
  y <- nonzero(y=y, min.p=min.p, sort=FALSE)$y.filter
  
  if (verbose) cat("Analyzing", NCOL(y), "responses: \n")
  fm <- y.one ~ .
  fm[[3]] <- fixed[[2]]
  fit <- vector(mode="list", length=NCOL(y))
  names(fit) <- colnames(y)
  for (j in 1:ncol(y)){
    y.one <- y[, j]
    data1 <- as.data.frame(cbind(y.one, data))
    tryCatch( {
      fit[[j]] <- mixed_model(fixed=fm, random=random, data=data1, family=family, 
                              zi_fixed=zi_fixed, zi_random=zi_random, penalized=penalized) 
      if (verbose) cat(j, "")
      
    }, error = function(e) {cat("\n", "y", j, " error: ", conditionMessage(e), sep="")} )
  } 
  fit <- fit[!sapply(fit, is.null)]
  
  responses <- names(fit)
  out <- summary(fit[[1]])
  variables <- list(dist=rownames(out$coef_table))
  if (!is.null(out$coef_table_zi)) 
    variables$zero <- rownames(out$coef_table_zi)
  
  res <- list(fit=fit, responses=responses, variables=variables, call=call)
  class(res) <- c("mms", class(fit[[1]]))
  
  stop.time <- Sys.time()
  minutes <- round(difftime(stop.time, start.time, units="min"), 3)
  if (verbose) 
    cat("\n Computational time:", minutes, "minutes \n")
  
  res
}

