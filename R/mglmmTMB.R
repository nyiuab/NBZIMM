

mglmmTMB <- function(y, formula, data, family = nbinom2(), 
                     zi = ~0, disp = ~1,
                     min.p = 0, verbose = TRUE)
{
  library(glmmTMB)
  start.time <- Sys.time()
  
  call <- match.call()
  y <- nonzero(y=y, min.p=min.p, sort=FALSE)$y.filter
  if (missing(data)) data <- NULL
  
  if (verbose) cat("Analyzing", NCOL(y), "responses: \n")
  fm <- y.one ~ .
  fm[[3]] <- formula[[2]]
  fit <- vector(mode="list", length=NCOL(y))
  names(fit) <- colnames(y)
  for (j in 1:NCOL(y)){
    y.one <- y[, j]
    data1 <- data.frame(cbind(y.one, data))
    tryCatch( {
      fit[[j]] <- glmmTMB(formula=fm, data=data1, family=family,
                      zi=zi, disp=disp) 
      if (verbose) cat(j, "")
    }, error = function(e) {message("\n", "y", j, " error: ", conditionMessage(e), sep="")} )
  } 
  fit <- fit[!sapply(fit, is.null)]
  responses <- names(fit)
  out <- summary(fit[[1]])$coefficients
  variables <- list(cond=rownames(out$cond))
  if (!is.null(out$zi)) 
    variables$zi <- rownames(out$zi)
  if (!is.null(out$disp)) 
    variables$disp <- rownames(out$disp)

  res <- list(fit=fit, responses=responses, variables=variables, call=call)
  
  class(res) <- "mglmmTMB"
  stop.time <- Sys.time()
  minutes <- round(difftime(stop.time, start.time, units="min"), 3)
  if (verbose) 
    cat("\n Computational time:", minutes, "minutes \n")
  
  res
}

fixed.mglmmTMB <- function(object) 
{
  fit <- object$fit
  out <- lapply(fit, summary)
  for (j in 1:length(out)) out[[j]] <- out[[j]]$coefficients
  
  res <- vector(mode="list", length=length(out[[1]][!sapply(out[[1]], is.null)]))
  names(res) <- names(out[[1]][!sapply(out[[1]], is.null)])
  for (k in 1:length(res)){
    responses <- NULL
    for (j in 1:length(out)){
      res[[k]] <- rbind(res[[k]], out[[j]][[k]])
      responses <- c(responses, rep(names(out)[j], nrow(out[[j]][[k]])))
    }
    variables <- rownames(res[[k]])
    res[[k]] <- data.frame(responses, variables, res[[k]])
    rownames(res[[k]]) <- paste(res[[k]][,1], "--", res[[k]][,2], sep="")
    colnames(res[[k]]) <- c("responses","variables","Estimate","Std.Error",
                            "z.value","pvalue")
  }
  
  for (k in 1:length(res)){
    res0 <- res[[k]]
    res0$padj <- res0$pvalue
    vars <- unique(res0[, 2])
    for(j in 1:length(vars))
    {
      p <- res0[res0[,2]==vars[j], "pvalue"]
      nam <- rownames(res0[res0[,2]==vars[j], ])
      res0[nam, "padj"] <- p.adjust(p, method="fdr")
    }
    res[[k]] <- res0
  }
  
  res
}



