

mgam <- function(y, formula, data, family=nb(), paraPen=NULL,
                 min.p=0, verbose=TRUE)
{
  library(mgcv)
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
      fit[[j]] <- gam(formula=fm, data=data1, family=family, method="REML", paraPen=paraPen) 
      if (verbose) cat(j, "")
    }, error = function(e) {message("\n", "y", j, " error: ", conditionMessage(e), sep="")} )
  } 
  fit <- fit[!sapply(fit, is.null)]
  responses <- names(fit)
  out <- summary(fit[[1]])
  p.variables <- rownames(out$p.table)
  s.variables <- rownames(out$s.table)
  variables <- list(p.variables=p.variables, s.variables=s.variables)

  res <- list(fit=fit, responses=responses, variables=variables, 
              call=call)
  
  class(res) <- "mgam"
  stop.time <- Sys.time()
  minutes <- round(difftime(stop.time, start.time, units="min"), 3)
  if (verbose) 
    cat("\n Computational time:", minutes, "minutes \n")
  
  res
}

fixed.mgam <- function(object) 
{
  fit <- object$fit
  resp <- object$responses
  p.vars <- object$variables$p.variables
  s.vars <- object$variables$s.variables
  out <- lapply(fit, summary)
  
  # parametric part
  p.table <- p.resp <- NULL
  if(!is.null(p.vars)) 
  {
    for (j in 1:length(out))
    {
      p.table <- rbind(p.table, out[[j]]$p.table)
      p.resp <- c(p.resp, rep(resp[j], length(p.vars)))
    }
    p.table <- data.frame(p.resp, rownames(p.table), p.table)
    rownames(p.table) <- paste(p.table[,1], "--", p.table[,2], sep="")
    colnames(p.table) <- c("responses","variables","Estimate",
                          "Std.Error","z.value","pvalue")

    p.table$padj <- p.table$pvalue
    for(j in 1:length(p.vars))
    {
      p <- p.table[p.table[,"variables"]==p.vars[j], "pvalue"]
      nam <- rownames(p.table[p.table[,"variables"]==p.vars[j], ])
      p.table[nam, "padj"] <- p.adjust(p, method="fdr")
    }
  }
  
# smooth part
  s.table <- s.resp <- NULL
  if(!is.null(s.vars)) 
  {
    for (j in 1:length(out))
    {
      s.table <- rbind(s.table, out[[j]]$s.table)
      s.resp <- c(s.resp, rep(resp[j], length(s.vars)))
    }
    s.table <- data.frame(s.resp, rownames(s.table), s.table)
    rownames(s.table) <- paste(s.table[,1], "--", s.table[,2], sep="")
    colnames(s.table) <- c("responses","variables","edf",
                           "Ref.df","Chi.sq","pvalue")
     
    s.table$padj <- s.table$pvalue
    for(j in 1:length(s.vars))
    {
      p <- s.table[s.table[,"variables"]==s.vars[j], "pvalue"]
      nam <- rownames(s.table[s.table[,"variables"]==s.vars[j], ])
      s.table[nam, "padj"] <- p.adjust(p, method="fdr")
    }
  }
  
  out <- list(parametric_terms = p.table, smooth_terms = s.table)
  return(out)
}




