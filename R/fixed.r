

fixed <- function(object)
{
  if ("nbmm" %in% class(object)) 
    res <- fixed.nb(object)
  if (all(class(object)=="lme")) 
    res <- fixed.lme(object)
  if ("zinbmm" %in% class(object) | "zigmm" %in% class(object)) 
    res <- fixed.zi(object)
  if ("mms" %in% class(object)) 
    res <- fixed.mms(object)
  if ("mgam" %in% class(object))
    res <- fixed.mgam(object)
  if ("mglmmTMB" %in% class(object))
    res <- fixed.mglmmTMB(object)
  
  res
}

fixed.nb <- function(object)
{
  if (!"nbmm" %in% class(object)) stop("only for glmm.nb()")
  res <- summary(object)$tTable[, c(1,2,5), drop = FALSE]
  colnames(res) <- c("Estimate", "Std.Error", "pvalue")
  res <- list(dist = res)
}

fixed.lme <- function(object)
{
  if (any(class(object)!="lme")) stop("only for lme()")
  res <- summary(object)$tTable[, c(1,2,5), drop = FALSE]
  colnames(res) <- c("Estimate", "Std.Error", "pvalue")
  res <- list(dist = res)
}

fixed.zi <- function(object)
{
  if (!"zinbmm" %in% class(object) & !"zigmm" %in% class(object)) 
    stop("only for glmm.zinb() or lme.zig()")
  
  dist <- summary(object)$tTable[, c(1,2,5), drop = FALSE]
  colnames(dist) <- c("Estimate", "Std.Error", "pvalue")
  
  if(is.na(object$zi.fit[1])) 
    zero <- matrix(NA, nrow=1, ncol=3)
  else{
    zi.random <- any(class(object$zi.fit)=="lme")
    if (zi.random) zero <- summary(object$zi.fit)$tTable 
    else zero <- summary(object$zi.fit)$coefficients 
    zero <- zero[, c(1, 2, ncol(zero)), drop = FALSE]
  }
  colnames(zero) <- c("Estimate", "Std.Error", "pvalue")
  
  res <- list(dist = dist, zi = zero)
  res
}

fixed.mms <- function(object) 
{
  fit <- object$fit
  fit <- fit[!sapply(fit, is.null)]
  
  if (class(object)[2]=="nb") out <- lapply(fit, fixed.nb)
  if (class(object)[2]=="lme") out <- lapply(fit, fixed.lme)
  if (class(object)[2] %in% c("zinb","zig")) out <- lapply(fit, fixed.zi)
  
  res <- vector(mode="list", length = length(out[[1]]))
  names(res) <- names(out[[1]])
  for (k in 1:length(res)){
    responses <- NULL
    for (j in 1:length(out)){
      res[[k]] <- rbind(res[[k]], out[[j]][[k]])
      responses <- c(responses, rep(names(out)[j], nrow(out[[j]][[k]])))
    }
    variables <- rownames(res[[k]])
    res[[k]] <- data.frame(responses, variables, res[[k]])
    rownames(res[[k]]) <- paste(res[[k]][,1], "--", res[[k]][,2], sep="")
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

