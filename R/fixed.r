

fixed <- function(object)
{
  if ("nbmm" %in% class(object)) res <- fixed.nb(object)
  if (all(class(object)=="lme")) res <- fixed.lme(object)
  if ("zinbmm" %in% class(object) | "zigmm" %in% class(object)) res <- fixed.zi(object)
  if ("mms" %in% class(object)) res <- fixed.mms(object)
  res
}

fixed.nb <- function(object)
{
  if (!"nbmm" %in% class(object)) stop("only for glmm.nb()")
  
  res <- summary(object)$tTable[, c(1,2,5), drop = FALSE]
  res[, c(1,2)] <- round(res[, c(1,2)], digits = 3)
  res[, 3] <- signif(res[, 3], 2)
  colnames(res) <- c("Estimate", "Std.Error", "pvalue")
  
  res <- list(dist = res)
}

fixed.lme <- function(object)
{
  if (any(class(object)!="lme")) stop("only for lme()")
  
  res <- summary(object)$tTable[, c(1,2,5), drop = FALSE]
  res[, c(1,2)] <- round(res[, c(1,2)], digits = 3)
  res[, 3] <- signif(res[, 3], 2)
  colnames(res) <- c("Estimate", "Std.Error", "pvalue")
  
  res <- list(dist = res)
}

fixed.zi <- function(object)
{
  if (!"zinbmm" %in% class(object) & !"zigmm" %in% class(object)) 
    stop("only for glmm.zinb() or lme.zig()")
  
  dist <- summary(object)$tTable[, c(1,2,5), drop = FALSE]
  dist[, c(1,2)] <- round(dist[, c(1,2)], digits = 3)
  dist[, 3] <- signif(dist[, 3], 2)
  colnames(dist) <- c("Estimate", "Std.Error", "pvalue")
  
  if(is.na(object$fit.zero[1])) {
    zero <- matrix(NA, nrow=ncol(object$xz), ncol=3)
    rownames(zero) <- colnames(object$xz)
  }
  else{
    zi.random <- any(class(object$fit.zero)=="lme")
    if (zi.random) zero <- summary(object$fit.zero)$tTable 
    else zero <- summary(object$fit.zero)$coefficients 
    zero <- zero[, c(1, 2, ncol(zero)), drop = FALSE]
    zero[, c(1,2)] <- round(zero[, c(1,2)], digits = 3)
    zero[, 3] <- signif(zero[, 3], 2)
  }
  colnames(zero) <- c("Estimate", "Std.Error", "pvalue")
  
  res <- list(dist = dist, zero = zero)
  res
}

fixed.mms <- function(object) 
{
  if (all(class(object)!="mms")) stop("only for mms()")
  fit <- object$fit
  fit <- fit[!sapply(fit, is.null)]
  
  if (any(class(object)=="nb")) out <- lapply(fit, fixed.nb)
  if (any(class(object)=="lme")) out <- lapply(fit, fixed.lme)
  if (any(class(object) %in% c("zinb", "zig"))) out <- lapply(fit, fixed.zi)

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
  
  res
}

