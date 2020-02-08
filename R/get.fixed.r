

get.fixed <- function(object, part=c("dist", "zero"), vr.name, sort=FALSE)
{
  part <- part[1]
  if (class(object)[1] != "mms") stop("only for object from 'mms'")  
  res <- object$responses
  var <- object$variables
  if (!vr.name %in% c(res, unlist(var))) stop("wrong name given")
  ss <- fixed(object)[[part]]
  if (vr.name %in% res) {
    out <- ss[ss[, 1]==vr.name, ]
    rownames(out) <- out[, 2]  
  }
  else{
    out <- ss[ss[, 2]==vr.name, ]
    rownames(out) <- out[, 1]
  }
  out <- out[, 3:6]
  out <- out[rownames(out)!="(Intercept)", ]
  out <- as.matrix(out)
  if (sort) out <- out[names(sort(out[, "pvalue"])), ]
  
  out
}


