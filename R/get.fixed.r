

get.fixed <- function(object, part=c("dist", "zero"), vr.name, adj.p=FALSE, sort.p=FALSE)
{
  part <- part[1]
  if (class(object)[1] != "mms") stop("only for object from 'mms'")  
  if (!(class(object)[2] %in% c("zinb", "zig"))) part <- "dist"
  res <- object$responses
  var <- object$variables
  if (!vr.name %in% c(res, unlist(var))) stop("wrong name given")
  ss <- fixed(object, adj.p=adj.p)[[part]]
  if (vr.name %in% res) {
    out <- ss[ss[, 1]==vr.name, ]
    rownames(out) <- out[, 2]  
  }
  else{
    out <- ss[ss[, 2]==vr.name, ]
    rownames(out) <- out[, 1]
  }
  out <- out[, 3:5]
  out <- out[rownames(out)!="(Intercept)", ]
  out <- as.matrix(out)
  if (sort.p) out <- out[names(sort(out[, "pvalue"])), ]
  
  out
}


