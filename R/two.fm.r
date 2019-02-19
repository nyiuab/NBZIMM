

two.fm <- function (formula, data, na.action, subset, weights) 
{
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "weights", "offset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  if (length(formula[[3]]) > 1 && identical(formula[[3]][[1]], as.name("|"))) {
    ff <- formula
    formula[[3]][1] <- call("+")
    mf$formula <- formula
    ffc <- . ~ .
    ffz <- ~.
    ffc[[2]] <- ff[[2]]
    ffc[[3]] <- ff[[3]][[2]]
    ffz[[3]] <- ff[[3]][[3]]
    ffz[[2]] <- NULL
  }
  else {
    ffz <- ffc <- ff <- formula
    ffz[[2]] <- NULL
  }
  if (inherits(try(terms(ffz), silent = TRUE), "try-error")) {
    ffz <- eval(parse(text = sprintf(paste("%s -", deparse(ffc[[2]])), deparse(ffz))))
  }
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  mtX <- terms(ffc, data = data)
  X <- model.matrix(mtX, mf)
  mtZ <- terms(ffz, data = data)
  mtZ <- terms(update(mtZ, ~.), data = data)
  Z <- model.matrix(mtZ, mf)
  Y <- model.response(mf, "numeric")
  weights <- model.weights(mf)
  if (is.null(weights)) weights <- rep(1, length(Y))
  names(weights) <- rownames(mf)
  offsetx <- model_offset_2(mf, terms = mtX, offset = FALSE)
  if (is.null(offsetx)) offsetx <- rep(0, length(Y))
  offsetz <- model_offset_2(mf, terms = mtZ, offset = FALSE)
  if (is.null(offsetz)) offsetz <- rep(0, length(Y))
  
  list(mf=mf, fmc=ffc, fmz=ffz, Y=Y, X=X, Z=Z, weights=weights, 
       offsetx=offsetx, offsetz=offsetz)
}
