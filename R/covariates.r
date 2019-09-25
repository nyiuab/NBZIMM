
# standardize covariates

covariates <- function (x.con, x.cat, con.rescale = TRUE, cat.center = FALSE, 
                        fill.missing = TRUE, ind.group = NULL)
{
# x1: continuous variables; x2: categorical variables
  if (missing(x.con)) x.con <- NULL
  if (missing(x.cat)) x.cat <- NULL
  x1 <- x.con
  x2 <- x.cat
  if (!is.null(x1)) {
    x1 <- as.matrix(x1)
    if(is.null(colnames(x1))) {
      warning("some variables have no names.", call. = FALSE)
      colnames(x1) <- paste("x.con", 1:ncol(x1), sep = "")
    }
  }
  if (!is.null(x2)) {
    x2 <- as.matrix(x2)
    if(is.null(colnames(x2))) {
      warning("some variables have no names.", call. = FALSE)
      colnames(x2) <- paste("x.cat", 1:ncol(x2), sep = "")
    }
  }
  
  x <- cbind(x1, x2)
  func <- function(d){
    na.prop <- length(which(is.na(d)))/length(d)
    return(na.prop)
  }
  na.prop <- apply(x, 2, func)
  if(any(na.prop > 20/100)) {
    d <- which(na.prop > 20/100)
    warning("the following ", length(d), " variables have more than 20% missing values:", call. = FALSE)
    for(j in 1:length(d)) warning(names(d)[j], call. = FALSE)
  }
  
  if (!is.null(x1) & con.rescale) x1 <- scale(x1)
  
  if (!is.null(x2)) {
    X2 <- NULL
    for (i in 1:ncol(x2)) {
      x <- x.new <- x2[, i]
      x <- as.factor(x)
      f <- ~ x - 1
      mf <- model.frame(f, na.action = na.pass)
      x.new <- model.matrix(f, mf)[, -1, drop = FALSE]
      colnames(x.new) <- paste(colnames(x2)[i], levels(x)[-1], sep = "_")
      if (ncol(x.new) > 1) {
        rsums <- rowSums(x.new)
        if (all(rsums[!is.na(rsums)] == 1)) x.new <- x.new[, -1, drop = FALSE]
      }
      X2 <- cbind(X2, x.new)
    }
    if (cat.center) {
      x2 <- scale(X2, scale = FALSE)
    }
    else x2 <- X2
  }
  
  x <- cbind(x1, x2)
  
  w <- apply(x, 2, var, na.rm = TRUE)
  if (length(is.na(w))){
    w <- w[!is.na(w)]
    x <- x[, names(w), drop = FALSE]
  }
  if (length(which(w == 0))) {
    x <- x[, which(w != 0), drop = FALSE]
    d <- which(w == 0)
    warning(length(d), " variables with no-variation are removed!", call. = FALSE )
  }
  
  if (fill.missing & any(is.na(x))) {
    if(!is.null(ind.group)) {
      ind.group <- as.factor(ind.group)
      if(nrow(x) != length(ind.group)) {
        warning("x and ind.group have different obs. Cannot use group information!", call. = FALSE, immediate. = TRUE)
        ind.group <- NULL
      }
    }
    warning(sum(is.na(x)) , " missing values have been filled!", call. = FALSE)
    if(is.null(ind.group)) {
      func = function(w){
        na.index <- which(is.na(w))
        w[na.index] <- mean(w, na.rm = TRUE)
        return(w)
      }
      x <- apply(x, 2, func)
    }
    if(!is.null(ind.group)){
      func <- function(w){
        group.means <- tapply(w, ind.group, mean, na.rm = TRUE)
        group.means <- ifelse(is.na(group.means), mean(w, na.rm = TRUE), group.means)
        for(k in 1:length(group.means)){
          na.index <- which(is.na(w) & ind.group == names(group.means)[k])
          w[na.index] <- mean(w[ind.group == names(group.means)[k]], na.rm = TRUE)
        }
        return(w)
      }
      x <- apply(x, 2, func)
    }
  }
  
  X <- data.frame(x)
  colnames(X) <- colnames(x)  
  X
}


#*******************************************************************************

