# Quantification type I
qt1 <- function (dat, # data matrix 
                 y, # dependent variable
                 func.name = c ("solve", "ginv")) # Selection of function name taking inverse matrix
{
  vname <- colnames(dat) # variable name
  vname.y <- deparse(substitute(y))
  cname <- unlist(sapply(dat, levels)) # category name
  dat <- data.frame(dat, y)
  dat <- subset(dat, complete.cases (dat)) # except cases with missing values
  p <- ncol(dat) # The rightmost is the dependent variable, the rest are the categorical variables
  ncat <- p-1 # The rightmost end is the dependent variable, the rest are the categorical variables
  stopifnot(all (sapply (dat [, 1: ncat], is.factor))) # All independent variables must be factor
  dat [, 1: ncat] <-lapply (dat [, 1: ncat, drop = FALSE], as.integer)
  nc <- nrow (dat) # number of columns
  mx <- sapply (dat [, 1: ncat, drop = FALSE], max) # Maximum value of each category variable
  start <- c (0, cumsum (mx) [-ncat]) # total sequence number
  nobe <- sum(mx) # Total number of categories
  
  # Convert to data matrix using dummy variables
  x <- t(apply (dat, 1, function (obs) {
    zeros <- numeric(nobe)
    zeros [start + obs [1: ncat]] <- 1
    c (zeros [-start-1], obs [ncat + 1])
  }))
  
  # Solve as multiple regression analysis
  a <- cov(x)
  ndim <- nobe-ncat
  if(match.arg(func.name) == "solve") {
    inverse <- solve
    B <- inverse(a[1: ndim, 1: ndim], a[ndim + 1, 1: ndim])
  } else {
    
    library (MASS)
    inverse <- ginv
    B <- inverse(a[1: ndim, 1: ndim]) %*% a[ndim + 1, 1: ndim]
  }
  
  m <- colMeans (x)
  const <- m[ndim + 1] - sum (B * m [1: ndim])
  prediction <- x[, 1: ndim] %*% as.matrix(B) + const
  observed <- x[, ndim + 1]
  prediction <- cbind(observed, prediction, observed-prediction)
  
  # Quantification Convert to solution as class I
  ncase <- nrow(dat)
  s <- colSums(x)
  name <- coef <- NULL
  en <- 0
  for (i in 1: ncat) {
    st <- en + 1
    en <- st + mx [i]-2
    target <- st:en
    temp.mean <- sum(s[target] * B[target]) / ncase
    const <- const + temp.mean
    coef <- c(coef, -temp.mean, B [target] -temp.mean)
  }
  coef <- c(coef, const)
  names(coef) <- c(paste(rep(vname, mx), cname, sep = "."), "constant term")
  
  # Partial correlation coefficient between item variable and dependent variable
  par <- matrix (0, nrow = nc, ncol = ncat)
  for(j in 1: nc) {
    en <- 0
    for (i in 1: ncat) {
      st <- en + 1
      en <- st + mx [i]-2
      target <- st:en
      par[j, i] <- crossprod(x[j, target], B[target])
    }
  }
  
  par <- cbind(par, observed)
  r <- cor(par)
  print(vname)
  i <- inverse(r)
  d <- diag(i)
  partial.cor <- (-i / sqrt (outer (d, d))) [ncat + 1, 1: ncat]
  partial.t <- abs(partial.cor) * sqrt ((nc-ncat-1) / (1-partial.cor ^ 2))
  partial.p <- pt(partial.t, nc-ncat-1, lower.tail = FALSE) * 2
  partial <- cbind (partial.cor, partial.t, partial.p)
  
  coef <- as.matrix (coef)
  colnames(coef) <- "category score"
  colnames(prediction) <- c ("observed value", "predicted value", "residual")
  colnames(partial) <- c("partial correlation coefficient", "t value", "P value")
  rownames(prediction) <- paste ("#", 1: nc, sep = "")
  rownames(partial) <- vname
  rownames(r) <- colnames(r) <- c(vname, vname.y)
  return(structure (list (coefficients = as.matrix (coef),
                          r = r, partial = partial, prediction = prediction), class = "qt1"))
}
# print method
print.qt1 <- function(obj, # object returned by qt1
                      digits = 5) # Number of displayed digits
{
  print (round (obj $ coefficients, digits = digits))
}
# summary method
summary.qt1 <-function (obj, # qt1 returns the object
                        digits = 5) # Number of displayed digits
{
  print.default (obj, digits = digits)
}
# plot method
plot.qt1 <-function (obj, # qt1 returns object
                     which = c("category.score", "fitness"), ...)# Choose whether to show category scores or show observation and forecast values # Arguments passed to barplot, plot
{
  if (match.arg (which) == "category.score") {
    coefficients <- obj $ coefficients [-length (obj $ coefficients),]
    coefficients <- rev (coefficients)
    cname <- names(coefficients)
    names(coefficients) <-NULL
    barplot(coefficients, horiz = TRUE, xlab = "category score", ...)
    text(0, 1.2 * (1: length (cname) -0.5), cname, pos = ifelse (coefficients> 0, 2, 4))
  }
  else {
    result <-obj $ prediction
    plot (result [, 2], result [, 1], xlab = "predicted value", ylab = "observed value", asp = 1, ...)
    abline (c (0,1))
  }
}