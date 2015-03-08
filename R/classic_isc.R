#' Calculate item reliabilty correlation from full data.
#' Note that this calculates the correlation of each column against the others.
#'
#' @param M data.frame, matrix, or lavaan fit object.
#'        if a data.frame or matrix is given, calculate item-total correlation
#'        (using one against the others). If a lavaan fit object is given, return
#'        a data.frame with mean factor loadings by factor.
#' @param is.corr boolean indicating if M is a corraltion matrix
#' @return array with item-total correlation for each column
#' @export
item.rel = function(M, is.corr=FALSE){
  if (is.corr){
    # M is a correlation matrix
    covttl = sum(M)
    covcolsums = colSums(M)
    (covcolsums-1) / sqrt(covttl - 2*covcolsums + 1)   # +1 because subtract own var=1 twice
  }
  else if (class(M) == 'lavaan') {
    # M is lavaan object (probably more R-saavy ways to do this)
    pars = parameterEstimates(M, standardized = TRUE)
    loadings = subset(pars, op == '=~') #only factor loadings
    newnames = c(std.all='rel', lhs='factor', rhs='item')
    rename(loadings, newnames)[,newnames]
    
    #ddply(loadings, .(lhs), summarize, isc = mean(est))  
  }
  else {
    # use means to calculate item-total correlation
    ttl = rowSums(M)
    out = laply(colnames(M), function(name){
      #data.frame(sub=name, rel=cor(M[,name], ttl - M[,name]))
      rel=cor(M[,name], ttl-M[,name])
    })
    names(out) = colnames(M)
    out
  }
}

#' Calculate inter-item correlation (internal consistency).
#' Note that this returns all of the relevant correlations (then you can average them).
#' 
#' @param M data.frame or matrix
#' @param M2 optional data.frame or matrix to cross-correlate
#' @return array of all inter-item correlations
#' @export

item.ic = function(M, M2=NULL, lower=TRUE){
  # inter-item correlation (average off-diag corr)
  if (class(M) == 'lavaan'){
    lam = inspect(M)$lambda         # dataframe of indices to items for each factor
    Sigma = fitted(M)$cov
    ic = alply(lam, 2, function(arr) {
      indx = which(arr != 0)
      Sigma[indx,indx]
      })
    names(ic) = attr(ic, 'split_labels')$X1
    if (lower) return(lapply(ic, function(M) M[lower.tri(M)]))
    else return(ic)
  }
  if (is.null(M2)){
    cors = cor(M)
    if (lower) return(cors[lower.tri(cors)])
    else return(cors)
  }
  else {
    return(c(cor(M, M2)))
  }
}

#' Estimate the reliability for latent factors underlying items
#' Note that if it is a lavaan object, it will return a ___ with reliability 
#' for each factor. If it is a data.frame, it will assume a single factor model,
#' and approximate reliability from item-total correlations.
group.rel = function(M1){
  if (class(M1) == 'lavaan'){
    # latent variable reliability
    reliability(M1)
  }
}

#' Estimate the similarity between two groups
#' 
#' @param M1 First group. May be data.frame, matrix, or lavaan fit object.
#' @param M2 Second group. (uneccessary if using lavaan object)
#' @export
group.cor = function(M1, M2){
  if (class(M1) == 'lavaan'){
    # latent variable correlations
    pars = parameterEstimates(M1, standardize=TRUE)
    latent.vars = unique(pars[pars$op == "=~", 'lhs'])
    latent.pars = subset(pars, op == '~~' & lhs %in% latent.vars)    # get only latent variables
    # make sure variance on latent variables is 1
    stopifnot(
      all(subset(latent.pars, lhs == rhs)$est == 1)
      )
    return(subset(latent.pars, lhs != rhs)
           )
  }
  else {
    # mean-mean correlation
    return(cor(rowMeans(M1), rowMeans(M2)))
  }
}