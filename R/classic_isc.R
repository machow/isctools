#' Calculate item-total correlation from full data.
#' Note that this calculates the correlation of each column against the others.
#'
#' @param M data.frame or matrix
#' @param is.corr boolean indicating if M is a corraltion matrix
#' @return array with item-total correlation for each column
#' @export
item.ttlcor = function(M, is.corr=FALSE){
  if (is.corr){
    # M is a correlation matrix
    covttl = sum(M)
    covcolsums = colSums(M)
    (covcolsums-1) / sqrt(covttl - 2*covcolsums + 1)   # +1 because subtract own var=1 twice
  }
  else {
    # use means to calculate item-total correlation
    ttl = rowSums(M)
    out = rep(NA, ncol(M))
    for (ii in 1:ncol(M)){
        out[ii] = cor(M[,ii], ttl - M[,ii])
    }
    return(out)
  }
}

#' Calculate inter-item correlation (internal consistency).
#' Note that this returns all of the relevant correlations (then you can average them).
#' 
#' @param M data.frame or matrix
#' @param M2 optional data.frame or matrix to cross-correlate
#' @return array of all inter-item correlations
#' @export

item.ic = function(M, M2=NULL){
  # inter-item correlation (average off-diag corr)
  if (is.null(M2)){
    cors = cor(M)
    return( cors[lower.tri(cors)] )
  }
  else {
    return(c(cor(M, M2)))
  }
}

