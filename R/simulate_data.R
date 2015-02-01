#' Simulate participant data from latent timecourse.
#' 
#' @param N_part number of participants
#' @param N_time number of timepoints for each participants
#' @param var_rat signal variance to noise ratio
#' @param optional latent timecourse to draw participants from
#' @return list with two entries. latent is the latent timecourse. measured are the outputs.
#' @export
gen.meas = function(N_part, N_time, var_rat, latent=NULL){
  if (is.null(latent)) latent = rnorm(N_time) #mean 0, sd 1
  var_latent = var(latent)
  is_neg = sign(var_rat)
  var_rat = abs(var_rat)
  if (var_rat == 0) {
    is_neg = 1
    err_var = 0
  }
  else err_var = (var_latent * (1 - var_rat)) / var_rat
  M = matrix(nrow=N_time, ncol=N_part)
  
  for (ii in 1:N_part) {
    M[,ii] = is_neg * latent + rnorm(N_time, sd=sqrt(err_var))
  }
  return(list(latent=latent, measured=M))
}

#' Expected correlation between mean timecourses for two perfectly
#' correlated latent variables.
#'
#' @param eta1 proportion of variance that is true signal
#' @param eta2 same but for second latent variable
#' @param N number of measures taken for each latent variable
#' @export
lat.cor = function(eta1, eta2, N){
  var1 = eta1  +  (1-eta1) / N
  var2 = eta2  +  (1-eta2) / N
  sqrt(eta1*eta2) / sqrt(var1*var2)
}



#' Use spectral decomposition to create data with exact correlation matrix
#'
#' @param lambda loadings (correlation with be lambda^2)
#' @param nsubs number of participants
#' @param npoints number of observations per participant
#' @return data.frame with subs columns and npoints rows
#' @export
gen.cordata = function(lambda, nsubs, npoints){
  if (abs(lambda) == 1) return(matrix(rep(rnorm(npoints), nsubs), ncol=nsubs))
  C = diag(1, ncol=nsubs, nrow=nsubs)
  C[lower.tri(C)] = lambda
  C[upper.tri(C)] = lambda
  
  M = matrix(rnorm(npoints*nsubs), ncol=nsubs)
  scores = scale(princomp(M)$scores)
  
  eig = eigen(C)
  Lam = diag(eig$values)
  S = eig$vectors
  scores %*% sqrt(Lam) %*% t(S)
}


#' Return correlation matrix using latent variable model
#' @param g1 indexes for columns of corrmat-to-between. EG g1 = 1:4, g2 = 5:8
#' @param g2 see g1
#' @param lam1 loadings of g1 columns on first latent variable
#' @param lam2 see lam1
#' @param rho correlation between latent variables
#' @export
gen.corrmat = function(g1, g2, lam1, lam2, rho){
  C = matrix(nrow=length(c(g1,g2)), ncol=length(c(g1,g2)))
  C[g1,g1] = lam1^2
  C[g2,g2] = lam2^2
  C[g1,g2] = rho*lam1*lam2
  C[g2,g1] = rho*lam1*lam2
  diag(C) = 1
  C
}
