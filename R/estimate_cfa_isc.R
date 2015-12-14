orthogonalize <- function(factors){
  # factors is character vector
  M = combn(factors, 2)       # 2 x n_combinations matrix
  var_char <- apply(M, 2, paste, collapse=' ~~ 0*')
  paste(var_char, collapse='\n')
}

#' Return a list where each entry is an array of item names corresponding to a factor.
#' @export
group_items <- function(i_names, f.names){
  list(f1 = grep(f.names[1], i_names, value=TRUE),
       f2 = grep(f.names[2], i_names, value=TRUE))
}

#' Return lavaan model string with factor loadings derived from each input
#' @param ... named argument with character array.
#' @export
#' @examples
#' isc.model(f1 = c('a','b'),   # f1 =~ a + b
#'           f2 = c('x','z'))   # f2 =~ x + z
isc.model <- function(...){
  inputs = list(...)
  factors = lapply(names(inputs), function(k){
    paste0(k, ' =~ ', paste(inputs[[k]], collapse=' + '))
  })
  
  paste(factors, collapse='\n')
}

#' Create lavaan syntax for various model types
#' @export
build.model <- function(..., type=c('1f', 'hierarchical', 'oblique', 'bifactor')){
  groups = list(...)
  models <- list()
  
  lav.loadings <- isc.model(...)
    
  if ('1f' %in% type){
    models$f1 <- isc.model(f1 = unlist(groups))
  }
  
  if ('hierarchical' %in% type){
    lvl2 <- isc.model(g=names(groups))      # make first level factors load on g
    models$hierarchical <- paste(lav.loadings, lvl2, sep="\n")
  }
  
  if ('oblique' %in% type){
    models$oblique <- lav.loadings
  }
  
  # prepare for CCA models
  if ('bifactor' %in% type){
    lav.orthog   <- orthogonalize(c('g', names(groups)))
    lav.g <- isc.model(g = unlist(groups))
    models$bifactor <- paste(lav.loadings, lav.g, lav.orthog, sep='\n')
  }
  return(models)
}

#' Fit confirmatory factor analysis model with two groups
#'
#' @param df data.frame with names indicating group membership
#' @param f.names a 2-item vector used to grep the 1st and 2nd factor columns
#' @export
#' @examples
#' d = data.frame(gen.cordata(.5, 6, 100))
#' names(d) <- paste0(c(rep(c('f1', 'f2'), each=3)), 1:6)
#' fit.isc(d)
fit.isc = function(df, f.names = c('f1', 'f2')){
  groups = list(f1 = grep(f.names[1], names(df), value=TRUE),
                f2 = grep(f.names[2], names(df), value=TRUE))
  
  model.1f = isc.model(f1 = unlist(groups))
  model.2f = do.call(isc.model, groups)
  
  out = list()
  out$fit = list(f1 = lavaan::cfa(model.1f, df, std.ov=TRUE, std.lv=TRUE),
                 f2 = lavaan::cfa(model.2f, df, std.ov=TRUE, std.lv=TRUE))
  
  out$diff = lavaan::anova(out$fit$f1, out$fit$f2)
  fit.indices = c("chisq", "df", "cfi", "rmsea", "srmr", "mfi", "aic", "bic")
  out$fit.indices = do.call(rbind, lapply(out$fit, lavaan::fitMeasures, fit.indices))
  
  out                                        
}