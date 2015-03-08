#' Return lavaan model string with factor from each input
#' @param ... named argument with character array.
#' @export
#' @examples
#' isc.model(f1 = c('a','b'),   # f1 =~ a + b
#'           f2 = c('x','z'))   # f2 =~ x + z
isc.model = function(...){
  inputs = list(...)
  factors = lapply(names(inputs), function(k){
    paste0(k, ' =~ ', paste(inputs[[k]], collapse=' + '))
  })
  
  paste(factors, collapse='\n')
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
  model.1f = isc.model(f1 = names(df))
  
  model.2f = isc.model(f1 = grep(f.names[1], names(df), value=TRUE), 
                       f2 = grep(f.names[2], names(df), value=TRUE)
  )
  
  out = list()
  out$fit = list(f1 = lavaan::cfa(model.1f, df, std.ov=TRUE, std.lv=TRUE),
                 f2 = lavaan::cfa(model.2f, df, std.ov=TRUE, std.lv=TRUE))
  
  out$diff = lavaan::anova(out$fit$f1, out$fit$f2)
  fit.indices = c("chisq", "df", "cfi", "rmsea", "srmr", "mfi", "aic", "bic")
  out$fit.indices = do.call(rbind, lapply(out$fit, lavaan::fitMeasures, fit.indices))
  
  out                                        
}