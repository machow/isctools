% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/estimate_cfa_isc.R
\name{isc.model}
\alias{isc.model}
\title{Return lavaan model string with factor loadings derived from each input}
\usage{
isc.model(...)
}
\arguments{
\item{...}{named argument with character array.}
}
\description{
Return lavaan model string with factor loadings derived from each input
}
\examples{
isc.model(f1 = c('a','b'),   # f1 =~ a + b
          f2 = c('x','z'))   # f2 =~ x + z
}

