% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/classic_isc.R
\name{item.rel}
\alias{item.rel}
\title{Calculate item reliabilty correlation from full data.
Note that this calculates the correlation of each column against the others.}
\usage{
item.rel(M, is.corr = FALSE)
}
\arguments{
\item{M}{data.frame, matrix, or lavaan fit object.
if a data.frame or matrix is given, calculate item-total correlation
(using one against the others). If a lavaan fit object is given, return
a data.frame with mean factor loadings by factor.}

\item{is.corr}{boolean indicating if M is a corraltion matrix}
}
\value{
array with item-total correlation for each column
}
\description{
Calculate item reliabilty correlation from full data.
Note that this calculates the correlation of each column against the others.
}

