\name{getSignal}
\alias{getSignal}
\title{Handle the output of the \code{pileup} and \code{coverage} functions.}
\usage{
getSignal(pu, n)
}
\arguments{
  \item{pu}{output of a \code{pileup} function call}

  \item{n}{index of the desired range in the GenomicRanges
  object used in the pileup function call.}
}
\description{
Extract the signal corresponding to a specified range
maintaining the correct formatting.
}

