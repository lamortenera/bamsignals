\name{coverage}
\alias{coverage}
\title{Compute read coverage from a bam file.}
\usage{
  coverage(gr, bampath, mapqual = 0, format = T)
}
\arguments{
  \item{gr}{GenomicRanges object used to specify the
  regions}

  \item{bampath}{path to the bam file storing the read. The
  file must be indexed.}

  \item{mapqual}{discard reads with mapping quality
  strictly lower than this parameter. The value 0 ensures
  that no read will be discarded, the value 254 that only
  reads with the highest possible mapping quality will be
  considered.}

  \item{format}{if the ranges have all the same width this
  method will return a matrix.}
}
\value{
  a list with the following arguments: \item{counts}{the
  vector containing the read counts. This will be formatted
  into a matrix depending on whether the ranges have all
  the same length.} \item{starts, ends}{Vectors defining
  the boundaries of the count vector.  To extract counts
  relative to the i-th range, use
  \code{as.numeric(counts)[starts[i]:ends[i]]}, or the
  \code{getSignal} function to preserve the formatting.}
}
\description{
  Compute read coverage in the regions specified by a
  GenomicRanges object.
}

