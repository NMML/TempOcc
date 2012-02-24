\name{plot_temp_occ}
\alias{plot_temp_occ}
\title{Plot temporal occupancy data for inidivual "sites" using ggplot2}
\usage{
  plot_temp_occ(Dat, Time, Ind, plot.ncol = 1, DayMoYear,
    Label.times, First, Last, min.bar.length = 1)
}
\arguments{
  \item{Dat}{matrix object with number of rows equal to the
  number of sites (individuals) and the number of columns
  equal to the number of occasions (analagous to an
  encounter history).  NAs give missing, 0 gives not hauled
  out, 1 gives hauled out}

  \item{Time}{a vector giving times that each occasion is
  associated with; for evenly spaced occasions this could
  simply be c(1:ncol(Dat))}

  \item{Ind}{a vector of individual IDs (can be a character
  or integer vector)}

  \item{plot.ncol}{an integer specifying the number of
  columns desired in the figure}

  \item{DayMoYear}{a vector giving the day, month, and
  year, of each element in Time (formatted using the
  date.mmddyy function in the 'date' package)}

  \item{Label.times}{a vector giving the elements of
  DayMoYear to plot as labels on the graph}

  \item{First}{an integer vector giving the indices of the
  first observations of animals in the study; assumed to be
  in same order as "Ind"}

  \item{Last}{an integer vector givign the indices of the
  last observations of animals}

  \item{min.bar.length}{relative number of hours each
  haulout "block" takes up on the graph}
}
\description{
  Plot temporal occupancy data for inidivual "sites" using
  ggplot2
}
\author{
  Paul Conn
}
\keyword{occupancy}
\keyword{plot,}
\keyword{temporal}

