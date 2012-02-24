\name{sim_data_probit}
\alias{sim_data_probit}
\title{Simulate temporal occupancy data}
\usage{
  sim_data_probit(n.indiv, n.obs, q.window = 2, tau, P)
}
\arguments{
  \item{n.indiv}{integer giving number of sites
  (individuals)}

  \item{n.obs}{integer giving number of observation per
  site (assumed equal for every individual)}

  \item{q.window}{is an integer giving the number of
  replicates for each eta (must be >1; default is 2)}

  \item{tau}{gives the precision (1/variance) of the CAR
  process}

  \item{P}{is a vector with 4 detection probabilities: (1)
  state_t=not occupied, detected at t-1, (2) state_t =
  occupied, detected at t-1, (3) state_t=not occupied, not
  detected at t-1, and (4) state_t = occupied, not detected
  at t-1}
}
\value{
  Stacked data.frame object with "Response" as the first
  column (0=missing, 1=Not occupied, 2=Occupied), "Hr"
  (timestep) as second column, and "Ind" (individual
  identifier) as third column
}
\description{
  Simulate temporal occupancy data
}
\author{
  Paul Conn
}
\keyword{model,}
\keyword{occupancy,}
\keyword{probit}
\keyword{simulation}

