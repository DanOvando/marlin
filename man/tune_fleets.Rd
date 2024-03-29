% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tune_fleets.R
\name{tune_fleets}
\alias{tune_fleets}
\title{tune_fleets tunes parameters of the fleet model to achieve desired initial conditions.
note that this is not exact: post-tuning values will not perfectly match inputs since
for example some tuning steps depend on prior tuning step, making it difficult to tune everything
at once.}
\usage{
tune_fleets(
  fauna,
  fleets,
  years = 50,
  tune_type = "explt",
  tune_costs = TRUE,
  fine_tune_costs = TRUE
)
}
\arguments{
\item{fauna}{a fauna object}

\item{fleets}{a fleet object}

\item{years}{the number of years to tune over}

\item{tune_type}{one of 'explt' or 'depletion' to tune to an exploitation rate or a target depletion (B/B0)}

\item{tune_costs}{TRUE or FALSE to tune costs to a target cost to revenue ratio}

\item{fine_tune_costs}{TRUE or FALSE}
}
\value{
tuned fleet object
}
\description{
tune_fleets tunes parameters of the fleet model to achieve desired initial conditions.
note that this is not exact: post-tuning values will not perfectly match inputs since
for example some tuning steps depend on prior tuning step, making it difficult to tune everything
at once.
}
