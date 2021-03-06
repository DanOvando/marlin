% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create-critter.R
\name{create_critter}
\alias{create_critter}
\title{Create Critter}
\usage{
create_critter(
  common_name = "white seabass",
  scientific_name = NA,
  critter_type = "fish",
  seasonal_habitat = list(),
  season_blocks = list(),
  recruit_habitat = NA,
  seasons = 1,
  rec_form = 1,
  fec_form = "power",
  adult_movement = 0,
  adult_movement_sigma = 2,
  recruit_movement = 0,
  recruit_movement_sigma = 10,
  fished_depletion = 0.4,
  init_explt = 0.1,
  explt_type = "f",
  burn_years = 50,
  weight_a = NA,
  ...
)
}
\arguments{
\item{common_name}{}

\item{scientific_name}{}

\item{critter_type}{}

\item{seasonal_habitat}{}

\item{season_blocks}{}

\item{seasons}{}

\item{rec_form}{}

\item{adult_movement}{}

\item{adult_movement_sigma}{}

\item{fished_depletion}{}

\item{init_explt}{initial annual exploitation rate}

\item{explt_type}{f or fmsy}

\item{...}{}

\item{rec_habitat}{}
}
\description{
Creates a critter object. If only a scientific name is provided
create_critter will try and look up relevant life history from
FishLife
}
\details{
Critical inputs are adult_movement, adult_movement_sigma, and resolution
}
