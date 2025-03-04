% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{process_diann_report}
\alias{process_diann_report}
\title{Processs the imported DIA-NN precursor data}
\usage{
process_diann_report(
  data_raw,
  sample_names,
  sample_groups,
  peptide.id = "Stripped.Sequence",
  quant_column = "Precursor.Normalised",
  min_val_per_group = 1,
  norm.method = "none",
  proteotypic.only = F,
  Q.Value = 1,
  PG.Q.Value = 1,
  Lib.Q.Value = 1,
  Lib.PG.Q.Value = 1,
  protein.q = 1,
  gg.q = 1
)
}
\arguments{
\item{data_raw}{}

\item{sample_names}{}

\item{sample_groups}{}

\item{peptide.id}{}

\item{quant_column}{}

\item{min_val_per_group}{}

\item{norm.method}{}

\item{proteotypic.only}{}

\item{Q.Value}{}

\item{PG.Q.Value}{}

\item{Lib.Q.Value}{}

\item{Lib.PG.Q.Value}{}

\item{protein.q}{}

\item{gg.q}{}
}
\description{
Processs the imported DIA-NN precursor data
}
