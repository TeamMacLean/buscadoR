% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/buscador.R
\name{buscar}
\alias{buscar}
\title{check and combine the search files and perform the protein annotation}
\usage{
buscar(
  deeptmhmm = NULL,
  hmmer = NULL,
  blast = NULL,
  fasta = NULL,
  hmmer_eval_cutoff = 20
)
}
\arguments{
\item{deeptmhmm}{results file from deeptmhmm (as .gz or plain text)}

\item{hmmer}{results file from hmmer --domtblout (as .gz or plain text)}

\item{blast}{results file from blastp -outfmt 6 (as. gz or plain text)}

\item{fasta}{original fasta file used in searches}

\item{hmmer_eval_cutoff}{evalue cutoff for hmmer search}
}
\value{
busc object
}
\description{
check and combine the search files and perform the protein annotation
}
