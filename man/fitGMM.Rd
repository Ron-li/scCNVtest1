% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/supp_functions_1.4.R
\name{fitGMM}
\alias{fitGMM}
\title{Fit Gaussian mixture model}
\usage{
fitGMM(fc, SNN_results, GMMpara, null_dist, project)
}
\arguments{
\item{fc}{Fold change matrix}

\item{SNN_results}{A dataframe with group information for every cell}

\item{GMMpara}{Pretrained parameters for Gaussian mixture model}

\item{null_dist}{Pre-calculated statistics distribution of a normal cohort}

\item{project}{names for the file to be saved}
}
\value{
creates a dataframe with cell names, neutral proportion, neutral probability and p-value
}
\description{
This function calculates some statistics for each of the cell to help identify the tumor and normal group.
}
