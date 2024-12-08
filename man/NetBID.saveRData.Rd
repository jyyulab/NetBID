% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pipeline_functions.R
\name{NetBID.saveRData}
\alias{NetBID.saveRData}
\title{Save Data Produced by Corresponding NetBID2 Pipeline Step.}
\usage{
NetBID.saveRData(network.par = NULL, analysis.par = NULL, step = "exp-load")
}
\arguments{
\item{network.par}{list, stores all related datasets from network construction pipeline step.}

\item{analysis.par}{list, stores all related datasets from driver analysis pipeline step.}

\item{step}{character, name of the pipeline step decided by user for easier reference.}
}
\description{
\code{NetBID.saveRData} is a function to save complicated list object generated by certain steps of NetBID2's pipeline
(e.g. load gene expression file from GEO, 'exp-load').
This function is not essential, but it is highly suggested for easier pipeline step checkout and reference.
}
\details{
There are two important steps in the NetBID2 pipeline, network construction and driver analysis.
User can save these two complicated list objects, network.par and analysis.par.
Assigning the \code{step} name to save the RData for easier reference.
Calling \code{NetBID.loadRData} to load the corresponding step RData, users can avoid repeating the former steps.
}
\examples{
\dontrun{
analysis.par <- list()
analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
NetBID.saveRData(analysis.par=analysis.par,step='ms-tab_test')
}
}
