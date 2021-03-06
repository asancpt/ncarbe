% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/beNCA.R
\name{beNCA}
\alias{beNCA}
\title{Performs a statistical analysis of 2x2 bioequivalence study}
\usage{
beNCA(concData, SUBJ = "SUBJ", GRP = "GRP", PRD = "PRD", TRT = "TRT",
  method = "kbe", ...)
}
\arguments{
\item{SUBJ}{Subject ID, any data type}

\item{GRP}{column name in which information of "RT" or "TR" exists.}

\item{PRD}{column name in which information of 1 or 2 exists.}

\item{TRT}{column name in which information of "R" or "T" exists.}

\item{method}{\code{kbe} by authors or \code{nlme} package uploaded on CRAN}
}
\value{
returns text results of statistical analysis of 2x2 bioequivalence study including and \code{beNCAdataset.csv} file in the working directory which can be run in SAS.
}
\description{
\code{beNCA} returns text output of a statistical analysis of 2x2 bioequivalence study. Analysis of variance, least square means, confidence interval, and sample size will be calculated for AUClast and Cmax. Wilcoxon Signed-Rank Test and Hodges-Lehmann Estimate will be performed for Tmax.
}
\details{
Reference: Shein-Chung Chow, Jen-pei Liu. Design and Analysis of Bioavailability and Bioequivalence Studies, 3rd ed. 2008. (ISBN:9781584886686)
}
\examples{
file <- system.file('example', 'beConc.csv', package = 'ncarbe')
concData <- read.csv(file, as.is = TRUE)
beNCA(concData)
}
