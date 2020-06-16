% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mutcorcell.R
\name{gene2cellsummary}
\alias{gene2cellsummary}
\title{gene2cellsummary}
\usage{
gene2cellsummary(gene, method = "xCell", mutcell)
}
\arguments{
\item{gene}{Somatic mutant gene name}

\item{method}{Method must be one of "xCell","ssGSEA" and "CIBERSORT".}

\item{mutcell}{The result of `mutcorcell` funtion.}
}
\value{
A matrix shows the short name, full name, pvalue, fdr of the cells driven by a somatic mutation
}
\description{
Function `gene2cellsummary` is a generic function used to produce result summaries of the immune cells driven by a somatic mutation.
}
\examples{
mutcell<-GetExampleData("mutcell") # The result of `mutcorcell` funtion.
genecellsummary<-gene2cellsummary(gene="TP53",mutcell=mutcell)
}