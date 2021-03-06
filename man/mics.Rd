% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mics.R
\name{mics}
\alias{mics}
\title{Testing cell-type-specific mediation effects from bulk methylation data}
\usage{
mics(
  meth_data,
  S,
  X,
  Y,
  cell_prop = NA,
  meth_ref = NA,
  MCP.type = "FDR",
  cpg_ind = NA,
  thd = 0,
  maxp_sq = FALSE
)
}
\arguments{
\item{meth_data}{The observed methylation matrix. Rows represent CpG sites and columns are samples.
Its CpG site names are required in order to automatically align with the CpG site names of the methylation reference when estimating cellular compositions.}

\item{S}{The exposure vector. Each component corresponds to a sample.}

\item{X}{The confounding factor matrix. Rows are samples and columns represent confounding factors.
Note that if there is only one confounding factor, X also needs to be in a matrix form and has only one column.}

\item{Y}{The outcome vector. Each component corresponds to a sample.}

\item{cell_prop}{The cellular proportion matrix, where rows correspond to cell types and columns are samples. The default is NA, in which case the cell proportions are estimated from the methylation reference.}

\item{meth_ref}{The methylation reference matrix. If cell_prop is NA, meth_ref must be specified, where rows represent CpG sites and columns are cell type annotations.
Its CpG site names are required in order to automatically align with the CpG site names of the observed methylation matrix when estimating cellular compositions.}

\item{MCP.type}{The type of multiple comparison procedure when implementing Heller/Sampson method. It is either FWER or FDR. The default is FDR.}

\item{cpg_ind}{The index vector of the informative cpg sites. The default is NA, which means that all cpg sites are used to estimate
the cellular compositions. Otherwise, only the informative cpg sites are used.}

\item{thd}{The threshold used to filter out cell types with a small proportion. After estimating the cellular composition for each sample,
we filter out the cell types with its proportions' 0.75 quantile across all samples less than the threshold. The default is 0, which means that no cell type would be removed.}

\item{maxp_sq}{A Boolean type argument indicating whether outputs combined p-values using maximum and squaring procedure. The default is FALSE.}
}
\value{
A list containing the following items. \item{pval_mediator}{The cell-type-specific p-value matrix by regressing the mediator to the exposure, where rows represent CpG sites and columns are cell types.}
\item{pval_outcome}{The cell-type-specific p-value matrix by regressing the outcome to the mediator and the exposure, where rows represent CpG sites and columns are cell types.}
\item{pval_joint_MultiMed}{The joint cell-type-specific p-value matrix using the Heller/Sampson procedure, where rows represent CpG sites and columns are cell types.}
\item{pval_joint_sq}{Only when the argument maxp_sq is TRUE, the returned value is the joint cell-type-specific p-value matrix by maximizing and squaring, where rows represent CpG sites and columns are cell types. Otherwise, it is NA.}
\item{P_matr}{The cell proportion matrix, where rows are cell types and columns correspond to samples.}
}
\description{
Testing cell-type-specific mediation effects from bulk methylation data
}
\examples{
rm(list = ls())
data(example_data)
out <- mics(meth_data = Ometh, S = S, X = X, Y = Y, cell_prop = P_matr, 
                       MCP.type = "FDR", maxp_sq = TRUE)

pval_MultiMed <- out$pval_joint_MultiMed

fdr_thred <- 0.2

ind1 <- which(pval_MultiMed[,1] < fdr_thred)

ind2 <- which(pval_MultiMed[,2] < fdr_thred)

ind3 <- which(pval_MultiMed[,3] < fdr_thred)

#detected CpG sites in cell type 1
ind1

#detected CpG sites in cell type 2
ind2

#detected CpG sites in cell type 3
ind3
}
