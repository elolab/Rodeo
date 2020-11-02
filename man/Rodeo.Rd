\name{Rodeo}
\alias{Rodeo}
\title{
RObust DEcOnvolution
}

\description{
This function estimates how stongly each cell type expresses each gene based on bulk expression data and cell type proportions. The estimates are calculated using robust linear regression and they benefit from large sample size.
}

\usage{
Rodeo(E, C)
}

\arguments{
  \item{E}{ a bulk expression matrix (rows: genes, cols: samples) }
  \item{C}{ a cell type proportion matrix (rows: cell types, cols: samples) }
}

\value{ A matrix S, which indicates how stongly each cell type (col) expresses each gene (row) when pure. The cell types are those present in input matrix C (cell type proportions) and genes are those from input matrix E (bulk expression).}

\references{
%% ~put references to the literature/web site here ~
}

\author{
Maria K. Jaakkola and Laura L. Elo
}

\note{
Using this function requires package MASS
}

\examples{

# library(MASS)

# Generate artificial cell type proportion matrix C for 10 samples (cols) and 3 cell types (rows) 
C = matrix(runif(30, min=0, max=100), nrow=3)
C = t(t(C)/colSums(C))
colnames(C) = paste("sample", 1:10, sep="")
rownames(C) = paste("celltype", 1:3, sep="")

# Generate artifical bulk matrix E (this is measured expression data in real applications)
Sknown = matrix(runif(9, min=0, max=100), nrow=3)
rownames(Sknown) = paste("gene", 1:3, sep="")
E = Sknown \%*\% C
rownames(E) = rownames(Sknown)
colnames(E) = colnames(C)

# Run Rodeo to estimate how strongly each of the 3 cell types expresses the 3 genes
# Sestimated = Rodeo(E, C)

## The function is currently defined as
function (E, C) 
{
    args = CheckInput(E, C)
    E = args$E
    C = args$C
    S = EstimateS(E, C)
    return(S)
  }
}

