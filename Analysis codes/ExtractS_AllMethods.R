# Input: 'E' is a bulk expression matrix (rows: genes, cols: samples),
#        'C' is a cell type proportion matrix (rows: cell types, cols: samples in the same order than in 'E'),
#        'groups' (only for LRCDE) is a vector of 1's (cases) and 2's (controls) indicating sample groups.
# Output: All functions return a matrix S with columns as cell types and rows as genes.

library(Rodeo) # Version 1.0
GetS_Rodeo = function(E, C){
  res = DetectS(E, C, groups=NA, iter=500, seed=252525)
  S = res$S
  return(S)
}

library(CellMix) # Version 1.6.2
GetS_Deconf = function(E, C){
  res = ged(object=E, x=nrow(C), 'deconf')
  S = basis(res)
  return(S)
}

GetS_csQProg = function(E, C){
  res = ged(E, C, 'cs-qprog')
  S = basis(res)
  return(S)
}

GetS_csLsFit = function(E, C){
  res = ged(E, C, 'cs-lsfit')
  S = basis(res)
  return(S)
}

library(csSAM) # Version 1.2.4
GetS_csSAM = function(E, C){
  
  res = csfit(cc=t(C), G=t(E))
  S = t(res$ghat)
  rownames(S) = rownames(E)
  return(S)
}

library(lrcde) # Version 1.0.1.0000
# NOTE: LRCDE recommends log transformed data (skip for SimBulk, use for others)
GetS_LRCDE = function(E, C, groups){
  
  # Define input
  E = t(as.matrix(E))
  C = t(as.matrix(C))
  
  # Make an empty template for estimated S
  Stemplate = mat.or.vec(ncol(E), ncol(C))
  colnames(Stemplate) = colnames(C)
  rownames(Stemplate) = colnames(E)
  
  # Run LRCDE and break results into cell type data frames
  res = lrcde(E, C, groups, "RemoveMe.csv")
  ctres = split(res[[1]], res[[1]]$cell)
  
  # Extract S (general S is a weighted mean of case and control S, weights based on sample sizes)
  S_case = Stemplate
  for(ct in colnames(S_case)){ S_case[, ct] = ctres[[ct]][, "case"]}
  S_control = Stemplate
  for(ct in colnames(S_control)){ S_control[, ct] = ctres[[ct]][, "base"]}
  S_comb = (sum(groups==2)/length(groups) * S_case) + (sum(groups==1)/length(groups) * S_control)
  
  return(S_comb)
}

library(linseed) # Version 0.99.2
# NOTE: C is only used for the number of cell types, the method doesn't really need it
GetS_LinSeed = function(E, C){
  
  #
  lo = LinseedObject$new(E, samples=1:ncol(E), topGenes=nrow(E))
  lo$calculatePairwiseLinearity()
  lo$calculateSpearmanCorrelation()
  #lo$proportions = C
  lo$setCellTypeNumber(nrow(C))
  lo$project("full")
  lo$smartSearchCorners(dataset="full", error="norm")
  lo$deconvolveByEndpoints()
  lo$selectGenes(nrow(E)/nrow(C))
  
  S = lo$signatures
  S = S[rownames(E)[rownames(E) %in% rownames(S)], ]
  return(S)
}

