# Input: 'E' is a bulk expression matrix (rows: genes, cols: samples),
#        'C' is a cell type proportion matrix (rows: cell types, cols: samples)
# Output: This function checks that all inputs are valid and returns them in controlled format
CheckInput = function(E, C){
  
  # Check that 'E' and 'C' are some sort of tables and convert them into matrices
  if(!any( class(E) %in% c("matrix", "data.frame") )) stop("'E' should be a matrix or a data frame.")
  if(!any( class(C) %in% c("matrix", "data.frame") )) stop("'C' should be a matrix or a data frame.")
  E = as.matrix(E)
  C = as.matrix(C)
  
  # 'E' and 'C' should have the same number of samples
  if(ncol(C) != ncol(E)) stop("Bulk matrix 'E' and cell type proportion matrix 'C' should have equal number of columns (i.e. samples).")
  
  # Remove genes that are unexressed in at least 40% of samples (in both sample groups)
  # if(is.na(groups)){
  #   limit = floor(0.4*ncol(E))
  #   dropgenes = which(apply(E, 1, function(g){sum(g==0)}) >= limit)
  #   #if(length(dropgenes) > 0) E = E[-dropgenes,]
  # } else{
  #   limit1 = floor(0.4*sum(groups==1))
  #   limit2 = floor(0.4*sum(groups==2))
  #   dropgenes1 = which(apply(E[,groups==1], 1, function(g){sum(g==0)}) >= limit1)
  #   dropgenes2 = which(apply(E[,groups==2], 1, function(g){sum(g==0)}) >= limit2)
  #   dropgenes = intersect(dropgenes1, dropgenes2)
  #   if(length(dropgenes) > 0) E = E[-dropgenes,]
  # }
  
  
  return(list(E=E, C=C))
}

