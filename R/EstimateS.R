# Input: 'signature' is a signature matrix of average expressions in pure cell populations (cols: cell types, rows: genes),
#        'data' is a bulk expression matrix (cols: samples, rows: genes)
# Output: cell type proportion matrix (cols: samples, rows: cell types)
EstimateS = function(data, C){
  
  # Initialize Sraw
  Sraw = mat.or.vec(nrow(data), nrow(C))
  colnames(Sraw) = rownames(C)
  rownames(Sraw) = rownames(data)
  
  # Detect mean cell type expressions for each gene at time (i.e. j = gene)
  for (j in 1:nrow(data)) {
    
    sampledata = as.data.frame(cbind(data[j, ], t(C))) # rows: samples, cols: bulk expression of the current gene + cell types
    colnames(sampledata)[1] = "BulkExpression"
    
    # Fit regression and extract mean cell type expressions
    rlm_model = rlm(BulkExpression ~ 0+., data=sampledata) # , psi=psi.bisquare
    indices = 1:nrow(C)
    coefs = rlm_model$coefficients
    
    # Re-fit without cell types with negative coefficient
    while(any(coefs < 0, na.rm=TRUE)){
      indices = indices[which(coefs > 0)]
      rlm_model = rlm(BulkExpression ~ 0+., data=sampledata[, c(1, indices+1)])
      coefs = rlm_model$coefficients
    }
    
    Sraw[j, indices] = coefs
  }
  
  Sraw[Sraw < 0] = 0
  
  return(Sraw)
}
