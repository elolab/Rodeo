# Input: 'Sk' is a known matrix S (cols: cell types, rows: genes),
#        'Se' is an estimated matrix S (cols: cell types, rows: genes),
#        'toself' is either TRUE or FALSE.
# Output: Returns a matrix telling correlations between known (rows) and 
#         estimated (cols) cell type type specific gene expression profiles
#         (i.e. columns in matrices 'Sk' and 'Se'). If 'todelf' is set to TRUE,
#         only the diagonal of the correlation matrix (i.e. correlations between
#         known and estimated same cell types) is returned.
CalculateCellTypeCorrelations = function(Sk, Se, toself=TRUE){
  
  # Ensure that same genes and cell types are in the same order
  genes = intersect(rownames(Sk), rownames(Se))
  message = paste(length(genes), "common genes", sep=" ")
  print(message)
  celltypes = intersect(colnames(Sk), colnames(Se))
  Sk = Sk[genes, celltypes]
  Se = Se[genes, celltypes]
  
  # Calculate cell type correlations
  cors = cor(Sk, Se)
  toreturn = cors
  if(toself) toreturn = diag(cors)
  
  return(toreturn)
}
