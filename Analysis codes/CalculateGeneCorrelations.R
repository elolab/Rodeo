
# Input: 'Sk' is a known matrix S (cols: cell types, rows: genes),
#        'Se' is an estimated matrix S (cols: cell types, rows: genes),
# Output: Returns a named vector of correlations between gene profiles
#         (i.e. rows in 'Sk' and 'Se').
CalculateGeneCorrelations = function(Sk, Se){
  
  # Ensure that same genes and cell types are in the same order
  genes = intersect(rownames(Sk), rownames(Se))
  message = paste(length(genes), "common genes", sep=" ")
  print(message)
  celltypes = intersect(colnames(Sk), colnames(Se))
  Sk = Sk[genes,celltypes]
  Se = Se[genes,celltypes]
  
  # Calculate gene correlations
  cors = NULL
  for(g in 1:nrow(Sk)) cors = c(cors, cor(as.numeric(Sk[g,]), as.numeric(Se[g,])))
  names(cors) = genes
  
  return(cors)
}