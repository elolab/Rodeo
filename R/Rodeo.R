# Input: 'E' is a bulk expression matrix (rows: genes, cols: samples),
#        'C' is a cell type proportion matrix (rows: cell types, cols: samples).
# output: returns a matrix S, which indicates how stongly each cell type (col) expresses each gene (row) when pure. 
Rodeo = function(E, C){
  
  # Check that input is ok
  args = CheckInput(E, C)
  E = args$E
  C = args$C
  
  # Calculate S using robust linear model
  S = EstimateS(E, C)
  
  return(S)
  
}







