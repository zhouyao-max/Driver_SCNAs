######To obtain expression profile
######expression in at least 10% samples
######log2(x+1)
gene_deal <- function(exp_mat){
  genes <- row.names(exp_mat)
  n1 <- c();
  num_1 <- apply(exp_mat,1,function(x){n1 <- c(length(which(x==0)),n1)})
  gene_2 <- genes[which(num_1 < 0.9*ncol(exp_mat))]
  exp_matrix <- log2(as.matrix(exp_mat[gene_2,])+1)
  return(exp_matrix)
}