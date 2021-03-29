##########################################################
####STRING:PPI data
####H1:hallmark genesets
##m_hallmark:hallmark geneset list
####path_input：input data storage path
#####path_output:output results storage path
last_programe <- function(STRING=STRING,H1=H1,path_input,path_output,m_hallmark){
  cna_data <-  read.table(paste(path_input,"data_CNA.txt",sep = ""),sep = "\t",header = T,stringsAsFactors = F)
  exp_data <-  read.table(paste(path_input,"data_RNA_Seq_v2_expression_median.txt",sep = ""),sep = "\t",header = T,stringsAsFactors = F)
  result <- deal_data(cna_data,exp_data,thre=0.03)
  exp_profile <- result[[1]]
  cna_profile <- result[[2]] 
  candidata_gene_hallmark <- gene_hallmark_candidata_net(cna_profile,exp_profile,STRING,H1)
  pos <- which(candidata_gene_hallmark==1,arr.ind=T) 
  binary_net <- data.frame(hallmark = row.names(candidata_gene_hallmark)[pos[,1]],gene=colnames(candidata_gene_hallmark)[pos[,2]])
  hall_mat <- dys_hallmark_profile(exp_profile,m_hallmark) #####hallmark dysregulation profile
  last_net <- plrs_find_factor(binary_net=binary_net,cna_profile=cna_profile,hall_mat_one_zero=hall_mat,pthr = 0.05,path_output=path_output)  #### Constructing the active driver gene-hallmark network using partial least squares regression analysis(PLSR)
  return(last_net)
}