########Building the candidate gene-hallmark Network using random walk with restart(RWR)################
####cna_profile:copy number calls in each sample
####exp_profile: PCGs expression profile
####STRING: PPI data
####path_output:output results storage path
####H1:hallmark genesets
gene_hallmark_candidata_net <- function(cna_profile,exp_profile,STRING,H1){
  result_last <- c()
  STRING <-  STRING[which(STRING$protein1 %in% row.names(exp_profile)),] 
  STRING <-  STRING[which(STRING$protein2 %in% row.names(exp_profile)),]
  gene <- unique(c(STRING$protein1,STRING$protein2)) 
  CG <- intersect(unique(c(STRING$protein1,STRING$protein2)),row.names(cna_profile))
  rw_mat <- c() 
  weight_ppinet <- weight_ppi(exp_profile,STRING) 
  weight_ppinet <- weight_ppinet[complete.cases(weight_ppinet),] 
  ######Building the adjacency matrix
  weight_ppinet$protein1 <- factor(weight_ppinet$protein1,levels=gene)
  weight_ppinet$protein2 <- factor(weight_ppinet$protein2,levels=gene)
  ppi_net <- spread(weight_ppinet, key = protein1, value = cor_score, fill = 0,drop=F)
  row.names(ppi_net) <- ppi_net$protein2
  ppi_net <- ppi_net[,-1]
  ppi_net <- ppi_net + t(ppi_net)
  ppi_net <- abs(ppi_net)
  ppi_net <- as.matrix(sweep(ppi_net, 2, colSums(ppi_net), "/")) 
  ppi_net[which(is.nan(ppi_net),arr.ind = T)] <- 0  
  #######RWR##################
  for(i in CG ){
    score<-rw_weight(ppi_net = as.data.frame(ppi_net),gamma = 0.3,gene=i) 
    rw_mat <- cbind(rw_mat,score[gene,])
  }
  colnames(rw_mat) <- CG
  #save(rw_mat,file=paste(path,"rw_mat.Rdata")) 
  ########GSEA########################
  name <- unique(H1$Hallmark)
  id <- c()
  for (i in 1:ncol(rw_mat)){
    x <- rw_mat[,i]
    names(x) <- row.names(rw_mat)
    geneList_human <- sort(x,decreasing=T) 
    n <- which(geneList_human==0)[1]
    if( n < 500){
      id <- c(id,i)
      next}
    gsea_h1_human <- as.data.frame(GSEA(geneList_human, TERM2GENE = H1, verbose=FALSE, pvalueCutoff = 1)) 
    gsea_h1_human <- gsea_h1_human[name,]
    hall_sample <- rep(0,len=50)
    names(hall_sample) <- name 
    pos1 <- which(gsea_h1_human$NES > 0 & gsea_h1_human$qvalues < 0.05) #####NES>0，pvalue < 0.05
    hall_sample[pos1] <- 1
    result_last <- cbind(result_last,hall_sample)   
  }
  colnames(result_last) <- CG[-id]
  result_last <- result_last[,which(colSums(result_last) > 0)]
  return(result_last)
}