################Copy number genetic alteration profiles and gene expression profiles #######################
#####@cna_data copy number calls in each sample
#####@exp_data PCGs expression profile
#####@thre minimum number of samples with genetic alterations
deal_data <- function(cna_data,exp_data,thre){
result <- list()
exp_data <- exp_data[!duplicated(exp_data$Hugo_Symbol),] 
cna_data <- cna_data[!duplicated(cna_data$Hugo_Symbol),]
exp_data <- exp_data[complete.cases(exp_data),]
instect_sample <- intersect(colnames(cna_data)[-c(1,2)], colnames(exp_data)[-c(1,2)])
row.names(exp_data) <- exp_data$Hugo_Symbol
exp_data <- exp_data[,instect_sample] 
exp_result <- gene_deal(exp_data)
instect_gene <- cna_data$Hugo_Symbol[which(cna_data$Hugo_Symbol %in% row.names(exp_result))]
row.names(cna_data) <- cna_data$Hugo_Symbol
cna_data <- as.matrix(cna_data[instect_gene,instect_sample])
cna_data[which(cna_data==1)] <- 0
cna_data[which(cna_data==-1)] <- 0    
cna_result <- cna_deal(exp_result,cna_data,thre)
result[[1]] <- log2(exp_data+1)  
result[[2]] <- cna_result
return(result)	                                        
}
