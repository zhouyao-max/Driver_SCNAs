######To obtain copy number genetic alteration profile##############################################
###@with a dominant SCNAs type£¨amplification or deletion, p < 0.05, binomial test)
###@CRNA expression coincidence with SCNAs£¨ Wilcox test FDR < 0.05£©
###@SCNAs frequency >=3% 
cna_deal <- function(exp_mat,cna_mat,thre){
  cna_binom <- data.frame() 
  domin_gene <- c()
  for (i in row.names(cna_mat)){
    x <- as.integer(cna_mat[i,])
    x1 <- rep(0,length(x))
    len_loss <- length(which(x==-2)) 
    len_gain <- length(which(x==2)) 
    n <- len_loss+len_gain 
    n1 <- len_loss-len_gain
    if(n ==0 ){next;}
    if(n1>0 & len_gain!=0){ 
      p <- binom.test(len_loss,n,p=0.5,alternative = "greater")
      if (p$p.value < 0.05){
        x1[which(x==-2)]<- -1
        cna_binom <- rbind(cna_binom,x1)
        domin_gene <- c(domin_gene,i)  }}
    if(n1 < 0 & len_loss!=0){ 
      p <- binom.test(len_gain,n,p=0.5,alternative = "greater")
      if (p$p.value < 0.05){
        x1[which(x==2)]<- 1
        cna_binom <- rbind(cna_binom,x1)
        domin_gene <- c(domin_gene,i)
      }  
    }
    if ( len_loss > 0 & len_gain==0){ 
      x1[which(x==-2)] <- -1
      cna_binom <- rbind(cna_binom,x1)
      domin_gene <- c(domin_gene,i)
    }
    if ( len_gain > 0 & len_loss==0){
      x1[which(x==2)] <- 1
      cna_binom <- rbind(cna_binom,x1)
      domin_gene <- c(domin_gene,i)}
  
  }
  row.names(cna_binom) <- domin_gene; 
  colnames(cna_binom) <- colnames(cna_mat)
  exp <- exp_mat[domin_gene,] 
  wilcox_gene <- c()
  p_wilcox<- c()
  all_mut <- c()
  for (i in domin_gene ) {
    cna_sample <- colnames(cna_binom)[which(cna_binom[i,]!=0)] 
    t <- ncol(cna_binom)*thre
    if(length(cna_sample) >=t ){  
      a <- as.numeric(exp[i,cna_sample])
      b <- as.numeric(exp[i,colnames(cna_binom)[which(cna_binom[i,]==0)]])
      if(sum(cna_binom[i,] ) < 0){
        wil <- wilcox.test(a,b,alternative="less")
        p_wilcox <- c(p_wilcox,wil$p.value)
        wilcox_gene <- c(wilcox_gene,i)
      }else{ 
        wil <- wilcox.test(a,b,alternative="g")
        p_wilcox <- c(p_wilcox,wil$p.value)
        wilcox_gene <- c(wilcox_gene,i) 
      }
      
    }
    if(length(cna_sample)== length(colnames(cna_binom))){ 
      all_mut <- rbind(all_mut,cna_binom[i,])}
  }
  p_last <- p.adjust(p_wilcox,"BH") 
  gene_coin <- wilcox_gene[p_last < 0.05]
  cna_last <- cna_binom[gene_coin,]
  cna_last <- rbind(cna_last,all_mut)
  return(cna_last)
} 