# plrs»Ø¹é ------------------------------------------------------------------
###plrs_find_factor is a program to select the important regulatory factors foractivity of a hallmark
##binary_net candidate gene-hallmark network
###pthr is the p threshold
##path_output:output results storage path
##hall_mat_one_zero: hallmark dysregulation profile
###cna_profile: cna profile
plrs_find_factor <- function(binary_net,cna_profile,hall_mat_one_zero,pthr,path_output){
  last_binary_net <- c()
  hallmarks <- unique(binary_net$hallmark)
  for(i in hallmarks){
    genes <- binary_net[which(binary_net$hallmark %in% i),"gene"]
    if(length(genes)==1){ 
      fit.data <- as.data.frame(t(rbind(cna_profile[genes,],hall_mat_one_zero[i,])))
      colnames(fit.data) <- c(genes,"y")
      lm.D9 <- lm(as.formula(paste("y ~ ",genes)),data=fit.data) 
      variable_coef <- summary(lm.D9)$coefficients[2,1]
      pvalue <- summary(lm.D9)$coefficients[2,4]
      gene_hall <- dplyr::filter(binary_net,gene %in% genes & hallmark == i) 
      gene_hall <- cbind(cbind(gene_hall,variable_coef),pvalue)
      last_binary_net <- rbind(last_binary_net,gene_hall)
    }else{
      fit.data <- as.data.frame(t(rbind(cna_profile[genes,],hall_mat_one_zero[i,])))
      num_d <- dim(fit.data)[2]
      colnames(fit.data) <- c(genes,"y")
      x <- as.matrix(fit.data[,genes]) 
      y <- as.matrix(fit.data[,"y"])
      mod<- plsr(y ~ x, data = fit.data, validation = "LOO", jackknife = TRUE,scale=T)
      ncomp<-mod$ncomp ##extract the all ncomps
      ##calculate the adjcv for all ncomps, adjcv is a measure to evaluated the regression, the smaller the better
      adjcv<-(RMSEP(mod)$val)[2,1,2:ncomp]; ##find the number of coponents for plsr with the minmum RMSEP
      ##relapce the adjcv with the max adjcv where the real adjcv is NaN
      pos<-which(adjcv=="NaN")
      if(length(pos)){
        adjcv[pos] <- max(adjcv[-pos])
      }
      ## find the number of ncomp respoding the minmum adjcv
      min_adjcv <- min(adjcv)
      ncomp <- which(adjcv==min_adjcv)
      ##use the jack test for calculating the p value for coeffients 
      jtest<-jack.test(mod,ncomp = ncomp);
      variable_coef<-(jtest$coefficients)[1:(num_d-1),1,1];
      pvalue<-(jtest$pvalue)[1:(num_d-1),1,1];
      #pos<- which(pvalue<=pthr)
      #if(length(pos)){
      gene <- genes
      gene_hall <- dplyr::filter(binary_net,gene %in% gene & hallmark == i) 
      gene_hall <- cbind(cbind(gene_hall,variable_coef),pvalue)
      last_binary_net <- rbind(last_binary_net,gene_hall)
      #}
      
    }
  }
  last_binary_net <- filter(last_binary_net,pvalue < pthr)
  save(last_binary_net,file = paste(path_output,"last_net.RData",sep = ""))
  return(last_binary_net)
}
