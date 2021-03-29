# constructing hallmark dysregulation profile ---------------------------------------------------------
#exp_profile : expression profile
##m_hallmark:hallmark geneset list
dys_hallmark_profile <- function(exp_profile,m_hallmark){
  EXP <- exp_profile/rowMeans(exp_profile)
  hallmark_fc <- gsva(as.matrix(EXP),method="ssgsea",m_hallmark,verbose=FALSE) 
  last <- c()
  result_fc <- list()
  for(i in 1:1000){###1000
    s <- sample(row.names(exp_profile),length(row.names(exp_profile)))
    row.names(EXP) <- s
    result_fc[[i]] <- gsva(as.matrix(EXP),method="ssgsea",m_hallmark,verbose=FALSE)
  }
  for (j in names(m_hallmark)) {
    p <- c()
    m_mat <- c()
    for (z in 1:1000) {
      m_mat <- rbind(m_mat,result_fc[[z]][j,]) 
    }
    for (k in 1:ncol(m_mat)) {
      if( hallmark_fc[j,k] > median(m_mat[,k])){####activation
        z <- length(which(m_mat[,k] > hallmark_fc[j,k]))/1000
        p <- c(p,z)
      }
      if( hallmark_fc[j,k] < median(m_mat[,k])){####inactivation
        z <- length(which(m_mat[,k] <= hallmark_fc[j,k]))/1000
        p <- c(p,z)
      }
    }
    r <- hallmark_fc[j,]
    p.fdr <- p.adjust(p,method = "BH")
    r[which(p.fdr > 0.05)] <- 0
    last <- rbind(last,r)
  }
  row.names(last) <- names(m_hallmark)
  #save(result_fc,file=paste(path,"result_fc.RData",sep=""))
  return(last)
}
