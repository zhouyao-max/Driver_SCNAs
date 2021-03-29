#########RWR#####################
####weight_ppinet weighted PPI
####gamma=0.3 
####i seed nodes
rw_weight<-function(ppi_net,gamma,gene=i) {
  p0 <- rep(0,length(colnames(ppi_net)))
  names(p0) <- colnames(ppi_net)
  p0[gene] <- 1 
  p0<- as.matrix(p0)
  PT <- p0 
  delta <- 1
  ppi_net<-as.matrix(ppi_net)
  while(delta > 1e-8) {
    PT1 <- ((1-gamma)*ppi_net) %*% PT
    PT2 <- (gamma*p0) 
    PT3 <- PT1 + PT2
    delta <- sum(abs(PT3 - PT))
    PT <- PT3
  }
  return(PT)
}

