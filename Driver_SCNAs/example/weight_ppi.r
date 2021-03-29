###############building a weighted PPI network######################################################
##exp_profile:expression profile
##STRING:PPI data
weight_ppi <- function(exp_profile,STRING){
  x <- exp_profile[match(STRING$protein1,row.names(exp_profile)),]
  y <- exp_profile[match(STRING$protein2,row.names(exp_profile)),]
  d <- as.matrix(cbind(x,y))
  cor_score <- apply(d,1,function(x){cor(x[1:(dim(d)[2]/2)],x[(dim(d)[2]/2 +1):length(x)],method = "pearson")})
  result <- cbind(STRING ,cor_score)
  return(result)
}