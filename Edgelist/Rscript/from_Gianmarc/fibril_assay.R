#Function to calculate the yield of different fibril topologies
#
# For documentation on using orca package, see paper: https://www.jstatsoft.org/article/view/v071i10/0 
# 
# Gianmarc Grazioli, Ph.D.
#

fibril_assay<-function(net){
  require(orca)
  edgelist<-as.edgelist(net)
  #Get the orbit degrees that we'll need for classification
  od<-as.data.frame(count5(edgelist))
  #  test for 2-ribbon:
  isNode.in.2ribbon<-function(odRow){
    rowTestResult <-(((odRow["O0"]==2)&(odRow["O2"]==1)&(odRow["O8"]==1))|((odRow["O0"]==3)&(odRow["O7"]==1)&(odRow["O8"]==2)))
    rowTestResult
  }
  #  test for 1,2 2-ribbon:
  isNode.in.12.2ribbon<-function(odRow){
    rowTestResult <-(((odRow["O12"]==1) | (odRow["O13"]==1)) | (odRow["O61"]==1))
    rowTestResult
  }
  membership <-as.data.frame(matrix(0, nrow=net$gal$n,ncol=2))
  colnames(membership) <- c("2-ribbon", "1,2 2-ribbon")
  for (i in 1:nrow(od)) {
    if(isNode.in.2ribbon(od[i,])) membership[i, "2-ribbon"] = 1
    if(isNode.in.12.2ribbon(od[i,])) membership[i, "1,2 2-ribbon"] = 1
  }
  return(membership)
}

