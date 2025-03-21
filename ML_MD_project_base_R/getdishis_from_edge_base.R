getdishis_from_edge_base<-function(edges) {
  #do sna analysis
  compdist<-component.dist(edges)
  #return histogram
  table(cut(compdist$csize,breaks=seq(0,300+binwidth,binwidth)))
}