getdishis_from_edge<-function(edges) {
  #get cluster size histogram for all frame, and add them up, return the result directly
  foreach(index=1:length(edges),.combine="+") %dopar% (
    getdishis_from_edge_base(edges[[index]])
  )
 
}