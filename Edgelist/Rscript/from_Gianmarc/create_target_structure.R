##################################################################################################
# Function used to create target network structure.
#
# Currently builds all but 3-prism, I'll add that later. 
#
# Gianmarc Grazioli, Ph.D.
##################################################################################################

create.target.structure <- function(n, struct, mode = "path"){
  require(network)
  ### struct: {"1-ribbon", "2-ribbon", "3-prism", "1,2 2-ribbon", "double 1,2 2-ribbon"}
  ### mode: {"path", "cycle"}
  isEven <- (n %% 2 == 0)
  if((struct != "3-prism" & struct != "1-ribbon") & isEven == FALSE){
    print("Warning: use even n value for 2-ribbon or 1,2 2-ribbon")
  }
  isDivBy3 <- (n %% 3 == 0)
  if(struct == "3-prism" & isDivBy3 == FALSE){
    print("Warning: use n that is a multiple of 3 for 3-prism")
  }
  # initialize network object:
  net = network.initialize(n, directed = F)
  # build the desired network:
  # 1-RIBBON
  if(struct == "1-ribbon"){
    for (i in 1:(n-1)) {
      net[i, i+1] = 1
    }
    if(mode == "cycle"){net[1, n] = 1}
  # 2-RIBBON
  } else if(struct == "2-ribbon" | struct == "1,2 2-ribbon" | struct == "double 1,2 2-ribbon"){
    # start with two 1-ribbon chains:
    halfPt = floor(n/2)
    chain1 = 1:halfPt
    chain2 = n:(halfPt+1)
    # add ties across the chains
    for (i in 1:length(chain1)) {
      net[chain1[i], chain2[i]] = 1
    }
    # add ties along each chain
    for (i in 1:length(chain1)-1) {
      net[chain1[i], chain1[i+1]] = 1
      net[chain2[i], chain2[i+1]] = 1
    }
    # 1,2 2-RIBBON
    if(struct == "1,2 2-ribbon"){
    # for 1,2 2-ribbon, add 1,2 ties across chains
      for (i in 1:(length(chain1)-1)) {
        net[chain1[i], chain2[i+1]] = 1
      }
      if(mode == "cycle"){
        net[chain1[1], chain1[halfPt]] = 1
        net[chain2[1], chain2[halfPt]] = 1
        net[chain1[halfPt], n] = 1
      }
    }
    # DOUBLE 1,2 2-RIBBON
    if(struct == "double 1,2 2-ribbon"){
      # for 1,2 2-ribbon, add 1,2 ties across chains
      for (i in 1:length(chain1)-1) {
        net[chain1[i], chain2[i+1]] = 1
        net[chain2[i], chain1[i+1]] = 1
      }
      if(mode == "cycle"){
        net[chain1[1], chain1[halfPt]] = 1
        net[chain2[1], chain2[halfPt]] = 1
        net[chain1[halfPt], n] = 1
        net[1, halfPt+1] = 1
      }
    }
  }
  return(net)
}

