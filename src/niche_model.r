
### niche model experiment
require(igraph)
require(sna)

sar <- function(area){
  return(10*(area**0.27))
}

connectance <- function(S){
  return(0.8*(S**(-0.5)))
}

source('utils.r')

#areas <- append(seq(.1,1,.1), 1:80)
areas <- seq(1,100000,1)
sps <- as.integer(sar(areas))




species_area <- as.data.frame(cbind(areas,sps))

plot(areas,sps)
plot(log10(areas),log10(sps))

#sps <- unique(sps)
sps <- seq(20,189,1)
c_s <- connectance(sps)

simulations <- 50

output <- NULL

for(j in 1:simulations){ 
  print(j)
  load(paste('./networks/network-',j,sep=''))
  n_regional <- food_web #this is the regional network
  
  species <- c()
  area <- c()
  C <- c()
  links <- c()
#  mfcls <- c()
  indegree <- c()
  outdegree <- c()
  indegree_normalized <- c()
  outdegree_normalized <- c()
  basal <- c()
  top <- c()
  intermediate <- c()
  overlap <- c()
  sd_gen <- c()
  sd_vul <- c()
  omnivory <- c()
  S_top <- c()
  S_intermediate <- c()
  S_basal <- c()
  

  for(i in 1:length(sps)){ #this loop is to generate the local networks from the regional one based on SAR and C
    #this networks are the ones we want to analise. 
    print(paste('species',sps[i])) 
    
    connectance<-c_s[i] #Connectance associated to the given number of species 
    l<-floor(c_s[i]*sps[i]*(sps[i]-1)) #Number of links required to have the desired connectance
    
   cols <- sample(ncol(n_regional$M), sps[i]) # We pick randomly the number of species we want to create the local network from the regional one
   
   n <- n_regional$M[cols,cols] #this is the local network
 
   
# if we don't want to fix the S-C relationship, we comment this section.   
    while(sum(n) < l){
      
    cols = sample(ncol(n_regional$M), sps[i]) # We pick randomly the number of species we want to create the local network from the regional one
    #rows= sample(nrow(n_regional$M), species[i])
    n = n_regional$M[cols,cols] #this is the local network
    
    }

    graf <- igraph::graph.adjacency(n); # Making sure the network is connected
    y <- igraph::no.clusters(graf, mode='weak');
    while(y > 1){
      cols <- sample(ncol(n_regional$M), sps[i]) # We pick randomly the number of species we want to create the local network from the regional one
      n <- n_regional$M[cols,cols] #this is the local network
      graf <- igraph::graph.adjacency(n); # Making sure the network is connected
      y <- igraph::no.clusters(graf, mode='weak');
    }
    
   C <- append(C,  sum(n)/(dim(n)[1])^2)
   
#    mfcl <- tryCatch({
#      MeanFoodChainLength(n)
#    }, warning = function(w) {
#      'NA'
#    }, error = function(e) {
#      'NA'
#    }, finally = {
#      
#    })
   
   species <- append(species, dim(n)[1])
   area <- append(area, areas[i])
   links <- append(links, sum(n))
   #mfcls <- append(mfcls, mfcl)
   indegree <- append(indegree, Generality(n))
   outdegree <- append(outdegree, Vulnerability(n))
   outdegree_normalized <- append(outdegree_normalized, NormalisedVulnerability(n))
   indegree_normalized <- append(indegree_normalized, NormalisedGenerality(n))
   basal <- append(basal, FractionOfBasal(n))
   top <- append(top, FractionOfTop(n))
   intermediate <- append(intermediate, FractionOfIntermediate(n))
   sd_gen <- append(sd_gen, SDGenerality(n))
   sd_vul <- append(sd_vul, SDVulnerability(n))
   overlap <- append(overlap, CalculatePredatorOverlap(n))
   omnivory <- append(omnivory, Omnivory(n))
   S_basal <- append(S_basal, NumberOfBasal(n))
   S_top <- append(S_top, NumberOfTop(n))
   S_intermediate <- append(S_intermediate, NumberOfIntermediate(n))
  
   
  }
  
  cur_out <- data.frame(rep(i,length(species)), species, area, C, links, indegree, outdegree, indegree_normalized, outdegree_normalized, basal, top, intermediate, overlap, sd_gen, sd_vul, omnivory, S_intermediate, S_top, S_basal)
 #mfcls,  
  if(is.null(output)){
    output <- cur_out
  }else{
    output <- rbind(output, cur_out)
  }
  
  
    
   
}

write.csv(output, file='output_nar_withoutbasal_fixC2.csv')

#plot(species, mfcls)
#plot(areas, mfcls)
