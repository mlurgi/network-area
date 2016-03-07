

setwd('/Users/nuriagaliana/Documents/network_area_new')
source('NicheNetwork.R')
source('utils_nu.R')

require(cheddar)
require(igraph)


#parameters
c <- 0.3   #colonisation rate
extinction <- c(.1, .12, .15, .2, .3, .4, .5, .6, .7, .8, .95) #extinction rate


timesteps <- 1000
replicates <- 100

output <- NULL


for(e in extinction){
  for(j in 1:replicates){
    load(paste('./networks/network-',j,sep=''))
    S <- food_web$S
    presences <- array(0, S)
    
    print(j)
    for(i in 1:timesteps){
      #basal species
      basal_sps <- which(colSums(food_web$M) == 0)
      cur_com <- presences
      
      #this is what species are present
      present <- which(cur_com != 0)
      
      if(sum(cur_com) != 0){
        
        #this is the set of species in the network that have prey in the current community
        with_prey <- which(colSums(as.matrix(food_web$M[present, ])) > 0)
        #this are the potential colonisers
        potential <- append(basal_sps , with_prey)
        #to which we remove the ones that are already present
        potential <- setdiff(potential, present)
        
      }else{
        potential <- basal_sps
      }
      
      if(length(potential) != 0){
        colonisers <- potential[which(runif(length(potential), 0, 1) < c)]
        presences[colonisers] <- 1
      }
      
      extinctions <- FALSE
      if(length(present) != 0){
        extinct <- present[which(runif(length(present), 0, 1) < e)]
        if(length(extinct) != 0){
          presences[extinct] <- 0
          extinctions <- TRUE
        }
      }
      
      if(extinctions){
        cur_com <- presences
        #this is what species are present
        present <- which(cur_com != 0)
        #this is the set of species in the network that have prey in the current community
        with_prey <- which(colSums(as.matrix(food_web$M[present, ])) > 0)
        legal <- union(with_prey, basal_sps)
        illegal <- present[!(present %in% legal)]
        presences[illegal] <- 0  
      }
      
    }
    
    results <- presences
    net <- food_web
    
    current_species <- which(results > 0)
    
    if(length(current_species) < 2){
      next
    }
    
    M <- net$M[current_species, current_species]
    
    #we can check if the network is connected
    graf <- igraph::graph.adjacency(M); # Making sure the network is connected
    y <- igraph::no.clusters(graf, mode='weak');
    attempts <- 0
    #     if(y > 1 ){
    #       print('network not found')
    #       next
    #     }
    
    C <- sum(M)/((dim(M)[1])*(dim(M)[1]-1)) #Connectance if we prevent cannibalism
    
    mfcl <- tryCatch({
      MeanFoodChainLength(M)
    }, warning = function(w) {
      'NA'
    }, error = function(e) {
      'NA'
    }, finally = {
      
    })
    
    mfcls <- mfcl
    links <- sum(M)
    indegree <- Generality(M)
    normalised_indegree <- NormalisedGenerality(M)
    outdegree <- Vulnerability(M)
    normalised_outdegree <- NormalisedVulnerability(M)
    basal <- FractionOfBasal(M)
    top <- FractionOfTop(M)
    intermediate <- FractionOfIntermediate(M)
    sd_gen <- SDGenerality(M)
    sd_vul <- SDVulnerability(M)
    overlap <- CalculatePredatorOverlap(M)
    omnivory <- Omnivory(M)
    S_top <- NumberOfTop(M)
    S_intermediate <- NumberOfIntermediate(M)
    S_basal <- NumberOfBasal(M)
    ratio <- c/e
    colonization_rate <- c
    extinction_rate <- e
    
    cur_out <- data.frame(replicate=j, S=length(current_species), C, links,  mfcls, indegree, normalised_indegree, outdegree, normalised_outdegree, basal, top, intermediate, overlap, sd_gen, sd_vul, omnivory, S_top, S_intermediate, S_basal, ratio, colonization_rate, extinction_rate)
    
    if(is.null(output)){
      output <- cur_out
    }else{
      output <- rbind(output, cur_out)
    } 
    
  }
}

write.csv(output, file='output_ttib_rates..csv')




