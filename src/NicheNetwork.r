## ---------------------------
##
## Script name: simulations.r
##
## Purpose of script: This script implements the niche model (originally proposed by Williams and Martinez (Nature, 2000)
##
## Author: Dr Miguel Lurgi
## Lecturer in Biosciences (Computational Ecology)
## Computational Ecology Lab - Department of Biosciences
## Swansea University, UK
## 
## and
##
## Centre for Biodiversity Theory and Modelling
## Theoretical and Experimental Ecology Station, CNRS, France
##
## Date Created: November-2017
##
## Copyright (c) Miguel Lurgi, 2017-2021
## Email: miguel.lurgi@swansea.ac.uk
##
## ---------------------------
##
## Notes:
##
## This script is provided as supplementary material for the paper:
## Galiana, N., Lurgi, M., Claramunt-López, B. et al. The spatial scaling of species interaction networks. 
## Nature Ecology & Evolution 2, 782–790 (2018). https://doi.org/10.1038/s41559-018-0517-3
##
## ---------------------------

connectance <- function(S){
  return(0.8*(S**(-0.5)))
}

NicheNetwork <- function(S, C){
  require(igraph)
  require(sna)
  
  S <- S;
  x <- 0;
  y <- S;
  dup_sp <- FALSE;
  
  while(x < S | y > 1 | dup_sp){
    dup_sp <- FALSE;
    M <- matrix(0, S, S);
    
    # first we obtain the n value for all the species
    n <- sort(runif(S, 0, 1));
    
    #we then obtain the feeding range for each of the species, drawn from a beta
    #distribution with beta value = (1/(2*C)) - 1
    beta <- (1/(2*C)) - 1;
    r <- rbeta(S, 1, beta) * n;    
    
    #we enforce the species with the lowest niche value to be a basal species
    r[1] <- 0;
    
    #the centre of the feeding range is drawn from an uniform distribution
    #between ri/2 and min(ni, 1-ri/2)
    c <- numeric(S);
    for(i in 1:S){
      c[i] <- runif(1, r[i]/2, min(n[i], (1-r[i]/2)));
      
      #once we have the r's and the c's we can obtain the boundaries of the feeding range
      offset <- r[i]/2;
      upper_bound <- c[i] + offset;
      lower_bound <- c[i] - offset;
      
      for(j in 1:S){
        if(n[j] > lower_bound & n[j] < upper_bound)
          M[j,i] = 1;
      } 
    }
    
    #we verify that the network (i) is connected and (2) does not have cycles with no external energy input
    M_temp <- M;
    diag(M_temp) <- 0;
    
    graf <- igraph::graph.adjacency(M);
    y <- igraph::no.clusters(graf, mode='weak');
    if(y > 1){
      #print("NicheNetwork :: the network is not connected...");
      next;
    }
    
    clts <- igraph::clusters(graf, mode='strong');
    x <- clts$no;
    if(x < S){
      clts_nos <- which(clts$csize > 1);  
      cycles_no_input <- FALSE;
      
      for(c in clts_nos){
        members <- which(clts$membership == c);
        cluster_ok <- FALSE;
        for(m in members){
          prey <- neighbors(graf, m, mode='in');
          if( length(intersect(prey, members))  > length(members) + 1 ){
            ## if this happens, this cycle/cluster has external energy input
            cluster_ok <- TRUE;
            break;
          }
        }
        
        if(!cluster_ok){
          print("NicheNetwork :: the network has cycles with no external energy input...");
          cycles_no_input <- TRUE;
          break;
        }
      }
      
      if(cycles_no_input){
        next;
      }else{
        x <- S;
      }
      
    }
    
    #and we also verify that there are not duplicate species
    for(i in 1:S){
      if(dup_sp){
        break;
      }
      preys_i <- M[,i];
      predators_i <- M[i,];
      for(j in 1:S){
        if(i == j) next;
        sim_prey <- preys_i == M[,j];
        sim_preds <- predators_i == M[j,];
        
        if(sum(sim_prey) == S && sum(sim_preds) == S ){
          dup_sp <- TRUE;
          #print("NicheNetwork :: there is a duplicate species");
          break;
        }   
      }
      
      #as long as we are traversing the network we might as well check
      #for basal species with cannibalistic links...
      
      if(sum(preys_i) == 1 && M[i,i] == 1){
        print("NicheNetwork :: a cannibalistic link was found in a basal species... removing it");
        M[i,i] <- 0;
      } 
    }
  }  
  
  ########################################################################################
  
  return(list(S=dim(M)[1], C=sum(M)/((dim(M)[1])**2), M=M, niche=n, centre=c, radius=r));
  
}
