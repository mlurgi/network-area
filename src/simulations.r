

source('NicheNetwork.r')
source('utils.r')
source('space.r')

require(cheddar)
require(igraph)
require(betapart)

# Load function for multiple site measures of beta diversity

## Note that this function is capable of calculatig all three partitions (i.e. it includes 
## Baselga's multiple-site partition), but these must be specified using the index.family 
## argument of Baselga's beta.multi (i.e. 'carvalho' is default). This is directly 
## modified from the betapart package of Baselga et al. (2013), see also: 
## http://127.0.0.1:20382/library/betapart/html/betapart-package.html

beta.multi.carv <- function (x, index.family = "carvalho") 
{library(betapart)
  index.family <- match.arg(index.family, c("jaccard", "sorensen","carvalho"))
  if (!inherits(x, "betapart")) {
    x <- betapart.core(x)
  }
  maxbibj <- sum(x$max.not.shared[lower.tri(x$max.not.shared)])
  minbibj <- sum(x$min.not.shared[lower.tri(x$min.not.shared)])
  switch(index.family, sorensen = {
    beta.sim <- minbibj/(minbibj + x$a)
    beta.sne <- (x$a/(minbibj + x$a)) * ((maxbibj - minbibj)/((2 * 
                                                                 x$a) + maxbibj + minbibj))
    beta.sor <- (minbibj + maxbibj)/(minbibj + maxbibj + 
                                       (2 * x$a))
    multi <- list(beta.SIM = beta.sim, beta.SNE = beta.sne, 
                  beta.SOR = beta.sor)
  }, jaccard = {
    beta.jtu <- (2 * minbibj)/((2 * minbibj) + x$a)
    beta.jne <- (x$a/((2 * minbibj) + x$a)) * ((maxbibj - 
                                                  minbibj)/((x$a) + maxbibj + minbibj))
    beta.jac <- (minbibj + maxbibj)/(minbibj + maxbibj + 
                                       x$a)
    multi <- list(beta.JTU = beta.jtu, beta.JNE = beta.jne, 
                  beta.JAC = beta.jac)
  }, carvalho = {
    beta.3 <- (2 * minbibj)/(minbibj + maxbibj + x$a)
    beta.rich <- ((maxbibj - minbibj)/(minbibj + maxbibj + x$a))
    beta.cc <- (minbibj + maxbibj)/(minbibj + maxbibj + x$a)
    multi <- list(beta.3multi = beta.3, beta.RICH = beta.rich, 
                  beta.CC = beta.cc)
  })
  return(multi)
}

#shape of the community
islands <- 75
#S <- 50
#C <- connectance(S)
landscape <- space(islands, 0.3)[[2]]
landscape[landscape > 0] <- 1
XY <- space(islands, 0.3)[[1]]

#parameters
c <- 0.1    #colonisation rate
e <- 0.4    #extinction rate
#d <- 0.1    #dispersal rate


# extinction_values <- c(.2, .3, .4, .5, .6, .7, .8)
# colonisations <- c(.01, .05, .1, .15, .2, .25, .3, .35, .4)
dispersals <- c(0,.01,.05,.1,.2,.3,.4,.5)

timesteps <- 1000
replicates <- 5
max_attempts <- 100

output <- NULL
output_occupancy <- NULL

# for(c in colonisations){
#   print(c)
#   for(e in extinction_values){
#     print(e)
    for(d in dispersals){
      print(d)
      for(j in 1:replicates){
        print(j)
        load(paste('./networks/network-',j,sep=''))
        S <- food_web$S
        presences <- array(0, c(islands,S))
        for(i in 1:timesteps){
          #basal species
          basal_sps <- which(colSums(food_web$M) == 0)
          for(isl in 1:islands){
            #this is the current community
            cur_com <- presences[isl,]
            #this is what species are present
            present <- which(cur_com != 0)
            #this is the set of species in the network that have prey in the current community
            with_prey <- which(colSums(as.matrix(food_web$M[present, ])) > 0)
            #this are the potential colonisers
            potential <- append(basal_sps , with_prey)
            #to which we remove the ones that are already present
            potential <- setdiff(potential, present)
            
            if(length(potential) != 0){
              colonisers <- potential[which(runif(length(potential), 0, 1) < c)]
              presences[isl, colonisers] <- 1
            }
            
            extinctions <- FALSE
            if(length(present) != 0){
              extinct <- present[which(runif(length(present), 0, 1) < e)]
              if(length(extinct) != 0){
                presences[isl, extinct] <- 0
                extinctions <- TRUE
              }
            }
            
            if(extinctions){
              cur_com <- presences[isl,]
              #this is what species are present
              present <- which(cur_com != 0)
              #this is the set of species in the network that have prey in the current community
              with_prey <- which(colSums(as.matrix(food_web$M[present, ])) > 0)
              legal <- union(with_prey, basal_sps)
              illegal <- present[!(present %in% legal)]
              presences[isl, illegal] <- 0  
            }
            
          }
          
          
          ### dispersal after colonisation and extinction processes
          
          #### we took this code from the metacomm model. Check Nuria's notes to see what they do
          randCol = matrix(runif(islands*S,0,1),nr=islands,nc=S)
          ColMat = matrix(0,nr=islands,nc=S)
          
          ConPop = landscape%*%presences
          LocalCol = 1-(1-d)^ConPop
          
          LocalPrey = presences %*% food_web$M
          LocalPrey[LocalPrey>0] = 1
          LocalPrey[,basal_sps] = 1
          
          ColMat[presences == 0 & LocalPrey == 1 & randCol<LocalCol] = 1
          
          presences = presences + ColMat
          
        }
        
        
        
        results <- presences # This is the community matrix presences
        net <- food_web
        
        occupancy <- c()
        dietbreath <- c()
        dispersal <- c()
        occupancy <- append(occupancy, colMeans(results))
        dietbreath <- append(dietbreath, colSums(net$M))
        dispersal <- append(dispersal, d)
        
        cur_out_occu <- data.frame(occupancy,dietbreath,rep(dispersal,length(occupancy)))
      
      
      if(is.null(output_occupancy)){
        output_occupancy <- cur_out_occu
      }else{
        output_occupancy <- rbind(output_occupancy, cur_out_occu)
      } 
        
        
        beta_total <- c()
        beta_turn <- c()
        beta_rich <- c()
        beta_total <- append(beta_total, as.numeric(beta.multi.carv(results)[3]))
        beta_turn <- append(beta_turn, as.numeric(beta.multi.carv(results)[1]))
        beta_rich <- append(beta_rich, as.numeric(beta.multi.carv(results)[2]))
        
        #### after obtaining the results we aggregate local communities and obtain network attributes 
        species <- c()
        area <-c()
        C <- c()
        links <- c()
        mfcls <- c()
        indegree <- c()
        normalised_indegree <- c()
        outdegree <- c()
        normalised_outdegree <- c()
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
        SD_Gen <- c()
        SD_Vul <- c()
        
        ## Order patches across space
        euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
        focal_coordinate <- c(0.5,0.5)
        distance_focal <- c()
        for (i in 1:islands){
          distance_focal <- append(distance_focal, euc.dist(focal_coordinate,XY[i,]))
        }
        
        results_ordered <- results[order(distance_focal),] #this orders the metacommunity matrix in respect to the center of the lanscape (focal_coordinate)
        
        for(i in 1:islands){
          areas <- i
#           # aggregation random
#           patches_selected <- as.matrix(results[sample(nrow(results), i),])
#           
#           if(i > 1){
#             current_species_ran <- which(colSums(patches_selected[,]) > 0)
#             
#           }else{
#             current_species_ran <- which(patches_selected[,] > 0)
#             
#           }
#           if(length(unique(current_species_ran)) <= 1){
#             next
#           }
#           
#           M <- net$M[current_species_ran, current_species_ran]
#           
#           #we can check if the network is connected
#           graf <- igraph::graph.adjacency(M); # Making sure the network is connected
#           y <- igraph::no.clusters(graf, mode='weak');
#           attempts <- 0
#           while(y > 1 & attempts < max_attempts){
#             patches_selected <- as.matrix(results[sample(nrow(results), i),])
#             
#             if(i > 1){
#               current_species_ran <- which(colSums(patches_selected[,]) > 0)
#               
#             }else{
#               current_species_ran <- which(patches_selected[,] > 0)
#               print(i)
#             }    
#             M <- net$M[current_species_ran, current_species_ran]
#             graf <- igraph::graph.adjacency(M); # Making sure the network is connected
#             y <- igraph::no.clusters(graf, mode='weak');
#             
#             attempts <- attempts + 1
#           }
#           
#           if(length(unique(current_species_ran)) < 2 | attempts >= max_attempts){
#             print('network not found')
#             next
#           }
#           
#               if(nrow(M) < 2){
#                 next
#               }
          
#          write.table(M, paste(i ,"localmetacom_ran", ".txt", sep = ""), row.names=FALSE, col.names=FALSE) #Saving the local network
          
          ## aggregation ordering them
          if(i > 1){
            current_species_ran <-unique(which(colSums(results_ordered[1:i,]) > 0))
          }else{
            current_species_ran <- unique(which(results_ordered[1,] > 0))
          }
          if(length(unique(current_species_ran)) <= 1){
            next
          }
          M <- net$M[current_species_ran, current_species_ran]
          
          #we can check if the network is connected
          graf <- igraph::graph.adjacency(M); # Making sure the network is connected
          y <- igraph::no.clusters(graf, mode='weak');
          attempts <- 0
          while(y > 1 & attempts < max_attempts){
          patches_selected <- as.matrix(results[sample(nrow(results), i),])
  
          if(i > 1){
          current_species_ran <- which(colSums(patches_selected[,]) > 0)
    
          }else{
          current_species_ran <- which(patches_selected[,] > 0)
          print(i)
          }    
          M <- net$M[current_species_ran, current_species_ran]
         graf <- igraph::graph.adjacency(M); # Making sure the network is connected
         y <- igraph::no.clusters(graf, mode='weak');
  
         attempts <- attempts + 1
         }

         if(length(unique(current_species_ran)) < 2 | attempts >= max_attempts){
          print('network not found')
          next
          }

          if(nrow(M)<2){
            next
          }

          
          species <- append(species, length(current_species_ran))
          C <- append(C, sum(M)/((dim(M)[1])*(dim(M)[1]-1))) #Connectance if we prevent cannibalism
          
          mfcl <- tryCatch({
            MeanFoodChainLength(M)
          }, warning = function(w) {
            'NA'
          }, error = function(e) {
            'NA'
          }, finally = {
            
          })
          
          mfcls <- append(mfcls, mfcl)
          area <- append(area, areas)
          links <- append(links, sum(M))
          indegree <- append(indegree, Generality(M))
          normalised_indegree <- append(normalised_indegree, NormalisedGenerality(M))
          outdegree <- append(outdegree, Vulnerability(M))
          normalised_outdegree <- append(normalised_outdegree, NormalisedVulnerability(M))
          basal <- append(basal, FractionOfBasal(M))
          top <- append(top, FractionOfTop(M))
          intermediate <- append(intermediate, FractionOfIntermediate(M))
          sd_gen <- append(sd_gen, SDGenerality(M))
          sd_vul <- append(sd_vul, SDVulnerability(M))
          overlap <- append(overlap, CalculatePredatorOverlap(M))
          omnivory <- append(omnivory, Omnivory(M))
          S_top <- append(S_top, NumberOfTop(M))
          S_intermediate <- append(S_intermediate, NumberOfIntermediate(M))
          S_basal <- append(S_basal, NumberOfBasal(M))
          SD_Vul <- append(SD_Vul, SDVulnerability(M))
          SD_Gen <- append(SD_Gen, SDGenerality(M))


          cur_out <- data.frame(rep(c,length(species)), rep(e,length(species)), rep(d,length(species)), rep(j,length(species)), beta_total, beta_rich, beta_turn, species, C, links, area,  indegree, normalised_indegree, outdegree, mfcls, normalised_outdegree, basal, top, intermediate, overlap, sd_gen, sd_vul, omnivory, S_top, S_intermediate, S_basal, SD_Gen, SD_Vul)
 
          
        }
        
        if(is.null(output)){
          output <- cur_out
        }else{
          output <- rbind(output, cur_out)
        } 
        
      }
      
      
     }
#   }
# }
# 

write.csv(output, file='output_final_mfcl.csv')
write.csv(output_occupancy, file='output_occupancy_final.csv')

