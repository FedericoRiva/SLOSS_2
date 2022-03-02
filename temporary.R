library(betapart)
prova <- c(1,1,1,1,1,1,1,0,0,0,0,0)
beta.multi(prova)

object <- list_simulations_ten

#alpha_set <-function(object){ # object is one of the lists from SETS_PATCHES
  
  LIST_SIM <- list()
  for (i in 1:100) {LIST_SIM[[i]] <- object}
  
  # merge the simulated sets of patches with the simulated assemblages in each patch
  for (i in 1: length(LIST_SIM)) {
    for (j in 1: length(LIST_SIM[[i]])) {
      for (k in 1: length(LIST_SIM[[i]][[j]])){
        
        if (length(LIST_SIM[[i]][[j]]) > 1) {
          
          LIST_SIM[[i]][[j]][[k]] <- merge(LIST_SIM[[i]][[j]][[k]], 
                                           list_sampled_data[[i]][j][[1]],
                                           by = "study_and_patch")
          
        } else  {
          
          LIST_SIM[[i]][[j]] <- NA
        }
        
      }
    }
  }
  
  
  
  # calculate the average patch areas in each set of patches
  LIST_SIM_area <- LIST_SIM
  
  for (i in 1: length(LIST_SIM_area)) {
    for (j in 1: length(LIST_SIM_area[[i]])) {
      for (k in 1: length(LIST_SIM_area[[i]][[j]])){
        
        if (!is.na(LIST_SIM_area[[i]][[j]])) {
          
          LIST_SIM_area[[i]][[j]][[k]] <- mean(LIST_SIM_area[[i]][[j]][[k]][,ncol(LIST_SIM_area[[i]][[j]][[k]])])
          
        } 
        
      }
    }
  }
  # warnings appears because every list has more than one element, so the !is.na returns multiple TRUE statements when a list is not NA
  
  LIST_SIM_area_sum <- LIST_SIM
  
  for (i in 1: length(LIST_SIM_area_sum)) {
    for (j in 1: length(LIST_SIM_area_sum[[i]])) {
      for (k in 1: length(LIST_SIM_area_sum[[i]][[j]])){
        
        if (!is.na(LIST_SIM_area_sum[[i]][[j]])) {
          
          LIST_SIM_area_sum[[i]][[j]][[k]] <- sum(LIST_SIM_area_sum[[i]][[j]][[k]][,ncol(LIST_SIM_area_sum[[i]][[j]][[k]])])
          
        } 
        
      }
    }
  }
  
  
  
  # transform lists into tabels
  for (i in 1:length(LIST_SIM_area)){
    for ( j in 1: length(LIST_SIM_area[[i]])) {
      
      if (!is.na(LIST_SIM_area[[i]][[j]])) {
        LIST_SIM_area[[i]][[j]] <- do.call(rbind.data.frame, LIST_SIM_area[[i]][[j]])
        colnames(LIST_SIM_area[[i]][[j]]) <- "mean_patch_area"
      }
    }
  }
  
  for (i in 1:length(LIST_SIM_area_sum)){
    for ( j in 1: length(LIST_SIM_area_sum[[i]])) {
      
      if (!is.na(LIST_SIM_area_sum[[i]][[j]])) {
        LIST_SIM_area_sum[[i]][[j]] <- do.call(rbind.data.frame, LIST_SIM_area_sum[[i]][[j]])
        colnames(LIST_SIM_area_sum[[i]][[j]]) <- "total_patch_area"
      }
    }
  }
  
  
  # calculate the species richness in each patch
  LIST_SIM_richness <- LIST_SIM
  
  for (i in 1: length(LIST_SIM_richness)) {
    for (j in 1: length(LIST_SIM_richness[[i]])) {
      for (k in 1: length(LIST_SIM_richness[[i]][[j]])){
        
        if (!is.na(LIST_SIM_richness[[i]][[j]])) {
          # take the column of each table from column two to column n-1, take their sum, and count how many times their colsum is larger than 0
          LIST_SIM_richness[[i]][[j]][[k]] <- sum(colSums(LIST_SIM_richness[[i]][[j]][[k]][,2: (ncol(LIST_SIM_richness[[i]][[j]][[k]])-1)]) > 0)
        }
        
      }
    }
  }
  
  
  
  
  for (i in 1:length(LIST_SIM_richness)){
    for ( j in 1: length(LIST_SIM_richness[[i]])) {
      
      if (!is.na(LIST_SIM_richness[[i]][[j]])) {
        LIST_SIM_richness[[i]][[j]] <- do.call(rbind.data.frame, LIST_SIM_richness[[i]][[j]])
        colnames(LIST_SIM_richness[[i]][[j]]) <- "richness"
      }
    }
  }
  
  
  
  
  
  # calculate evenness in each patch
  LIST_SIM_evenness <- LIST_SIM
  
  for (i in 1: length(LIST_SIM_evenness)) {
    for (j in 1: length(LIST_SIM_evenness[[i]])) {
      for (k in 1: length(LIST_SIM_evenness[[i]][[j]])){
        
        if (!is.na(LIST_SIM_evenness[[i]][[j]])) {
          # use the HurlbertPIE function on every table, excluding patch name (column 1) and patch size (column nrow)
          LIST_SIM_evenness[[i]][[j]][[k]] <- paleotree::HurlbertPIE(colSums(LIST_SIM_evenness[[i]][[j]][[k]][,2: (ncol(LIST_SIM_evenness[[i]][[j]][[k]])-1)]))
        } 
        
      }
    }
  }
  
  
  
  for (i in 1:length(LIST_SIM_evenness)){
    for ( j in 1: length(LIST_SIM_evenness[[i]])) {
      
      if (!is.na(LIST_SIM_evenness[[i]][[j]])) {
        LIST_SIM_evenness[[i]][[j]] <- do.call(rbind.data.frame, LIST_SIM_evenness[[i]][[j]])
        colnames(LIST_SIM_evenness[[i]][[j]]) <- "evenness"
      }
    }
  }
  
  
  
  LIST_SIM_richness_scaled <- LIST_SIM_richness
  
  for (i in 1:length(LIST_SIM_richness_scaled)){
    for ( j in 1: length(LIST_SIM_richness_scaled[[i]])) {
      
      if (!is.na(LIST_SIM_richness_scaled[[i]][[j]])) {
      
          if (var(LIST_SIM_richness_scaled[[i]][[j]]$richness) == 0) {
              LIST_SIM_richness_scaled[[i]][[j]]$richness <- rep(0, length(LIST_SIM_richness_scaled[[i]][[j]]$richness))
              } else {
              LIST_SIM_richness_scaled[[i]][[j]]$richness <- ((LIST_SIM_richness_scaled[[i]][[j]]$richness - mean(LIST_SIM_richness_scaled[[i]][[j]]$richness))/sd(LIST_SIM_richness_scaled[[i]][[j]]$richness))
          }
      }
    }
  }
  
  

  LIST_SIM_evenness_scaled <- LIST_SIM_evenness
  
  for (i in 1:length(LIST_SIM_evenness_scaled)){
    for ( j in 1: length(LIST_SIM_evenness_scaled[[i]])) {
      
      if (!is.na(LIST_SIM_evenness_scaled[[i]][[j]])) {
        
        if (var(LIST_SIM_evenness_scaled[[i]][[j]]$evenness) == 0) {
          LIST_SIM_evenness_scaled[[i]][[j]]$evenness <- rep(0, length(LIST_SIM_evenness_scaled[[i]][[j]]$evenness))
        } else {
          LIST_SIM_evenness_scaled[[i]][[j]]$evenness <- ((LIST_SIM_evenness_scaled[[i]][[j]]$evenness - mean(LIST_SIM_evenness_scaled[[i]][[j]]$evenness))/sd(LIST_SIM_evenness_scaled[[i]][[j]]$evenness))
        }
      }
    }
  }
  
  
  
  LIST_SIM_combined <- LIST_SIM_area
  
  for (i in 1:length(LIST_SIM_combined)){
    for ( j in 1: length(LIST_SIM_combined[[i]])) {
      
      if (!is.na(LIST_SIM_combined[[i]][[j]])) {
        LIST_SIM_combined[[i]][[j]] <- cbind( LIST_SIM_area[[i]][[j]],
                                              LIST_SIM_area_sum[[i]][[j]],
                                              LIST_SIM_richness[[i]][[j]],  
                                              LIST_SIM_evenness[[i]][[j]],
                                              LIST_SIM_richness_scaled[[i]][[j]],
                                              LIST_SIM_evenness_scaled[[i]][[j]])
        LIST_SIM_combined[[i]][[j]]$simulation_number <- rep(i, # give a number from 1 to 100 
                                                             nrow(LIST_SIM_combined[[i]][[j]]))
        LIST_SIM_combined[[i]][[j]]$study <- rep(names(list_sampled_data[[i]][j]), 
                                                 nrow(LIST_SIM_combined[[i]][[j]]))
        
        LIST_SIM_combined[[i]][[j]]$patch_set_number = 1:nrow(LIST_SIM_combined[[i]][[j]])
        
      }
    }
  }
  
  
  ## combine all studies in a table for each of the 100 simulations
  for (i in 1:length(LIST_SIM_combined)){
    LIST_SIM_combined[[i]] <- do.call(rbind.data.frame, LIST_SIM_combined[[i]])
    
  }
  
  LIST_SIM_combined <- do.call(rbind, LIST_SIM_combined)
  LIST_SIM_combined <- na.omit(LIST_SIM_combined)
  
  final <- LIST_SIM_combined
  #final_twenty <- LIST_SIM_combined
#} # end of the function
  
  
  ## ADD IUCN DATA
