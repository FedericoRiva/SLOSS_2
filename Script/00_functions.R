## FUNCTION 1
# SAMPLEs ROWS ON A PER-PATCH BASIS TO HAVE # OF INDIVIDUALS PROPORTIONAL TO PATCH SIZE
get_rows <- function(df, rows) df[rows, , drop = FALSE]

##########################################################
##########################################################
##########################################################

## FUNCTION 2
# CONVERTS LIST OF VECTORS OF SPECIES ABUNDANCES INTO TABLE
to_table <- function(list_to_table) {
  table_prova <- rbindlist(lapply(list_to_table, function(x) as.data.frame.list(x)), fill=TRUE)
  table_prova[is.na(table_prova)] <- 0
  table_prova$study_and_patch <- names(list_to_table)
  list_to_table <- table_prova
}

##########################################################
##########################################################
##########################################################

## FUNCTION 3
# GENERATES SETS OF PATCHES AT EQUAL HABITAT AREA
SETS_PATCHES <- function(db, percent) {  # db= dataset, percent = habitat amount
  
  set.seed(111)
  n_db <- 100 # 100 attempts to simulate random sets of patches
  # total habitat amount 
  target_area <- sum(db$frag_size_num) * percent
  
  # check how big is the largest patch in a dataset in comparison to the target area
  target_area
  max(db$frag_size_num)
  
  ## randomly select a n_db sets of patches 
  list_db <- list()
  for (i in 1 : n_db){
    sampled_db <- dplyr::sample_n(db,1) #sample one row
    db2 <- db[db$frag_size_num != sampled_db$frag_size_num,] #remove that row from the table
    while( sum(sampled_db$frag_size_num) < target_area) {
      sampled_db2 <- dplyr::sample_n(db2,1)
      sampled_db <- rbind(sampled_db,sampled_db2) #sample rows until target area, creating a database
      db2 <- db2[db2$frag_size_num != c(sampled_db2$frag_size_num),] #remove that row from the table
    }
    list_db[[i]] <-  sampled_db
  }
  #
  
  for (i in 1:n_db) {list_db[i][[1]] <- list_db[i][[1]][order(list_db[i][[1]]$frag_size_num),]}
  
  list_db <- unique(list_db)
  list_db <- rlist::list.filter(list_db, sum(frag_size_num) < (target_area + target_area*0.05) ) # a tolerance of 5% of the target area
  
  
  for (i in 1: length(list_db)) {
    list_db[[i]] <- as.data.frame(list_db[[i]]$study_and_patch)
    colnames(list_db[[i]]) <- "study_and_patch"
  }
  
  print(list_db)
  
}

##########################################################
##########################################################
##########################################################

## FUNCTION 4
# CALCULATE ALPHA DIVERSITY METRICS ACROSS SETS OF PATCHES
alpha_set <-function(object){ # object is one of the lists from SETS_PATCHES
  
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
  
  #calculate the sum of all patch areas within a set of patches
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
  
  # scale areas
  LIST_SIM_area_scaled <- LIST_SIM_area
  
  for (i in 1:length(LIST_SIM_area_scaled)){
    for ( j in 1: length(LIST_SIM_area_scaled[[i]])) {
      
      if (!is.na(LIST_SIM_area_scaled[[i]][[j]])) {
        
        if (var(LIST_SIM_area_scaled[[i]][[j]]$mean_patch_area) == 0) {
          LIST_SIM_area_scaled[[i]][[j]]$mean_patch_area_sc <- rep(0, length(LIST_SIM_area_scaled[[i]][[j]]$mean_patch_area))
        } else {
          LIST_SIM_area_scaled[[i]][[j]]$mean_patch_area_sc <- ((LIST_SIM_area_scaled[[i]][[j]]$mean_patch_area - mean(LIST_SIM_area_scaled[[i]][[j]]$mean_patch_area))/sd(LIST_SIM_area_scaled[[i]][[j]]$mean_patch_area))
        }
      }
    }
  }
  
  # calculate species richness in each patch
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
  
  
  ## scale richness 
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
  
  
  ## scale evenness
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
  
  
  # subset IUCN species
  LIST_SIM_IUCN <- LIST_SIM
  
  for (i in 1: length(LIST_SIM_IUCN)) {
    for (j in 1: length(LIST_SIM_IUCN[[i]])) {
      for (k in 1: length(LIST_SIM_IUCN[[i]][[j]])){
        
        if (!is.na(LIST_SIM_IUCN[[i]][[j]])) {
          # subset only IUCN declining species
          LIST_SIM_IUCN[[i]][[j]][[k]] <- LIST_SIM_IUCN[[i]][[j]][[k]][,(names(LIST_SIM_IUCN[[i]][[j]][[k]]) %in% IUCN_declining), drop = FALSE]
        }
        
      }
    }
  }
  
  
  ## richness IUCN
  # calculate the species richness in each patch
  LIST_SIM_richness_IUCN <- LIST_SIM_IUCN
  
  for (i in 1: length(LIST_SIM_richness_IUCN)) {
    for (j in 1: length(LIST_SIM_richness_IUCN[[i]])) {
      for (k in 1: length(LIST_SIM_richness_IUCN[[i]][[j]])){
        
        if (!is.na(LIST_SIM_richness_IUCN[[i]][[j]])) {
          # take the column of each table from column two to column n-1, take their sum, and count how many times their colsum is larger than 0
          LIST_SIM_richness_IUCN[[i]][[j]][[k]] <- sum(colSums(LIST_SIM_richness_IUCN[[i]][[j]][[k]]) > 0)
        }
        
      }
    }
  }
  
  
  for (i in 1:length(LIST_SIM_richness_IUCN)){
    for ( j in 1: length(LIST_SIM_richness_IUCN[[i]])) {
      
      if (!is.na(LIST_SIM_richness_IUCN[[i]][[j]])) {
        LIST_SIM_richness_IUCN[[i]][[j]] <- do.call(rbind.data.frame, LIST_SIM_richness_IUCN[[i]][[j]])
        colnames(LIST_SIM_richness_IUCN[[i]][[j]]) <- "richness_IUCN"
      }
    }
  }
  
  
  ## evenness IUCN
  LIST_SIM_evenness_IUCN <- LIST_SIM_IUCN
  
  for (i in 1: length(LIST_SIM_evenness_IUCN)) {
    for (j in 1: length(LIST_SIM_evenness_IUCN[[i]])) {
      for (k in 1: length(LIST_SIM_evenness_IUCN[[i]][[j]])){
        
        if (!is.na(LIST_SIM_evenness_IUCN[[i]][[j]])) {
          # use the HurlbertPIE function on every table, excluding patch name (column 1) and patch size (column nrow)
          LIST_SIM_evenness_IUCN[[i]][[j]][[k]] <- paleotree::HurlbertPIE(colSums(LIST_SIM_evenness_IUCN[[i]][[j]][[k]]))
        } 
        
      }
    }
  }
  
  
  
  for (i in 1:length(LIST_SIM_evenness_IUCN)){
    for ( j in 1: length(LIST_SIM_evenness_IUCN[[i]])) {
      
      if (!is.na(LIST_SIM_evenness_IUCN[[i]][[j]])) {
        LIST_SIM_evenness_IUCN[[i]][[j]] <- do.call(rbind.data.frame, LIST_SIM_evenness_IUCN[[i]][[j]])
        colnames(LIST_SIM_evenness_IUCN[[i]][[j]]) <- "evenness_IUCN"
      }
    }
  }
  
  # scale richness and evenness IUCN
  LIST_SIM_richness_IUCN_scaled <- LIST_SIM_richness_IUCN
  
  for (i in 1:length(LIST_SIM_richness_IUCN_scaled)){
    for ( j in 1: length(LIST_SIM_richness_IUCN_scaled[[i]])) {
      
      if (!is.na(LIST_SIM_richness_IUCN_scaled[[i]][[j]])) {
        
        if (var(LIST_SIM_richness_IUCN_scaled[[i]][[j]]$richness_IUCN) == 0) {
          LIST_SIM_richness_IUCN_scaled[[i]][[j]]$richness_IUCN <- rep(0, length(LIST_SIM_richness_IUCN_scaled[[i]][[j]]$richness_IUCN))
        } else {
          LIST_SIM_richness_IUCN_scaled[[i]][[j]]$richness_IUCN <- ((LIST_SIM_richness_IUCN_scaled[[i]][[j]]$richness_IUCN - mean(LIST_SIM_richness_IUCN_scaled[[i]][[j]]$richness_IUCN))/sd(LIST_SIM_richness_IUCN_scaled[[i]][[j]]$richness_IUCN))
        }
      }
    }
  }
  
  
  
  LIST_SIM_evenness_IUCN_scaled <- LIST_SIM_evenness_IUCN
  
  for (i in 1:length(LIST_SIM_evenness_IUCN_scaled)){
    for ( j in 1: length(LIST_SIM_evenness_IUCN_scaled[[i]])) {
      
      if (!is.na(LIST_SIM_evenness_IUCN_scaled[[i]][[j]])) {
        
        if (var(LIST_SIM_evenness_IUCN_scaled[[i]][[j]]$evenness_IUCN) == 0) {
          LIST_SIM_evenness_IUCN_scaled[[i]][[j]]$evenness_IUCN <- rep(0, length(LIST_SIM_evenness_IUCN_scaled[[i]][[j]]$evenness_IUCN))
        } else {
          LIST_SIM_evenness_IUCN_scaled[[i]][[j]]$evenness_IUCN <- ((LIST_SIM_evenness_IUCN_scaled[[i]][[j]]$evenness_IUCN - mean(LIST_SIM_evenness_IUCN_scaled[[i]][[j]]$evenness_IUCN))/sd(LIST_SIM_evenness_IUCN_scaled[[i]][[j]]$evenness_IUCN))
        }
      }
    }
  }
  
  ## finalize
  LIST_SIM_combined <- LIST_SIM_area
  
  for (i in 1:length(LIST_SIM_combined)){
    for ( j in 1: length(LIST_SIM_combined[[i]])) {
      
      if (!is.na(LIST_SIM_combined[[i]][[j]])) {
        LIST_SIM_combined[[i]][[j]] <- cbind( LIST_SIM_area[[i]][[j]],
                                              LIST_SIM_area_scaled[[i]][[j]],
                                              LIST_SIM_area_sum[[i]][[j]],
                                              LIST_SIM_richness[[i]][[j]],  
                                              LIST_SIM_evenness[[i]][[j]],
                                              LIST_SIM_richness_scaled[[i]][[j]],
                                              LIST_SIM_evenness_scaled[[i]][[j]],
                                              LIST_SIM_richness_IUCN[[i]][[j]],
                                              LIST_SIM_evenness_IUCN[[i]][[j]],
                                              LIST_SIM_richness_IUCN_scaled[[i]][[j]],
                                              LIST_SIM_evenness_IUCN_scaled[[i]][[j]])
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
  final <- final[,-1]
  colnames(final) <- gsub(".1", "_scaled", colnames(final))
  final
  
} # end of the function
  
##########################################################
##########################################################
##########################################################

## FUNCTION 5
# CALCULATE BETA DIVERSITY METRICS ACROSS SETS OF PATCHES

beta_set <-function(object){
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
  
  
  
  for (i in 1: length(LIST_SIM)) {
    for (j in 1: length(LIST_SIM[[i]])) {
      for (k in 1: length(LIST_SIM[[i]][[j]])){
        
        if (length(LIST_SIM[[i]][[j]]) > 1) {
          
          LIST_SIM[[i]][[j]][[k]] <- LIST_SIM[[i]][[j]][[k]][,2:(ncol(LIST_SIM[[i]][[j]][[k]])-1)]
          
        } 
      }
    }
  }
  
  
  ## for each sets of patches, transform multiple rows into a single row based on their sum
  
  for (i in 1: length(LIST_SIM)) {
    for (j in 1: length(LIST_SIM[[i]])) {
      for (k in 1: length(LIST_SIM[[i]][[j]])){
        
        if (!is.na(LIST_SIM[[i]][[j]])) {
          LIST_SIM[[i]][[j]][[k]] <- as.data.frame(t(colSums(as.data.frame(LIST_SIM[[i]][[j]][k]))))
          
        }
      }
    }
  }
  
  # transform multiple sets of patches in a scenario into a single table
  
  for (i in 1: length(LIST_SIM)) {
    for (j in 1: length(LIST_SIM[[i]])) {
      
      if (!is.na(LIST_SIM[[i]][[j]])) {
        LIST_SIM[[i]][[j]] <- do.call(rbind, LIST_SIM[[i]][[j]])
        
      }
    }
  }
  
  
  
  # Ruzika beta div (abundance)
  LIST_SIM_RUZI <- LIST_SIM
  
  for (i in 1: length(LIST_SIM_RUZI)) {
    for (j in 1: length(LIST_SIM_RUZI[[i]])) {
      
      if (!is.na(LIST_SIM_RUZI[[i]][[j]])) {
        LIST_SIM_RUZI[[i]][[j]] <- beta.multi.abund(LIST_SIM[[i]][[j]], index.family="ruzicka")
      }
    }
  }
  
  
  # Jaccard beta diversity
  LIST_SIM_JAC <- LIST_SIM
  
  for (i in 1: length(LIST_SIM_JAC)) {
    for (j in 1: length(LIST_SIM_JAC[[i]])) {
      
      if (!is.na(LIST_SIM_JAC[[i]][[j]])) {
        LIST_SIM_JAC[[i]][[j]][LIST_SIM_JAC[[i]][[j]] > 0] <- 1 
        LIST_SIM_JAC[[i]][[j]] <- beta.multi(LIST_SIM_JAC[[i]][[j]], index.family="jaccard")    
        
      }
    }
  }
  
  
  ## convert into dataframes
  for (i in 1: length(LIST_SIM_RUZI)) {
    for (j in 1: length(LIST_SIM_RUZI[[i]])) {
      
      if (!is.na(LIST_SIM_RUZI[[i]][[j]])) {
        LIST_SIM_RUZI[[i]][[j]] <- as.data.frame(LIST_SIM_RUZI[[i]][[j]])
      }
    }
  }
  
  for (i in 1: length(LIST_SIM_JAC)) {
    for (j in 1: length(LIST_SIM_JAC[[i]])) {
      
      if (!is.na(LIST_SIM_JAC[[i]][[j]])) {
        LIST_SIM_JAC[[i]][[j]] <- as.data.frame(LIST_SIM_JAC[[i]][[j]])
      }
    }
  }
  
  ## add simulation and dataset number
  for (i in 1: length(LIST_SIM_RUZI)) {
    for (j in 1: length(LIST_SIM_RUZI[[i]])) {
      
      if (!is.na(LIST_SIM_RUZI[[i]][[j]])) {
        LIST_SIM_RUZI[[i]][[j]]$simulation <- i
        LIST_SIM_RUZI[[i]][[j]]$dataset <- levels(as.factor(data$dataset_label))[j]
        LIST_SIM_RUZI[[i]][[j]]$dataset_and_sim <- paste(LIST_SIM_RUZI[[i]][[j]]$dataset,LIST_SIM_RUZI[[i]][[j]]$simulation)
      }
    }
  }
  
  for (i in 1: length(LIST_SIM_JAC)) {
    for (j in 1: length(LIST_SIM_JAC[[i]])) {
      
      if (!is.na(LIST_SIM_JAC[[i]][[j]])) {
        LIST_SIM_JAC[[i]][[j]]$simulation <- i
        LIST_SIM_JAC[[i]][[j]]$dataset <- levels(as.factor(data$dataset_label))[j]
        LIST_SIM_JAC[[i]][[j]]$dataset_and_sim <- paste(LIST_SIM_JAC[[i]][[j]]$dataset,LIST_SIM_JAC[[i]][[j]]$simulation)
        
      }
    }
  }
  
  ## transform lists in tables
  for (i in 1: length(LIST_SIM_RUZI)) {
    LIST_SIM_RUZI[[i]] <- do.call(rbind, LIST_SIM_RUZI[[i]])
  }
  
  for (i in 1: length(LIST_SIM_JAC)) {
    LIST_SIM_JAC[[i]] <- do.call(rbind, LIST_SIM_JAC[[i]])
  }
  
  # merge all tables and remove NAs
  ruzi_table <- na.omit(do.call(rbind, LIST_SIM_RUZI))
  jac_table <- na.omit(do.call(rbind, LIST_SIM_JAC))
  
  table_merged <- merge(jac_table, ruzi_table)
  table_merged
}

##########################################################
##########################################################
##########################################################
## FUNCTION 6
# properties of sets of patches

properties_set <-function(object){ # object is one of the lists from SETS_PATCHES
  
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
          
          LIST_SIM_area[[i]][[j]][[k]] <- psych::describe(LIST_SIM_area[[i]][[j]][[k]][,ncol(LIST_SIM_area[[i]][[j]][[k]])])
          
        } 
        
      }
    }
  }
  # warnings appears because every list has more than one element, so the !is.na returns multiple TRUE statements when a list is not NA
  
  # transform lists into tabels
  for (i in 1:length(LIST_SIM_area)){
    for ( j in 1: length(LIST_SIM_area[[i]])) {
      
      if (!is.na(LIST_SIM_area[[i]][[j]])) {
        LIST_SIM_area[[i]][[j]] <- do.call(rbind.data.frame, LIST_SIM_area[[i]][[j]])
        #colnames(LIST_SIM_area[[i]][[j]]) <- "mean_patch_area"
      }
    }
  }
  
  
  
  ## finalize
  LIST_SIM_combined <- LIST_SIM_area
  
  for (i in 1:length(LIST_SIM_combined)){
    for ( j in 1: length(LIST_SIM_combined[[i]])) {
      

      if (!is.na(LIST_SIM_combined[[i]][[j]])) {
        LIST_SIM_combined[[i]][[j]] <- cbind( LIST_SIM_area[[i]][[j]] )
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
  LIST_SIM_combined <- LIST_SIM_combined[complete.cases(LIST_SIM_combined[,c("study","mean")]),]
  
  
  final <- LIST_SIM_combined
  #final <- final[,-1]
  #colnames(final) <- gsub(".1", "_scaled", colnames(final))
  final
  
  #final_twenty <- LIST_SIM_combined
} # end of the function


