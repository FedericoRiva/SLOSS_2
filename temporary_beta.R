library(betapart)

object <- list_simulations_twenty

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



# Ruzika beta div (abundance)
LIST_SIM_RUZI <- LIST_SIM

for (i in 1: length(LIST_SIM_RUZI)) {
  for (j in 1: length(LIST_SIM_RUZI[[i]])) {
    for (k in 1: length(LIST_SIM_RUZI[[i]][[j]])){
      
      if (!is.na(LIST_SIM_RUZI[[i]][[j]])) {
        # take the column of each table from column two to column n-1, take their sum, and count how many times their colsum is larger than 0
        LIST_SIM_RUZI[[i]][[j]][[k]] <- beta.multi.abund(LIST_SIM_RUZI[[i]][[j]][[k]][,2:(ncol(LIST_SIM_RUZI[[i]][[j]][[k]])-1)], index.family="ruzicka")

      }
      
    }
  }
}




# Jaccard beta diversity
LIST_SIM_JAC <- LIST_SIM

for (i in 1: length(LIST_SIM_JAC)) {
  for (j in 1: length(LIST_SIM_JAC[[i]])) {
    for (k in 1: length(LIST_SIM_JAC[[i]][[j]])){
      
      if (!is.na(LIST_SIM_JAC[[i]][[j]])) {
        # everything larger than 0 should be converted into a 1 for the species abundances
        LIST_SIM_JAC[[i]][[j]][[k]][,2:(ncol(LIST_SIM_JAC[[i]][[j]][[k]])-1)][LIST_SIM_JAC[[i]][[j]][[k]][,2:(ncol(LIST_SIM_JAC[[i]][[j]][[k]])-1)] > 0] <- 1 
        LIST_SIM_JAC[[i]][[j]][[k]] <- beta.multi(LIST_SIM_JAC[[i]][[j]][[k]][,2:(ncol(LIST_SIM_JAC[[i]][[j]][[k]])-1)], index.family="jaccard")
        
      }
    }
  }
}


#as.data.frame(t(do.call(rbind.data.frame, LIST_SIM_JAC[[1]][[5]][[1]])))

##
for (i in 1: length(LIST_SIM_RUZI)) {
  for (j in 1: length(LIST_SIM_RUZI[[i]])) {
    for (k in 1: length(LIST_SIM_RUZI[[i]][[j]])){
      
      if (!is.na(LIST_SIM_RUZI[[i]][[j]])) {
        LIST_SIM_RUZI[[i]][[j]][[k]] <- as.data.frame(t(do.call(rbind.data.frame, LIST_SIM_RUZI[[i]][[j]][[k]])))
      }
    }
  }
}









beta.multi.abund(LIST_SIM[[1]][[5]][[1]][,2:(ncol(LIST_SIM[[1]][[5]][[1]])-1)], index.family="ruzicka")

LIST_SIM[[1]][[5]][[1]][,2:(ncol(LIST_SIM[[1]][[5]][[1]])-1)]


# transform into presences and absences
LIST_SIM[[i]][[j]][[k]][,2:(ncol(LIST_SIM[[i]][[j]][[k]])-1)][LIST_SIM[[i]][[j]][[k]][,2:(ncol(LIST_SIM[[i]][[j]][[k]])-1)] > 0] <- 1 

