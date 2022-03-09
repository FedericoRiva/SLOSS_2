# Riva and Fahrig 
# Ecosystem decay across sets of patches
#
# prepare data
library(tidyr)
library(data.table)
library(dplyr)
library(purrr)
library(rlist)
library(splitstackshape)

# diversity estimation
library(iNEXT)
library(SpadeR)
library(vegan)
library(paleotree)
library(betapart)

# plotting
library(ggplot2)
library(ggthemes)
library(ggpubr)

# modeling
library(glmmTMB)

# remove scientific notation
options(scipen=999)
set.seed(654)

#################################################################
##### DATA PREPARATION
#################################################################

## open FragSAD dataset available at https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecy.2861
## .csv data and metadata provided for convenience
## open FragSAD dataset based on Chase et al. 2020 Nature
data = fread("C:\\Users\\feder\\OneDrive\\Desktop\\Riva Fahrig SLOSS #2\\data\\fragSAD_predicts_ewers.csv", header = TRUE)
metadata = fread("C:\\Users\\feder\\OneDrive\\Desktop\\Riva Fahrig SLOSS #2\\data\\new_meta_2_merge.csv", header = TRUE)

metadata$dataset_id <- as.factor(metadata$dataset_id)
data$dataset_label <- as.factor(data$dataset_label)

# open metadata from original paper, including effort
metadata_original = fread("C:\\Users\\feder\\OneDrive\\Desktop\\Riva Fahrig SLOSS #2\\data\\FragSAD_metadata_original.csv", header = TRUE)

# open table compiled by FR while re-assessing the papers
papers_reassessment = fread("C:\\Users\\feder\\OneDrive\\Desktop\\Riva Fahrig SLOSS #2\\data\\papers_reassessment.csv", header = TRUE)
papers_reassessment_filter <- papers_reassessment[,c(2,3)]

## data cleaning 
## remove 'continuous' patches for which we don't have real size
data <- subset(data, data$frag_size_char != 'continuous') 
## remove studies that I assessed on a one-by-one basis
data <- merge(data, papers_reassessment_filter, by = 'dataset_label')
data <- subset(data, data$included == 1) 

## remove patches with patch size assigned arbitrarily or with patch size larger than 90% of the total Habitat amount (8 patches); see "paper_reassessment"
## unless otherwise stated, removed only arbitrary sized large patches
## datasets where inspected by FR
data <- data[!(data$dataset_label == "Almeida-Gomes_2014" & data$frag_size_num > 7000),]  
data <- data[!(data$dataset_label == "Bell_2006_a" & data$frag_size_num > 1000),] 
data <- data[!(data$dataset_label == "Bell_2006_b" & data$frag_size_num > 1000),] 
data <- data[!(data$dataset_label == "Benedick_2006" & data$frag_size_num > 155000),] 
data <- data[!(data$dataset_label == "Berg_1997" & data$frag_size_num > 999),] 
data <- data[!(data$dataset_label == "Bernard_2007" & data$frag_size_num > 366),] 
data <- data[!(data$dataset_label == "Bolger_1997" & data$frag_size_num > 999),]  
data <- data[!(data$dataset_label == "Bragagnolo_2007" & data$frag_size_num > 999),] 
data <- data[!(data$dataset_label == "Cabrera-Guzman_2012" & data$frag_size_num > 500),] 
data <- data[!(data$dataset_label == "Davies_2003" & data$frag_size_num > 100),]  
data <- data[!(data$dataset_label == "Gavish_2012_a" & data$frag_size_num > 100),] 
data <- data[!(data$dataset_label == "Gavish_2012_b" & data$frag_size_num > 100),] 
data <- data[!(data$dataset_label == "Gavish_2012_c" & data$frag_size_num > 100),] 
data <- data[!(data$dataset_label == "Giladi_2011_a" & data$frag_size_num > 100),] 
data <- data[!(data$dataset_label == "Giladi_2011_c" & data$frag_size_num > 100),] 
data <- data[!(data$dataset_label == "Jung_2014" & data$frag_size_num > 999),] 
data <- data[!(data$dataset_label == "Manu_2007" & data$frag_size_num > 500),] 
data <- data[!(data$dataset_label == "Nogueira_2016" & data$frag_size_num > 200),] 
data <- data[!(data$dataset_label == "Struebig_2008" & data$frag_size_num > 15000),] 
data <- data[!(data$dataset_label == "Vasconcelos_2006" & data$frag_size_num > 999),] 
data <- data[!(data$dataset_label == "Savilaakso_2009" & data$frag_size_num > 250),] 

## remove site filtering column
data <- data[,-c(10)]

## label as factor
data$dataset_label <- as.factor(data$dataset_label)

## create a unique ID for study name and its patches
data$study_and_patch <- as.factor(paste(data$dataset_label, data$frag_id))

## in "patches", each row represent a patch in the original data
patches <- data[, c(1,2,10,5)]
patches<- unique(patches)

## in "assemblages" and "individuals", every row represents one taxon observed in a patch
assemblages <-  data[,c(1,2,10,8,9)] 
individuals <- assemblages[,-c(4)]

### create lists where each study is one element of a list and replace studies where 
### abundance data are reported as non-integers with that value multiplied by 100
list_individuals <- list()
for (i in 1 : 76){list_individuals[[i]] <- subset(individuals,dataset_label == levels(individuals$dataset_label)[i] )}
# Kappes_2009, Lambert_2003 and Dickman_1987 are expressed as density per area or trap; multiply by 100 to obtain integer values
for (i in 1 : 76){if (min(list_individuals[[i]]$abundance) < 1) {list_individuals[[i]]$abundance <- (list_individuals[[i]]$abundance * 100)}}  
individuals <- do.call(rbind.data.frame, list_individuals)

### number of individuals sampled in each study, in each patch, and total patch area sampled in each study
study_individuals <- aggregate(abundance ~ dataset_label, data=individuals, FUN=sum)
patch_individuals <- aggregate(abundance ~ study_and_patch, data=individuals, FUN=sum)
study_total_area <- aggregate(frag_size_num ~ dataset_label, data=patches, FUN=sum)

### add information to the dataset borrowing from original data
patch_individuals$dataset_label <- individuals$dataset_label[match(patch_individuals$study_and_patch,individuals$study_and_patch)]
patch_individuals$frag_id <- individuals$frag_id[match(patch_individuals$study_and_patch,individuals$study_and_patch)]
patch_individuals$total_abundance <- study_individuals$abundance[match(patch_individuals$dataset_label,study_individuals$dataset_label)]
patch_individuals$patch_area <- patches$frag_size_num[match(patch_individuals$study_and_patch,patches$study_and_patch)]
patch_individuals$total_patch_area <- study_total_area$frag_size_num[match(patch_individuals$dataset_label,study_total_area$dataset_label)]

#################################################################
##### CALCULATE SAMPLING EFFORT PER PATCH
#################################################################

# create a new column to split into list based on sample design
# "sample_id" was "plot_id" in FragSAD, defined as "plot_id: Plot-identifier of a sampling plot or a transect. If only one plot per site, 1 is given." in supplementary information
# note that in Owen_2008 the sample_id do not start from 1 in every patch
data$unique_sampling_event <- as.factor(paste(data$study_and_patch, data$sample_id))

# retain the effort (column 6), study_and_patch (column 10), and also the sampling event (column 11, created one line of code above) 
data_sampling_effort <- data[, c(6,10,11)]

# unique, because every sampling (column 11) event has the same effort (column 6)
data_sampling_effort <- unique(data_sampling_effort) # 1475 sampling events for 1222 patches, some have been sampled multiple times

# sum the total effort put in every patch
data_sampling_effort <- aggregate(sample_eff ~ study_and_patch, data = data_sampling_effort, FUN=sum) # by aggregating, we have a sum of the sampling events

# add a dataset_label from original data
data_sampling_effort$dataset_label <- data$dataset_label[match(data_sampling_effort$study_and_patch, data$study_and_patch)]

# transform sampling effort in values between 0 and 1, where 1 is the maximum sampling effort observed in a dataset
data_sampling_effort <- split(data_sampling_effort, data_sampling_effort$dataset_label)
for (i in 1: length(data_sampling_effort)) {
    data_sampling_effort[[i]]$sample_eff <- data_sampling_effort[[i]]$sample_eff/max(data_sampling_effort[[i]]$sample_eff)
}

data_sampling_effort <- do.call(rbind.data.frame, data_sampling_effort)

# match the relative sampling effort with information on every patch
patch_individuals$sampling_effort <- data_sampling_effort$sample_eff[match(patch_individuals$study_and_patch, data_sampling_effort$study_and_patch)]
# the number of individuals expected if all patches were sampled with the same effort
patch_individuals$abundance_per_effort <- patch_individuals$abundance / patch_individuals$sampling_effort


# calculate normalized # of individuals adjusted (i.e., at equal sampling effort) per patch 
patch_individuals <- split(patch_individuals, patch_individuals$dataset_label)
for (i in 1:length(patch_individuals)) {
  patch_individuals[[i]]$normalized_individuals_per_sampl_unit <- patch_individuals[[i]]$abundance_per_effort / max(patch_individuals[[i]]$abundance_per_effort)
}

patch_individuals <- do.call(rbind.data.frame, patch_individuals)

# normalized_individuals_per_sampl_unit2 is a transformation to avoind values equal to 1 or zero; necessary to fit beta regression
# transformation follows Smithson, M., & Verkuilen, J. (2006). A better lemon squeezer? Maximum-likelihood regression with beta-distributed dependent variables. Psychological Methods, 11(1), 54–71.
patch_individuals$normalized_individuals_per_sampl_unit2 <- (patch_individuals$normalized_individuals_per_sampl_unit * (nrow(patch_individuals)-1) + 0.5) / nrow(patch_individuals)
patch_individuals$log_patch_size <- log(patch_individuals$patch_area)

# c.lfs + (c.lfs | dataset_label)  we followed Chase et al. 2020 model structure
model <- glmmTMB(normalized_individuals_per_sampl_unit2 ~ log_patch_size  +  
                    (log_patch_size|dataset_label), 
                  data = patch_individuals, family = "gaussian")## family=beta_family(link="logit")) ## 
                  # gaussian distributio ndoes a better job than beta distribution 
summary(model)
AIC(model)


patch_individuals$predicted_normalized_individuals <- predict(model, type = "response")

# strong effect of dataset random effect
plot(patch_individuals$predicted_normalized_individuals ~ patch_individuals$dataset_label)

# some datasets have a positive area effect, some have a negative patch area effect
plot(predict(model, type = "response"), patch_individuals$log_patch_size)

# relationship between model predictions and observed number of individuals per sampling unit (normalized to a value between 0 and 1)
plot(predict(model, type = "response"), patch_individuals$normalized_individuals_per_sampl_unit2)
abline(0, 1, col = "red")
#################################################################
##### CALCULATE INDIVIDUALS PER RANDOMIZATIONS
#################################################################

# calculate the proportion of individuals per patch
patch_individuals$abundance_per_ha <- patch_individuals$abundance / patch_individuals$patch_area
patch_individuals$relative_patch_size <- patch_individuals$patch_area / patch_individuals$total_patch_area 
# calculate the minimum to define the # of individuals simulated 
min_abundance_per_ha_study <- aggregate(abundance_per_ha ~ dataset_label, data=patch_individuals, FUN=min)
colnames(min_abundance_per_ha_study)[2] <- "min_abundance_per_ha_study"

# merge with patch-level information
patch_individuals <- merge(patch_individuals, min_abundance_per_ha_study, by = 'dataset_label')
patch_individuals$simulation_constant_abundance <- patch_individuals$min_abundance_per_ha_study * patch_individuals$patch_area
#patch_individuals$simulations_adjusted_abundance <- patch_individuals$simulation_constant_abundance * patch_individuals$predicted_normalized_individuals

# increase the number of individuals to be resampled in each dataset by standardizing the predictions within each dataset
patch_individuals <- split(patch_individuals, patch_individuals$dataset_label)
for (i in 1:length(patch_individuals)) {
  patch_individuals[[i]]$predicted_normalized_individuals_adjusted <- patch_individuals[[i]]$predicted_normalized_individuals / max(patch_individuals[[i]]$predicted_normalized_individuals)
}
patch_individuals <- do.call(rbind.data.frame, patch_individuals)

patch_individuals$simulations_adjusted_abundance <- patch_individuals$simulation_constant_abundance * patch_individuals$predicted_normalized_individuals_adjusted



#################################################################
##### SIMULATE 100 TIMES THE NUMBER OF INDIVIDUALS TO BE RESAMPLED IN EACH PATCH
#################################################################

# create lists where each study is one element of a list
list_patches <- list()
for (i in 1 : 76){list_patches[[i]] <- subset(patches,dataset_label == levels(data$dataset_label)[i] )}

list_assemblages <- list()
for (i in 1 : 76){list_assemblages[[i]] <- subset(assemblages,dataset_label == levels(data$dataset_label)[i] )}

# if the smallest value of a list is < 1, multiply all the abundances of the species recorded as relative occurrences for 100 
# so that rows can be created proportionally to their incidence
for (i in 1 : 76){if (min(list_assemblages[[i]]$abundance) < 1) {list_assemblages[[i]]$abundance <- (list_assemblages[[i]]$abundance * 100)}}  

# extend so that each row has nrow = nindividuals
for (i in 1 : 76) {list_assemblages[[i]] <- splitstackshape::expandRows(list_assemblages[[i]], "abundance") }
list_assemblages2 <- rbindlist(list_assemblages)

# list of all species seen in all patches in the original datasets
list_assemblages_patches <- list()
for (i in 1 : nlevels(list_assemblages2$study_and_patch)) {list_assemblages_patches[[i]] <-  subset(list_assemblages2,study_and_patch == levels(list_assemblages2$study_and_patch)[i] )}


### create 100 number of individuals sampled per patch
list_sim_individuals <- list()

for (i in 1 : nrow(patch_individuals)){
  list_sim_individuals[[i]] <- replicate(100, (floor(patch_individuals$simulations_adjusted_abundance[[i]]) + rbinom(1, 1,(patch_individuals$simulations_adjusted_abundance[[i]] - floor(patch_individuals$simulations_adjusted_abundance[[i]])))))
}
## these are the numbers of individuals that will be selected in each patch; 100 randomizations



################################################################
####### create randomized tables
################################################################

sim_ind <- do.call(rbind.data.frame, list_sim_individuals)
hist(log10(rowSums(sim_ind)))

## TEST IF SIMULATIONS ARE CORRECT: check if simulations are bigger than the number of individuals observed in a patch
# list_test_sim <- rep( list(list()), 1186 )
# for (i in 1:100){
#   for (k in 1: 1186){
#     list_test_sim[[k]][i] <- isTRUE(sim_ind[k,i] > patch_individuals[k, 3])
#  }
# }
# 
# test_sim <- do.call(rbind.data.frame, list_test_sim)
# any(test_sim==TRUE)
## none of the simulated # of individuals in each patch is larger than the original number of individuals observed


### SAMPLE ROWS ON A PER-PATCH BASIS TO HAVE # OF INDIVIDUALS PROPORTIONAL TO PATCH SIZE
get_rows <- function(df, rows) df[rows, , drop = FALSE]

# split_data <- split(list_assemblages2, list_assemblages2$study_and_patch)
# group_sizes <- vapply(split_data, nrow, integer(1))
# n <-sim_ind[,3]
# sampled_obs <- mapply(sample, group_sizes, n)
# sampled_data <- mapply(get_rows, split_data, sampled_obs, SIMPLIFY = FALSE)
# sampled_data2 <- do.call(rbind, sampled_data)


# followed https://jennybc.github.io/purrr-tutorial/ls12_different-sized-samples.html
list_sampled_data <- list()
split_data <- split(list_assemblages2[,1:4], list_assemblages2$study_and_patch)
group_sizes <- vapply(split_data, nrow, integer(1))

for (i in 1:100){
  n <-sim_ind[,i]
  sampled_obs <- mapply(sample, group_sizes, n)
  sampled_data <- mapply(get_rows, split_data, sampled_obs, SIMPLIFY = FALSE)
  sampled_data2 <- do.call(rbind, sampled_data) 
  
  list_sampled_data[[i]] <- sampled_data2
}


# ## if needed, break in lists
for (i in 1:100){
  list_sampled_data[[i]] <- split(list_sampled_data[[i]], list_sampled_data[[i]][,3]) #
}





# convert from a series of rows each representing a sampled individuals, into a vectors of abundances per patch
# when a patch has zero simulated individuals, this still create a column that I will remove later to avoid losing the patch from the dataset
for (i in 1:length(list_sampled_data)){
  for (j in 1:length(list_sampled_data[[i]])){
    if(nrow(list_sampled_data[[i]][j][[1]]) == 0) {
      mydf <- (data.frame(Var1 = "no_species_present", Freq = 0))
      list_sampled_data[[i]][j][[1]]<- setNames(data.frame(t(mydf[,-1])), mydf[,1])
    } else {
      list_sampled_data[[i]][j][[1]] <- setNames(data.frame(t(as.data.frame(table(list_sampled_data[[i]][j][[1]]$species))[,-1])), as.data.frame(table(list_sampled_data[[i]][j][[1]]$species))[,1])
    }
  }
}



# split the 100 lists into 76 datasets
study_index <- patches$dataset_label
for (i in 1: length(list_sampled_data)){
  list_sampled_data[[i]] <- split(list_sampled_data[[i]], study_index)
}

# function to convert list of vectors of species abundances into table
to_table <- function(list_to_table) {
  table_prova <- rbindlist(lapply(list_to_table, function(x) as.data.frame.list(x)), fill=TRUE)
  table_prova[is.na(table_prova)] <- 0
  table_prova$study_and_patch <- names(list_to_table)
  list_to_table <- table_prova
}

# make each list (100 simulations per 76 studies) into a table of species per sites
for (i in 1:length(list_sampled_data)){
  for (j in 1: length(list_sampled_data[[i]])){
    list_sampled_data[[i]][[j]] <- to_table(list_sampled_data[[i]][[j]]) #rbindlist(lapply(list_sampled_data[[i]][[j]], function(x) as.data.frame.list(x)), fill=TRUE)
  }
}

# remove empty column that was necessary to retain patches with zero simulated individuals
for (i in 1:length(list_sampled_data)){
  for (j in 1: length(list_sampled_data[[i]])){
    
    if("no_species_present" %in% colnames(list_sampled_data[[i]][j][[1]])){
      list_sampled_data[[i]][j][[1]] <- select(list_sampled_data[[i]][j][[1]], -no_species_present)
    } 
    
  }
}

# TEST: did we remove the "no_species_present" column added to retain empty patches? Note that empty patches must be retained for the simulations 
# list_sampled_data[[1]][74][[1]] # does not have the "no_species_present" column anymore

# add patch size to each list
list_patches_unlisted <- do.call(rbind.data.frame, list_patches)
list_patches_unlisted <- list_patches_unlisted[,3:4]

for (i in 1:length(list_sampled_data)){
  for (j in 1: length(list_sampled_data[[i]])){
    
    list_sampled_data[[i]][j][[1]] <- merge(list_sampled_data[[i]][j][[1]], list_patches_unlisted, by = "study_and_patch")
    
  }
}


####################################################################################
######### SIMULATE SETS OF PATCHES OF EQUAL HABITAT AMOUNTS
####################################################################################



SETS_PATCHES <- function(db, percent) {  # db= dataset, percent = habitat amount
  
  set.seed(111)
  # if (nrow(db) < 10) {
  #   n_db <- 100 # 200 simulations when the number of patches in a dataset is < 10
  # } else {
  #   n_db <- 100 # 200 simulations when the number of patches in a dataset is > 10
  # }
  n_db <- 100
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
  #list_db <- rlist::list.filter(list_db, sum(frag_size_num) > (target_area - target_area*0.05) )
  
  # list_area <- list()
  # for (i in 1 : length(list_db)) {list_area[i] <- sum(list_db[[i]]$frag_size_num)}
  # list_area <- do.call(rbind.data.frame, list_area)
  
  for (i in 1: length(list_db)) {
    list_db[[i]] <- as.data.frame(list_db[[i]]$study_and_patch)
    colnames(list_db[[i]]) <- "study_and_patch"
  }
  
  print(list_db)
  # list_area
}

#SETS_PATCHES(list_patches[[1]], 0.3)



#
#tryCatch puts NAs when error in the function (i.e., when the function cannot find sets of patches)
# list_simulations_ten <- vector("list",length(list_patches))
# for (i in 1:length(list_patches)) {list_simulations_ten[[i]] <- tryCatch(SETS_PATCHES(list_patches[[i]], 0.1), error=function(e) print(NA))}

list_simulations_twenty <- vector("list",length(list_patches))
for (i in 1:length(list_patches)) {list_simulations_twenty[[i]] <- tryCatch(SETS_PATCHES(list_patches[[i]], 0.2), error=function(e) print(NA))}

list_simulations_thirty <- vector("list",length(list_patches))
for (i in 1:length(list_patches)) {list_simulations_thirty[[i]] <- tryCatch(SETS_PATCHES(list_patches[[i]], 0.3), error=function(e) print(NA))}

list_simulations_forty <- vector("list",length(list_patches))
for (i in 1:length(list_patches)) {list_simulations_forty[[i]] <- tryCatch(SETS_PATCHES(list_patches[[i]], 0.4), error=function(e) print(NA))}
 
list_simulations_fifty <- vector("list",length(list_patches))
for (i in 1:length(list_patches)) {list_simulations_fifty[[i]] <- tryCatch(SETS_PATCHES(list_patches[[i]], 0.5), error=function(e) print(NA))}

list_simulations_sixty <- vector("list",length(list_patches))
for (i in 1:length(list_patches)) {list_simulations_sixty[[i]] <- tryCatch(SETS_PATCHES(list_patches[[i]], 0.6), error=function(e) print(NA))}
 
list_simulations_seventy <- vector("list",length(list_patches))
for (i in 1:length(list_patches)) {list_simulations_seventy[[i]] <- tryCatch(SETS_PATCHES(list_patches[[i]], 0.7), error=function(e) print(NA))}

list_simulations_eighty <- vector("list",length(list_patches))
for (i in 1:length(list_patches)) {list_simulations_eighty[[i]] <- tryCatch(SETS_PATCHES(list_patches[[i]], 0.8), error=function(e) print(NA))}

list_simulations_ninety <- vector("list",length(list_patches))
for (i in 1:length(list_patches)) {list_simulations_ninety[[i]] <- tryCatch(SETS_PATCHES(list_patches[[i]], 0.9), error=function(e) print(NA))}


# list_simulations <- list(#list_simulations_ten,
#                          list_simulations_twenty,
#                          #list_simulations_thirty,
#                          list_simulations_forty,
#                          #list_simulations_fifty,
#                          list_simulations_sixty,
#                          #list_simulations_seventy,
#                          list_simulations_eighty) #,
#                          #list_simulations_ninety)




########################################################################################################
########################################################################################################
# open IUCN declining species
IUCN_declining = fread("C:\\Users\\feder\\OneDrive\\Desktop\\Riva Fahrig SLOSS #2\\data\\IUCN_simple_summary.csv", header = TRUE)
IUCN_declining <- subset(IUCN_declining, scientificName %in% data$species)

table((as.factor(IUCN_declining$className)))

IUCN_declining <- IUCN_declining$scientificName
IUCN_declining <- gsub(" ", ".", IUCN_declining)

########################################################################################################
########################################################################################################


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
  
  #final_twenty <- LIST_SIM_combined
} # end of the function




# # Optionally set colours using RColorBrewer
# library(RColorBrewer)
# cols = brewer.pal(4, "Blues")
# # Define colour pallete
# pal = colorRampPalette(c("blue", "red"))
# # Use the following line with RColorBrewer
# pal = colorRampPalette(cols)
# # Rank variable for colour assignment
# final$order = findInterval(final$mean_patch_area_sc, sort(final$mean_patch_area_sc))
# # Make plot
# plot(richness ~ total_patch_area, final, pch=19, col=pal(nrow(final))[final$order])





#alpha
final_twenty <- alpha_set(list_simulations_twenty)
final_thirty <- alpha_set(list_simulations_thirty)
final_forty <- alpha_set(list_simulations_forty)
final_fifty <- alpha_set(list_simulations_fifty)
final_sixty <- alpha_set(list_simulations_sixty)
final_seventy <- alpha_set(list_simulations_seventy)
final_eighty <- alpha_set(list_simulations_eighty)

#_beta
final_twenty_beta <- beta_set(list_simulations_twenty)
final_thirty_beta <- beta_set(list_simulations_thirty)
final_forty_beta <- beta_set(list_simulations_forty)
final_fifty_beta <- beta_set(list_simulations_fifty)
final_sixty_beta <- beta_set(list_simulations_sixty)
final_seventy_beta <- beta_set(list_simulations_seventy)
final_eighty_beta <- beta_set(list_simulations_eighty)

# add scenario
final_twenty$habitat_amount <- rep("twenty_percent", nrow(final_twenty))
final_thirty$habitat_amount <- rep("thirty_percent", nrow(final_thirty))
final_forty$habitat_amount <- rep("forty_percent", nrow(final_forty))
final_fifty$habitat_amount <- rep("fifty_percent", nrow(final_fifty))
final_sixty$habitat_amount <- rep("sixty_percent", nrow(final_sixty))
final_seventy$habitat_amount <- rep("seventy_percent", nrow(final_seventy))
final_eighty$habitat_amount <- rep("eighty_percent", nrow(final_eighty))

final_twenty_beta$habitat_amount <- rep("twenty_percent", nrow(final_twenty_beta))
final_thirty_beta$habitat_amount <- rep("thirty_percent", nrow(final_thirty_beta))
final_forty_beta$habitat_amount <- rep("forty_percent", nrow(final_forty_beta))
final_fifty_beta$habitat_amount <- rep("fifty_percent", nrow(final_fifty_beta))
final_sixty_beta$habitat_amount <- rep("sixty_percent", nrow(final_sixty_beta))
final_seventy_beta$habitat_amount <- rep("seventy_percent", nrow(final_seventy_beta))
final_eighty_beta$habitat_amount <- rep("eighty_percent", nrow(final_eighty_beta))


# bring together
table_analysis <- rbind(final_twenty, final_thirty,
                        final_forty, final_fifty,
                        final_sixty, final_seventy,
                        final_eighty)

colnames(table_analysis)[13] <- "dataset_id"
table_analysis <- merge(table_analysis, metadata, by = "dataset_id")

write.csv(table_analysis, "table_analysis.csv")
##
table_analysis_beta <- rbind(final_twenty_beta, final_thirty_beta,
                        final_forty_beta, final_fifty_beta,
                        final_sixty_beta, final_seventy_beta,
                        final_eighty_beta)

write.csv(table_analysis_beta, "table_analysis_beta.csv")

##
library(glmmTMB)
library(effects)
library(ggeffects)

# table_analysis$log10ric <- log10(table_analysis$richness)
# table_analysis$log10_mean_patch_area <- log10(table_analysis$mean_patch_area)

model <- glmmTMB(evenness_scaled ~ mean_patch_area_sc + #* habitat_amount + 
                   (1|habitat_amount) + 
                   (1|dataset_id/patch_set_number) + #patch set number is the simulated assemblage, from 1 to 100
                    (1|simulation_number), 
                  data = table_analysis, family = "gaussian")

table_analysis$log_rich <- log(table_analysis$richness + 1)
table_analysis$log_area <- log(table_analysis$total_patch_area)

model <- glmmTMB(log_rich ~ log_area *
                  mean_patch_area_sc + #* habitat_amount + 
                   (1|habitat_amount) + 
                   (1|dataset_id/patch_set_number) + #patch set number is the simulated assemblage, from 1 to 100
                   (1|simulation_number), 
                 data = table_analysis, family = "gaussian")

summary(model)

table_analysis$predicted_richness <- predict(model, type = "response")

plot(allEffects(model), type = "response")





plot_model <- ggpredict(model, 
                        c("log_area [all]", "mean_patch_area_sc [-1, 0, 1]" ),
                        type = "re"
) 

plot(plot_model) + 
  labs(
    x = "log (area)", 
    y = "log (richness + 1)", 
    title = ""  ) + 
  labs(colour = "SD from mean patch size in a scenario")




# check ggplot tutorial here https://aosmith.rbind.io/2018/11/16/plot-fitted-lines/

ggplot(table_analysis, aes(x=log10_mean_patch_area, y=log10ric, color=dataset_id, shape=habitat_amount)) +
  #geom_point() +
  geom_smooth(method = "lm")+
  #geom_line(aes(y = predicted_richness))+
  theme(legend.position="none")

# table_analysis$evenness2 <- (table_analysis$evenness * (nrow(table_analysis)-1) + 0.5) / nrow(table_analysis)
# 
# model <- glmmTMB(evenness2 ~ mean_patch_area * habitat_amount + taxa +
#                    (1|dataset_id/patch_set_number) +
#                    (1|simulation_number), 
#                  data = table_analysis, family=beta_family(link="logit"))
# 
# summary(model)
# AIC(model)


plot_model <- ggpredict(model, c("mean_patch_area", "taxa"), type = "re") 
plot(plot_model)






########################################################################################################
########################################################################################################
################################################################################################################################################################################################################
########################################################################################################
################################################################################################################################################################################################################
########################################################################################################
################################################################################################################################################################################################################
########################################################################################################
################################################################################################################################################################################################################
########################################################################################################
################################################################################################################################################################################################################
########################################################################################################
########################################################################################################
# 
# 
# ## ELE analysis
# #### sloss function from Lexiguel
# sloss <- function(table, env=data.frame(), area) {
#   if(!is.matrix(table)) table <- as.matrix(table)
#   area <- substitute(area)
#   area <- eval(area, env, parent.frame())
#   SLOSS <- list(SL=list(), LS=list())
#   # First the calculation from small to large
#   SLOSS$SL$area <- c(0, cumsum(area[order(area)]))
#   Flor <- apply(table[order(area),], 2, cumsum)
#   Flor[Flor > 0] <- 1
#   SLOSS$SL$species <- c(0, apply(Flor, 1, sum))
#   # Now the calculation from large to small
#   SLOSS$LS$area <- c(0, cumsum(area[order(area, decreasing=TRUE)]))
#   Flor <- apply(table[order(area, decreasing=TRUE),], 2, cumsum)
#   Flor[Flor > 0] <- 1
#   SLOSS$LS$species <- c(0, apply(Flor, 1, sum))
#   # Calculation of SLOSS index
#   SLOSS$Index <- with(SLOSS$SL, curve_area(area,
#                                            species))/with(SLOSS$LS, curve_area(area, species))
#   # Final object
#   class(SLOSS) <- c("SLOSS","list")
#   return(SLOSS)
# }
# 
# curve_area <- function(x, y, bottom=0) {
#   D1 <- c(diff(x))
#   D2 <- c(diff(y))
#   Area <- sum(D1*((y - bottom)[-length(y)]) + D1*D2/2, na.rm=TRUE)
#   return(Area)
# }
# 
# 
# # edited from original function in Lexiguel to include points
# plot.SLOSS <- function(x, y=NULL, sl.lty=2, sl.lwd=1, sl.col="black", ls.lty=1,
#                        ls.lwd=1, ls.col="black", show.index=TRUE, digits.index=2, cex.index=1,
#                        pos.index=c(0.05,0.95), show.legend=FALSE, pos.legend="bottomright",
#                        bty.legend="o", main="SLOSS curves",...) {
#   with(x$SL, plot(area, species, type="l", lty=sl.lty, lwd=sl.lwd, col=sl.col,
#                   main=main, ...))
#   with(x$LS, lines(area, species, lty=ls.lty, lwd=ls.lwd, col=ls.col))
#   
#   with(x$SL, points(area, species, pch = 18, cex = 1))
#   with(x$LS, points(area, species, pch = 18, cex = 1))
#   
#   if(show.legend) {
#     legend(pos.legend, lty=c(sl.lty,ls.lty), lwd=c(sl.lwd,ls.lwd),
#            legend=c("small to large","large to small"), bty=bty.legend)
#   }
#   if(show.index) {
#     with(x$SL, text(max(area)*pos.index[1], max(species)*pos.index[2],
#                     labels=paste("SLOSS-index =",
#                                  round(x$Index, digits.index)),
#                     cex=cex.index, pos=4))
#   }
# }
# 
# 
# 
# 
# # 
# INPUT <- list_sampled_data[[10]][74][[1]]
# 
# SL_OR_SS <- function(INPUT) {
#   
#   table_sloss <- as.matrix(INPUT)
#   area_sloss <-   table_sloss[,ncol(table_sloss)]
#   
#   table_sloss <- table_sloss[, - c(1, ncol(table_sloss))]
#   sloss_comparison <- sloss(table_sloss, area = area_sloss)
#   
#   # sl <- approx(sloss_comparison$SL$area, 
#   #              sloss_comparison$SL$species, 
#   #              xout = seq(from = as.numeric(min(area_sloss)), to = (max(sloss_comparison$SL$area) - as.numeric(max(area_sloss))), length.out = 100), 
#   #              method = "linear")
#   # 
#   # ls <- approx(sloss_comparison$LS$area, 
#   #              sloss_comparison$LS$species, 
#   #              xout = seq(from = as.numeric(min(area_sloss)), to = (max(sloss_comparison$SL$area) - as.numeric(max(area_sloss))), length.out = 100), 
#   #              method = "linear")
#   
#   sl <- approx(sloss_comparison$SL$area,
#                sloss_comparison$SL$species,
#                xout = seq(from = as.numeric(max(area_sloss)), to = (max(sloss_comparison$SL$area) - as.numeric(max(area_sloss))), length.out = 100),
#                method = "linear")
#   
#   ls <- approx(sloss_comparison$LS$area,
#                sloss_comparison$LS$species,
#                xout = seq(from = as.numeric(max(area_sloss)), to = (max(sloss_comparison$SL$area) - as.numeric(max(area_sloss))), length.out = 100),
#                method = "linear")
#   
#   plot.SLOSS(sloss_comparison, show.index = F)
#   points(ls, col = "blue")
#   points(sl, col = "red")
#   
#   comparison_sloss <- as.data.frame(t(rbind(ls$y, sl$y)))
#   
#   comparison_sloss$outcome <- comparison_sloss$V1 - comparison_sloss$V2 # V1 is LS
#   
#   outcome <- if (all(comparison_sloss$outcome > 0)){
#     print("SL > SS")
#   } else {
#     if (all(comparison_sloss$outcome < 0)){print("SS > SL")
#     } else {
#       print("SS = SL")
#     }
#   }
#   
#   #return(outcome)
# }
# 
# SL_OR_SS(list_sampled_data[[88]][7][[1]]) # first number simulations from 1 to 100, second number dataset from 1 to 75
# 
# ## ELE study
# 
# 
# ELE_data <- list_sampled_data
# 
# 
# for (i in 1:length(ELE_data)){
#   for (j in 1: length(ELE_data[[i]])){
#     
#     ELE_data[[i]][j][[1]] <- SL_OR_SS(ELE_data[[i]][j][[1]])
#     
#   }
# }
# 
# study_names <- names(ELE_data[[1]])
# 
# 
# for (i in 1:length(ELE_data)){
#   ELE_data[[i]] <- do.call(rbind.data.frame, ELE_data[[i]])
# }
# 
# for (i in 1:length(ELE_data)){
#   colnames(ELE_data[[i]]) <- "comparison"
# }
# 
# ele_data_analysis <- do.call(rbind.data.frame, ELE_data)
# ele_data_analysis$dataset_id <- rep(study_names, 100)
# 
# 
# ele_data_analysis <- merge(ele_data_analysis, metadata, by = "dataset_id")
# table(ele_data_analysis$comparison)
# table(ele_data_analysis$taxa)
# 
# # only 45 studyes give SL > SS
# 
# ele_data_analysis$SS <- as.factor(ele_data_analysis$comparison)
# levels(ele_data_analysis$SS)
# 
# levels(ele_data_analysis$SS)[match("SL > SS",levels(ele_data_analysis$SS))] <- "0"
# levels(ele_data_analysis$SS)[match("SS = SL",levels(ele_data_analysis$SS))] <- "0"
# levels(ele_data_analysis$SS)[match("SS > SL",levels(ele_data_analysis$SS))] <- "1"
# 
# 
# 
# 
# ele_data_analysis$SL <- as.factor(ele_data_analysis$comparison)
# levels(ele_data_analysis$SL)
# 
# levels(ele_data_analysis$SL)[match("SL > SS",levels(ele_data_analysis$SL))] <- "1"
# levels(ele_data_analysis$SL)[match("SS = SL",levels(ele_data_analysis$SL))] <- "0"
# levels(ele_data_analysis$SL)[match("SS > SL",levels(ele_data_analysis$SL))] <- "0"
# 
# 
# 
# ele_data_analysis$dataset_id <- as.factor(ele_data_analysis$dataset_id)
# ele_data_analysis$sphere.fragment <- as.factor(ele_data_analysis$sphere.fragment)
# ele_data_analysis$sphere.matrix <- as.factor(ele_data_analysis$sphere.matrix)
# ele_data_analysis$biome <- as.factor(ele_data_analysis$biome)
# ele_data_analysis$taxa <- as.factor(ele_data_analysis$taxa)
# ele_data_analysis$time.since.fragmentation <- as.factor(ele_data_analysis$time.since.fragmentation)
# ele_data_analysis$Matrix.category <- as.factor(ele_data_analysis$Matrix.category)
# 
# 
# # calculate patch size evenness for every study
# evenness <- list()
# for (i in 1: length(list_patches)){
#   evenness[[i]] <- (vegan::diversity(list_patches[[i]]$frag_size_num))/log(length(list_patches[[i]]$frag_size_num))
# }
# 
# evenness <- do.call(rbind.data.frame, evenness)
# evenness$dataset_id <- study_names
# colnames(evenness) <- c("evenness", "dataset_id")
# 
# ele_data_analysis <- merge(ele_data_analysis, evenness, by = "dataset_id")
# 
# # calculate species richness for every study
# richness <- list()
# for (i in 1: length(list_patches)){
#   richness[[i]] <- length(table(list_assemblages[[i]]$species))
# }
# 
# richness <- do.call(rbind.data.frame, richness)
# richness$dataset_id <- study_names
# colnames(richness) <- c("richness", "dataset_id")
# 
# ele_data_analysis <- merge(ele_data_analysis, richness, by = "dataset_id")
# 
# # add ectothermic vs endothermic
# ele_data_analysis$temperature <- ele_data_analysis$taxa
# 
# levels(ele_data_analysis$temperature)[match("amphibians & reptiles",levels(ele_data_analysis$temperature))] <- "ectotherm"
# levels(ele_data_analysis$temperature)[match("birds",levels(ele_data_analysis$temperature))] <- "endotherm"
# levels(ele_data_analysis$temperature)[match("invertebrates",levels(ele_data_analysis$temperature))] <- "ectotherm"
# levels(ele_data_analysis$temperature)[match("mammals",levels(ele_data_analysis$temperature))] <- "endotherm"
# levels(ele_data_analysis$temperature)[match("plants",levels(ele_data_analysis$temperature))] <- "ectotherm"
# 
# 
# # add verts vs inv
# ele_data_analysis$vert <- ele_data_analysis$taxa
# 
# levels(ele_data_analysis$vert)[match("amphibians & reptiles",levels(ele_data_analysis$vert))] <- "vert"
# levels(ele_data_analysis$vert)[match("birds",levels(ele_data_analysis$vert))] <- "vert"
# levels(ele_data_analysis$vert)[match("invertebrates",levels(ele_data_analysis$vert))] <- "invert"
# levels(ele_data_analysis$vert)[match("mammals",levels(ele_data_analysis$vert))] <- "vert"
# levels(ele_data_analysis$vert)[match("plants",levels(ele_data_analysis$vert))] <- "plant"
# 
# # add whether a taxon flies or not
# ele_data_analysis$group <- c("bees", "bees", "lizards","orthoptera", "lepidoptera",
#                              "frogs","lizards","lepidoptera", "birds", "bats", 
#                              "mammals", "spiders", "bees", "bees", "amphibians",
#                              "plants", "plants", "birds", "ants", "termites",
#                              "mammals", "birds", "birds", "mammals", "beetles",
#                              "beetles","mammals", "spiders", "spiders","spiders",
#                              "mammals", "plants", "plants", "plants", "birds",
#                              "bats", "spiders", "bees", "beetles", "birds",
#                              "spiders", "snails", "beetles", "spiders", "mammals",
#                              "amphibians","amphibians","reptiles","birds","birds",
#                              "birds", "bats", "beetles", "bees", "spiders",
#                              "orthoptera","beetles", "beetles","frogs", "snails",
#                              "plants","bees", "plants","plants","bees",
#                              "lepidoptera","mammals","lepidoptera","bats","birds",
#                              "lepidoptera", "ants", "beetles", "reptiles", "lepidoptera")
# 
# 
# ele_data_analysis$fly <- c("fly", "fly", "no_fly","fly", "fly",
#                            "no_fly","no_fly","fly", "fly", "fly", 
#                            "no_fly", "no_fly", "fly", "fly", "no_fly",
#                            "plants", "plants", "fly", "fly", "fly",
#                            "no_fly", "fly", "fly", "no_fly", "fly",
#                            "fly","no_fly", "no_fly", "no_fly","no_fly",
#                            "no_fly", "plants", "plants", "plants", "fly",
#                            "fly", "no_fly", "fly", "fly", "fly",
#                            "no_fly", "no_fly", "fly", "no_fly", "no_fly",
#                            "no_fly","no_fly","no_fly","fly","fly",
#                            "fly", "fly", "fly", "fly", "no_fly",
#                            "fly","fly", "fly","no_fly", "no_fly",
#                            "plants","fly", "plants","plants","fly",
#                            "fly","no_fly","fly","fly","fly",
#                            "fly", "fly", "fly", "no_fly", "fly")
# 
# ele_data_analysis$animal <- c("animal", "animal", "animal","animal", "animal",
#                               "animal","animal","animal", "animal", "animal", 
#                               "animal", "animal", "animal", "animal", "animal",
#                               "plants", "plants", "animal", "animal", "animal",
#                               "animal", "animal", "animal", "animal", "animal",
#                               "animal","animal", "animal", "animal","animal",
#                               "animal", "plants", "plants", "plants", "animal",
#                               "animal", "animal", "animal", "animal", "animal",
#                               "animal", "animal", "animal", "animal", "animal",
#                               "animal","animal","animal","animal","animal",
#                               "animal", "animal", "animal", "animal", "animal",
#                               "animal","animal", "animal","animal", "animal",
#                               "plants","animal", "plants","plants","animal",
#                               "animal","animal","animal","animal","animal",
#                               "animal", "animal", "animal", "animal", "animal")
# 
# ele_data_analysis$log10rich <- log10(ele_data_analysis$richness)
# 
# # Edwards 2010 is BIRDS, not PLANTS
# 
# # models
# library(glmmTMB)
# library(effects)
# 
# 
# model <- glmmTMB(SS ~ 
#                    #animal +
#                    #vert +
#                    #taxa +
#                    #group +
#                    #fly +
#                    log10rich +
#                    #richness +
#                    #temperature +
#                    evenness + 
#                    #sphere.fragment + 
#                    #sphere.matrix + 
#                    #time.since.fragmentation + 
#                    #Matrix.category + 
#                    #biome + 
#                    #(1|temperature/taxa)+
#                    #(1|taxa) +
#                    (1|dataset_id), 
#                  data = ele_data_analysis, family = "binomial")
# summary(model)
# AIC(model)
# 
# #plot(allEffects(model), type = "response")
# plot_model <- ggpredict(model, 
#                         c("evenness [all]", "log10rich [1, 1.5, 2]" ),
#                         type = "re"
# ) 
# 
# plot(plot_model) + 
#   labs(
#     x = "Patch size evenness", 
#     y = "Probability of SS > SL", 
#     title = ""  ) + 
#   labs(colour = "Log10(species richness)")
# 
# table(ele_data_analysis$comparison, ele_data_analysis$taxa)

# boxplot <- ele_data_analysis %>%
#   ggplot( aes(x=ele_data_analysis$taxa, y= ele_data_analysis$log10rich, fill= ele_data_analysis$taxa)) +
#   #geom_violin() +
#   scale_fill_manual(values=c("gainsboro", "cyan2", "blue", "red", "yellow")) +
#   #geom_jitter(color="black", size=1, alpha=0.9) +
#   theme_bw()+
#   theme(
#     legend.position="none",
#     plot.title = element_text(size=11)
#   ) +
#   #ggtitle("True diversity of complexity themes in a publication") +
#   geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7) +
#   geom_boxplot(width=0.1, outlier.shape = NA, alpha = 0.8) +
#   #ylim(0,1)+
#   ylab("True diversity")+
#   xlab("") 
#   #geom_signif(comparisons = list(c("Control", "Ecological complexity")), 
#   #            map_signif_level=TRUE)
# boxplot

##
## data ready for analysis. 
##



## DON'T REMEMBER WHAT THIS IS; OCT 16 2021
# #nrow(list_patches[[i]]) ## this gives you how to 
# 
# for (i in 1:length(list_sampled_data)) {
#   for(j in 1:1186){
#   #list_sampled_data[[i]] <-subset(list_sampled_data[[i]], select = c("dataset_label", "frag_id")) 
#   list_sampled_data[[i]] <- add_count(list_sampled_data[[i]], species, study_and_patch)
#   list_sampled_data[[i]] <- unique(list_sampled_data[[i]])
#   list_sampled_data[[i]] <- reshape(list_sampled_data[[i]],                                 # Applying reshape function
#                            idvar = "study_and_patch",
#                            timevar = "species",
#                            direction = "wide")
# 
#   }
# }
# prova <- list_sampled_data[[1]]
# 
# 
# prova <- list_sampled_data[[1]][2]
# prova <- (unlist(prova))
# 
# prova <- add_count(prova, species, study_and_patch)
# prova <- unique(prova)
# prova <- reshape(prova,                                 # Applying reshape function
#                                   idvar = "study_and_patch",
#                                   timevar = "species",
#                                   direction = "wide")


































########################################################################################################
################################################################################################################################################################################################################
################################################################################################################################################################################################################
################################################################################################################################################################################################################
# ################################################################################################################################################################################################################
# ################################################################################################################################################################################################################
# ########################################################################################################
# list_sim_twenty <- list()
# for (i in 1:100) {list_sim_twenty[[i]] <- list_simulations_twenty}
# 
# # merge the simulated sets of patches with the simulated assemblages in each patch
# for (i in 1: length(list_sim_twenty)) {
#   for (j in 1: length(list_sim_twenty[[i]])) {
#     for (k in 1: length(list_sim_twenty[[i]][[j]])){
#       
#       if (length(list_sim_twenty[[i]][[j]]) > 1) {
#         
#         list_sim_twenty[[i]][[j]][[k]] <- merge(list_sim_twenty[[i]][[j]][[k]], 
#                                                 list_sampled_data[[i]][j][[1]],
#                                                 by = "study_and_patch")
#         
#       } else  {
#         
#         list_sim_twenty[[i]][[j]] <- NA
#       }
#       
#     }
#   }
# }
# 
# 
# 
# # calculate the average patch areas in each set of patches
# list_sim_twenty_area <- list_sim_twenty
# 
# for (i in 1: length(list_sim_twenty_area)) {
#   for (j in 1: length(list_sim_twenty_area[[i]])) {
#     for (k in 1: length(list_sim_twenty_area[[i]][[j]])){
#       
#       if (!is.na(list_sim_twenty_area[[i]][[j]])) {
#         
#         list_sim_twenty_area[[i]][[j]][[k]] <- mean(list_sim_twenty_area[[i]][[j]][[k]][,ncol(list_sim_twenty_area[[i]][[j]][[k]])])
#         
#       } 
#       
#     }
#   }
# }
# # warnings appears because every list has more than one element, so the !is.na returns multiple TRUE statements when a list is not NA
# 
# # transform lists into tabels
# for (i in 1:length(list_sim_twenty_area)){
#   for ( j in 1: length(list_sim_twenty_area[[i]])) {
#     
#     if (!is.na(list_sim_twenty_area[[i]][[j]])) {
#       list_sim_twenty_area[[i]][[j]] <- do.call(rbind.data.frame, list_sim_twenty_area[[i]][[j]])
#       colnames(list_sim_twenty_area[[i]][[j]]) <- "mean patch area"
#     }
#   }
# }
# 
# 
# 
# # calculate the species richness in each patch
# list_sim_twenty_richness <- list_sim_twenty
# 
# for (i in 1: length(list_sim_twenty_richness)) {
#   for (j in 1: length(list_sim_twenty_richness[[i]])) {
#     for (k in 1: length(list_sim_twenty_richness[[i]][[j]])){
#       
#       if (!is.na(list_sim_twenty_richness[[i]][[j]])) {
#         # take the column of each table from column two to column n-1, take their sum, and count how many times their colsum is larger than 0
#         list_sim_twenty_richness[[i]][[j]][[k]] <- sum(colSums(list_sim_twenty_richness[[i]][[j]][[k]][,2: (ncol(list_sim_twenty_richness[[i]][[j]][[k]])-1)]) > 0)
#       } 
#       
#     }
#   }
# }
# 
# 
# for (i in 1:length(list_sim_twenty_richness)){
#   for ( j in 1: length(list_sim_twenty_richness[[i]])) {
#     
#     if (!is.na(list_sim_twenty_richness[[i]][[j]])) {
#       list_sim_twenty_richness[[i]][[j]] <- do.call(rbind.data.frame, list_sim_twenty_richness[[i]][[j]])
#       colnames(list_sim_twenty_richness[[i]][[j]]) <- "richness"
#     }
#   }
# }
# 
# 
# 
# 
# 
# # calculate evenness in each patch
# list_sim_twenty_evenness <- list_sim_twenty
# 
# for (i in 1: length(list_sim_twenty_evenness)) {
#   for (j in 1: length(list_sim_twenty_evenness[[i]])) {
#     for (k in 1: length(list_sim_twenty_evenness[[i]][[j]])){
#       
#       if (!is.na(list_sim_twenty_evenness[[i]][[j]])) {
#         # use the HurlbertPIE function on every table, excluding patch name (column 1) and patch size (column nrow)
#         list_sim_twenty_evenness[[i]][[j]][[k]] <- paleotree::HurlbertPIE(colSums(list_sim_twenty_evenness[[i]][[j]][[k]][,2: (ncol(list_sim_twenty_evenness[[i]][[j]][[k]])-1)]))
#       } 
#       
#     }
#   }
# }
# 
# 
# 
# for (i in 1:length(list_sim_twenty_evenness)){
#   for ( j in 1: length(list_sim_twenty_evenness[[i]])) {
#     
#     if (!is.na(list_sim_twenty_evenness[[i]][[j]])) {
#       list_sim_twenty_evenness[[i]][[j]] <- do.call(rbind.data.frame, list_sim_twenty_evenness[[i]][[j]])
#       colnames(list_sim_twenty_evenness[[i]][[j]]) <- "evenness"
#     }
#   }
# }
# 
# 
# 
# 
# list_sim_twenty_combined <- list_sim_twenty_area
# 
# for (i in 1:length(list_sim_twenty_combined)){
#   for ( j in 1: length(list_sim_twenty_combined[[i]])) {
#     
#     if (!is.na(list_sim_twenty_combined[[i]][[j]])) {
#       list_sim_twenty_combined[[i]][[j]] <- cbind( list_sim_twenty_area[[i]][[j]],  
#                                                    list_sim_twenty_richness[[i]][[j]],  
#                                                    list_sim_twenty_evenness[[i]][[j]])
#        list_sim_twenty_combined[[i]][[j]]$simulation_number <- rep(i, # give a number from 1 to 100 
#                                                                    nrow(list_sim_twenty_combined[[i]][[j]]))
#        list_sim_twenty_combined[[i]][[j]]$study <- rep(names(list_sampled_data[[i]][j]), 
#                                                        nrow(list_sim_twenty_combined[[i]][[j]]))
#     }
#   }
# }
# 
# 
# ## combine all studies in a table for each of the 100 simulations
# for (i in 1:length(list_sim_twenty_combined)){
#       list_sim_twenty_combined[[i]] <- do.call(rbind.data.frame, list_sim_twenty_combined[[i]])
# 
#       }
# 
# list_sim_twenty_combined <- do.call(rbind, list_sim_twenty_combined)
# list_sim_twenty_combined <- na.omit(list_sim_twenty_combined)
# 
# final_twenty <- list_sim_twenty_combined
# rm(list_sim_twenty)
# 
# ########################################################################################################
# ########################################################################################################
# ########################################################################################################
# 
# 
# 
# list_sim_forty <- list()
# for (i in 1:100) {list_sim_forty[[i]] <- list_simulations_forty}
# 
# # merge the simulated sets of patches with the simulated assemblages in each patch
# for (i in 1: length(list_sim_forty)) {
#   for (j in 1: length(list_sim_forty[[i]])) {
#     for (k in 1: length(list_sim_forty[[i]][[j]])){
#       
#       if (length(list_sim_forty[[i]][[j]]) > 1) {
#         
#         list_sim_forty[[i]][[j]][[k]] <- merge(list_sim_forty[[i]][[j]][[k]], 
#                                                 list_sampled_data[[i]][j][[1]],
#                                                 by = "study_and_patch")
#         
#       } else  {
#         
#         list_sim_forty[[i]][[j]] <- NA
#       }
#       
#     }
#   }
# }
# 
# 
# 
# # calculate the average patch areas in each set of patches
# list_sim_forty_area <- list_sim_forty
# 
# for (i in 1: length(list_sim_forty_area)) {
#   for (j in 1: length(list_sim_forty_area[[i]])) {
#     for (k in 1: length(list_sim_forty_area[[i]][[j]])){
#       
#       if (!is.na(list_sim_forty_area[[i]][[j]])) {
#         
#         list_sim_forty_area[[i]][[j]][[k]] <- mean(list_sim_forty_area[[i]][[j]][[k]][,ncol(list_sim_forty_area[[i]][[j]][[k]])])
#         
#       } 
#       
#     }
#   }
# }
# # warnings appears because every list has more than one element, so the !is.na returns multiple TRUE statements when a list is not NA
# 
# # transform lists into tabels
# for (i in 1:length(list_sim_forty_area)){
#   for ( j in 1: length(list_sim_forty_area[[i]])) {
#     
#     if (!is.na(list_sim_forty_area[[i]][[j]])) {
#       list_sim_forty_area[[i]][[j]] <- do.call(rbind.data.frame, list_sim_forty_area[[i]][[j]])
#       colnames(list_sim_forty_area[[i]][[j]]) <- "mean patch area"
#     }
#   }
# }
# 
# 
# 
# # calculate the species richness in each patch
# list_sim_forty_richness <- list_sim_forty
# 
# for (i in 1: length(list_sim_forty_richness)) {
#   for (j in 1: length(list_sim_forty_richness[[i]])) {
#     for (k in 1: length(list_sim_forty_richness[[i]][[j]])){
#       
#       if (!is.na(list_sim_forty_richness[[i]][[j]])) {
#         # take the column of each table from column two to column n-1, take their sum, and count how many times their colsum is larger than 0
#         list_sim_forty_richness[[i]][[j]][[k]] <- sum(colSums(list_sim_forty_richness[[i]][[j]][[k]][,2: (ncol(list_sim_forty_richness[[i]][[j]][[k]])-1)]) > 0)
#       } 
#       
#     }
#   }
# }
# 
# 
# for (i in 1:length(list_sim_forty_richness)){
#   for ( j in 1: length(list_sim_forty_richness[[i]])) {
#     
#     if (!is.na(list_sim_forty_richness[[i]][[j]])) {
#       list_sim_forty_richness[[i]][[j]] <- do.call(rbind.data.frame, list_sim_forty_richness[[i]][[j]])
#       colnames(list_sim_forty_richness[[i]][[j]]) <- "richness"
#     }
#   }
# }
# 
# 
# 
# 
# 
# # calculate evenness in each patch
# list_sim_forty_evenness <- list_sim_forty
# 
# for (i in 1: length(list_sim_forty_evenness)) {
#   for (j in 1: length(list_sim_forty_evenness[[i]])) {
#     for (k in 1: length(list_sim_forty_evenness[[i]][[j]])){
#       
#       if (!is.na(list_sim_forty_evenness[[i]][[j]])) {
#         # use the HurlbertPIE function on every table, excluding patch name (column 1) and patch size (column nrow)
#         list_sim_forty_evenness[[i]][[j]][[k]] <- paleotree::HurlbertPIE(colSums(list_sim_forty_evenness[[i]][[j]][[k]][,2: (ncol(list_sim_forty_evenness[[i]][[j]][[k]])-1)]))
#       } 
#       
#     }
#   }
# }
# 
# 
# 
# for (i in 1:length(list_sim_forty_evenness)){
#   for ( j in 1: length(list_sim_forty_evenness[[i]])) {
#     
#     if (!is.na(list_sim_forty_evenness[[i]][[j]])) {
#       list_sim_forty_evenness[[i]][[j]] <- do.call(rbind.data.frame, list_sim_forty_evenness[[i]][[j]])
#       colnames(list_sim_forty_evenness[[i]][[j]]) <- "evenness"
#     }
#   }
# }
# 
# 
# 
# 
# list_sim_forty_combined <- list_sim_forty_area
# 
# for (i in 1:length(list_sim_forty_combined)){
#   for ( j in 1: length(list_sim_forty_combined[[i]])) {
#     
#     if (!is.na(list_sim_forty_combined[[i]][[j]])) {
#       list_sim_forty_combined[[i]][[j]] <- cbind( list_sim_forty_area[[i]][[j]],  
#                                                    list_sim_forty_richness[[i]][[j]],  
#                                                    list_sim_forty_evenness[[i]][[j]])
#       list_sim_forty_combined[[i]][[j]]$simulation_number <- rep(i, # give a number from 1 to 100 
#                                                                   nrow(list_sim_forty_combined[[i]][[j]]))
#       list_sim_forty_combined[[i]][[j]]$study <- rep(names(list_sampled_data[[i]][j]), 
#                                                       nrow(list_sim_forty_combined[[i]][[j]]))
#     }
#   }
# }
# 
# 
# ## combine all studies in a table for each of the 100 simulations
# for (i in 1:length(list_sim_forty_combined)){
#   list_sim_forty_combined[[i]] <- do.call(rbind.data.frame, list_sim_forty_combined[[i]])
#   
# }
# 
# list_sim_forty_combined <- do.call(rbind, list_sim_forty_combined)
# list_sim_forty_combined <- na.omit(list_sim_forty_combined)
# 
# final_forty <- list_sim_forty_combined
# rm(list_sim_forty)
# 
# ########################################################################################################
# ########################################################################################################
# ########################################################################################################
# 
# ########################################################################################################
# ########################################################################################################
# ########################################################################################################
# 
# 
# 
# list_sim_sixty <- list()
# for (i in 1:100) {list_sim_sixty[[i]] <- list_simulations_sixty}
# 
# # merge the simulated sets of patches with the simulated assemblages in each patch
# for (i in 1: length(list_sim_sixty)) {
#   for (j in 1: length(list_sim_sixty[[i]])) {
#     for (k in 1: length(list_sim_sixty[[i]][[j]])){
#       
#       if (length(list_sim_sixty[[i]][[j]]) > 1) {
#         
#         list_sim_sixty[[i]][[j]][[k]] <- merge(list_sim_sixty[[i]][[j]][[k]], 
#                                                list_sampled_data[[i]][j][[1]],
#                                                by = "study_and_patch")
#         
#       } else  {
#         
#         list_sim_sixty[[i]][[j]] <- NA
#       }
#       
#     }
#   }
# }
# 
# 
# 
# # calculate the average patch areas in each set of patches
# list_sim_sixty_area <- list_sim_sixty
# 
# for (i in 1: length(list_sim_sixty_area)) {
#   for (j in 1: length(list_sim_sixty_area[[i]])) {
#     for (k in 1: length(list_sim_sixty_area[[i]][[j]])){
#       
#       if (!is.na(list_sim_sixty_area[[i]][[j]])) {
#         
#         list_sim_sixty_area[[i]][[j]][[k]] <- mean(list_sim_sixty_area[[i]][[j]][[k]][,ncol(list_sim_sixty_area[[i]][[j]][[k]])])
#         
#       } 
#       
#     }
#   }
# }
# # warnings appears because every list has more than one element, so the !is.na returns multiple TRUE statements when a list is not NA
# 
# # transform lists into tabels
# for (i in 1:length(list_sim_sixty_area)){
#   for ( j in 1: length(list_sim_sixty_area[[i]])) {
#     
#     if (!is.na(list_sim_sixty_area[[i]][[j]])) {
#       list_sim_sixty_area[[i]][[j]] <- do.call(rbind.data.frame, list_sim_sixty_area[[i]][[j]])
#       colnames(list_sim_sixty_area[[i]][[j]]) <- "mean patch area"
#     }
#   }
# }
# 
# 
# 
# # calculate the species richness in each patch
# list_sim_sixty_richness <- list_sim_sixty
# 
# for (i in 1: length(list_sim_sixty_richness)) {
#   for (j in 1: length(list_sim_sixty_richness[[i]])) {
#     for (k in 1: length(list_sim_sixty_richness[[i]][[j]])){
#       
#       if (!is.na(list_sim_sixty_richness[[i]][[j]])) {
#         # take the column of each table from column two to column n-1, take their sum, and count how many times their colsum is larger than 0
#         list_sim_sixty_richness[[i]][[j]][[k]] <- sum(colSums(list_sim_sixty_richness[[i]][[j]][[k]][,2: (ncol(list_sim_sixty_richness[[i]][[j]][[k]])-1)]) > 0)
#       } 
#       
#     }
#   }
# }
# 
# 
# for (i in 1:length(list_sim_sixty_richness)){
#   for ( j in 1: length(list_sim_sixty_richness[[i]])) {
#     
#     if (!is.na(list_sim_sixty_richness[[i]][[j]])) {
#       list_sim_sixty_richness[[i]][[j]] <- do.call(rbind.data.frame, list_sim_sixty_richness[[i]][[j]])
#       colnames(list_sim_sixty_richness[[i]][[j]]) <- "richness"
#     }
#   }
# }
# 
# 
# 
# 
# 
# # calculate evenness in each patch
# list_sim_sixty_evenness <- list_sim_sixty
# 
# for (i in 1: length(list_sim_sixty_evenness)) {
#   for (j in 1: length(list_sim_sixty_evenness[[i]])) {
#     for (k in 1: length(list_sim_sixty_evenness[[i]][[j]])){
#       
#       if (!is.na(list_sim_sixty_evenness[[i]][[j]])) {
#         # use the HurlbertPIE function on every table, excluding patch name (column 1) and patch size (column nrow)
#         list_sim_sixty_evenness[[i]][[j]][[k]] <- paleotree::HurlbertPIE(colSums(list_sim_sixty_evenness[[i]][[j]][[k]][,2: (ncol(list_sim_sixty_evenness[[i]][[j]][[k]])-1)]))
#       } 
#       
#     }
#   }
# }
# 
# 
# 
# for (i in 1:length(list_sim_sixty_evenness)){
#   for ( j in 1: length(list_sim_sixty_evenness[[i]])) {
#     
#     if (!is.na(list_sim_sixty_evenness[[i]][[j]])) {
#       list_sim_sixty_evenness[[i]][[j]] <- do.call(rbind.data.frame, list_sim_sixty_evenness[[i]][[j]])
#       colnames(list_sim_sixty_evenness[[i]][[j]]) <- "evenness"
#     }
#   }
# }
# 
# 
# 
# 
# list_sim_sixty_combined <- list_sim_sixty_area
# 
# for (i in 1:length(list_sim_sixty_combined)){
#   for ( j in 1: length(list_sim_sixty_combined[[i]])) {
#     
#     if (!is.na(list_sim_sixty_combined[[i]][[j]])) {
#       list_sim_sixty_combined[[i]][[j]] <- cbind( list_sim_sixty_area[[i]][[j]],  
#                                                   list_sim_sixty_richness[[i]][[j]],  
#                                                   list_sim_sixty_evenness[[i]][[j]])
#       list_sim_sixty_combined[[i]][[j]]$simulation_number <- rep(i, # give a number from 1 to 100 
#                                                                  nrow(list_sim_sixty_combined[[i]][[j]]))
#       list_sim_sixty_combined[[i]][[j]]$study <- rep(names(list_sampled_data[[i]][j]), 
#                                                      nrow(list_sim_sixty_combined[[i]][[j]]))
#     }
#   }
# }
# 
# 
# ## combine all studies in a table for each of the 100 simulations
# for (i in 1:length(list_sim_sixty_combined)){
#   list_sim_sixty_combined[[i]] <- do.call(rbind.data.frame, list_sim_sixty_combined[[i]])
#   
# }
# 
# list_sim_sixty_combined <- do.call(rbind, list_sim_sixty_combined)
# list_sim_sixty_combined <- na.omit(list_sim_sixty_combined)
# 
# final_sixty <- list_sim_sixty_combined
# rm(list_sim_sixty)
# 
# ########################################################################################################
# ########################################################################################################
# ########################################################################################################
# 
# 
# ########################################################################################################
# ########################################################################################################
# ########################################################################################################
# 
# 
# 
# list_sim_eighty <- list()
# for (i in 1:100) {list_sim_eighty[[i]] <- list_simulations_eighty}
# 
# # merge the simulated sets of patches with the simulated assemblages in each patch
# for (i in 1: length(list_sim_eighty)) {
#   for (j in 1: length(list_sim_eighty[[i]])) {
#     for (k in 1: length(list_sim_eighty[[i]][[j]])){
#       
#       if (length(list_sim_eighty[[i]][[j]]) > 1) {
#         
#         list_sim_eighty[[i]][[j]][[k]] <- merge(list_sim_eighty[[i]][[j]][[k]], 
#                                                list_sampled_data[[i]][j][[1]],
#                                                by = "study_and_patch")
#         
#       } else  {
#         
#         list_sim_eighty[[i]][[j]] <- NA
#       }
#       
#     }
#   }
# }
# 
# 
# 
# # calculate the average patch areas in each set of patches
# list_sim_eighty_area <- list_sim_eighty
# 
# for (i in 1: length(list_sim_eighty_area)) {
#   for (j in 1: length(list_sim_eighty_area[[i]])) {
#     for (k in 1: length(list_sim_eighty_area[[i]][[j]])){
#       
#       if (!is.na(list_sim_eighty_area[[i]][[j]])) {
#         
#         list_sim_eighty_area[[i]][[j]][[k]] <- mean(list_sim_eighty_area[[i]][[j]][[k]][,ncol(list_sim_eighty_area[[i]][[j]][[k]])])
#         
#       } 
#       
#     }
#   }
# }
# # warnings appears because every list has more than one element, so the !is.na returns multiple TRUE statements when a list is not NA
# 
# # transform lists into tabels
# for (i in 1:length(list_sim_eighty_area)){
#   for ( j in 1: length(list_sim_eighty_area[[i]])) {
#     
#     if (!is.na(list_sim_eighty_area[[i]][[j]])) {
#       list_sim_eighty_area[[i]][[j]] <- do.call(rbind.data.frame, list_sim_eighty_area[[i]][[j]])
#       colnames(list_sim_eighty_area[[i]][[j]]) <- "mean patch area"
#     }
#   }
# }
# 
# 
# 
# # calculate the species richness in each patch
# list_sim_eighty_richness <- list_sim_eighty
# 
# for (i in 1: length(list_sim_eighty_richness)) {
#   for (j in 1: length(list_sim_eighty_richness[[i]])) {
#     for (k in 1: length(list_sim_eighty_richness[[i]][[j]])){
#       
#       if (!is.na(list_sim_eighty_richness[[i]][[j]])) {
#         # take the column of each table from column two to column n-1, take their sum, and count how many times their colsum is larger than 0
#         list_sim_eighty_richness[[i]][[j]][[k]] <- sum(colSums(list_sim_eighty_richness[[i]][[j]][[k]][,2: (ncol(list_sim_eighty_richness[[i]][[j]][[k]])-1)]) > 0)
#       } 
#       
#     }
#   }
# }
# 
# 
# for (i in 1:length(list_sim_eighty_richness)){
#   for ( j in 1: length(list_sim_eighty_richness[[i]])) {
#     
#     if (!is.na(list_sim_eighty_richness[[i]][[j]])) {
#       list_sim_eighty_richness[[i]][[j]] <- do.call(rbind.data.frame, list_sim_eighty_richness[[i]][[j]])
#       colnames(list_sim_eighty_richness[[i]][[j]]) <- "richness"
#     }
#   }
# }
# 
# 
# 
# 
# 
# # calculate evenness in each patch
# list_sim_eighty_evenness <- list_sim_eighty
# 
# for (i in 1: length(list_sim_eighty_evenness)) {
#   for (j in 1: length(list_sim_eighty_evenness[[i]])) {
#     for (k in 1: length(list_sim_eighty_evenness[[i]][[j]])){
#       
#       if (!is.na(list_sim_eighty_evenness[[i]][[j]])) {
#         # use the HurlbertPIE function on every table, excluding patch name (column 1) and patch size (column nrow)
#         list_sim_eighty_evenness[[i]][[j]][[k]] <- paleotree::HurlbertPIE(colSums(list_sim_eighty_evenness[[i]][[j]][[k]][,2: (ncol(list_sim_eighty_evenness[[i]][[j]][[k]])-1)]))
#       } 
#       
#     }
#   }
# }
# 
# 
# 
# for (i in 1:length(list_sim_eighty_evenness)){
#   for ( j in 1: length(list_sim_eighty_evenness[[i]])) {
#     
#     if (!is.na(list_sim_eighty_evenness[[i]][[j]])) {
#       list_sim_eighty_evenness[[i]][[j]] <- do.call(rbind.data.frame, list_sim_eighty_evenness[[i]][[j]])
#       colnames(list_sim_eighty_evenness[[i]][[j]]) <- "evenness"
#     }
#   }
# }
# 
# 
# 
# 
# list_sim_eighty_combined <- list_sim_eighty_area
# 
# for (i in 1:length(list_sim_eighty_combined)){
#   for ( j in 1: length(list_sim_eighty_combined[[i]])) {
#     
#     if (!is.na(list_sim_eighty_combined[[i]][[j]])) {
#       list_sim_eighty_combined[[i]][[j]] <- cbind( list_sim_eighty_area[[i]][[j]],  
#                                                   list_sim_eighty_richness[[i]][[j]],  
#                                                   list_sim_eighty_evenness[[i]][[j]])
#       list_sim_eighty_combined[[i]][[j]]$simulation_number <- rep(i, # give a number from 1 to 100 
#                                                                  nrow(list_sim_eighty_combined[[i]][[j]]))
#       list_sim_eighty_combined[[i]][[j]]$study <- rep(names(list_sampled_data[[i]][j]), 
#                                                      nrow(list_sim_eighty_combined[[i]][[j]]))
#     }
#   }
# }
# 
# 
# ## combine all studies in a table for each of the 100 simulations
# for (i in 1:length(list_sim_eighty_combined)){
#   list_sim_eighty_combined[[i]] <- do.call(rbind.data.frame, list_sim_eighty_combined[[i]])
#   
# }
# 
# list_sim_eighty_combined <- do.call(rbind, list_sim_eighty_combined)
# list_sim_eighty_combined <- na.omit(list_sim_eighty_combined)
# 
# final_eighty <- list_sim_eighty_combined
# rm(list_sim_eighty)

########################################################################################################
########################################################################################################
########################################################################################################