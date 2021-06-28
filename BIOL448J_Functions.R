# Ilan Rubin
# October 17, 2019
# BIOL 448J Function

# This file includes all the functions I have written for BIOL 448J through week 5.


###################################
# Biodiversity Measures Functions #
###################################

### All of these functions take a vector 'count'.
### This vectors is simply a list of the number of each species/OTU/etc.
###   found in a single sample.
### e.g., if, in one sample, we found 3 of species A, 12, of species B,
###   0 of species C, and 1 of species D, the vector count = c(3,12,0,1).

### Observed Richness ###
# A simple count of the number of species found.
# Usually denoted by R or S
richness = function(abundance_table){
  if (is.vector(abundance_table)){
    abundance_table = matrix(abundance_table,1,length(abundance_table))
  }
  count = colSums(abundance_table,na.rm=TRUE)
  sum(count>0)
}

### Estimated Richness ###
# We will likely go over this next week in lab.
# Basically, we need to fit some funcion to our collector's curve and find it's asymptote.

### Simpson's Diversity Index ###
# Usually denoted by lambda (or l)
# Equals the probability any two randomly chosen individuals in the population are the same.
# Therefore, actually a similarity index. To get diversity, D = 1/l
# Measure of both richness and evenness.
simpsons = function(abundance_table){
  if (is.vector(abundance_table)){
    abundance_table = matrix(abundance_table,1,length(abundance_table))
  }
  count = colSums(abundance_table,na.rm=TRUE)
  p = count/sum(count)
  sum(p^2)
}

### Shannon's Diversity Index (entropy) ###
# Usualy denoted by H'
# Does not have a specific mathematical explanation, but increases as diversity increases.
# Assumes infinite, well-mixed population.
# Also a measure of both richness and evenness.
shannons = function(abundance_table){
  if (is.vector(abundance_table)){
    abundance_table = matrix(abundance_table,1,length(abundance_table))
  }
  count = colSums(abundance_table,na.rm=TRUE)
  count = count[count>0]
  p = count/sum(count)
  -sum(p*log(p^2))
}


### Chao1 ###
# Chao1 is a non-parametric estimator of richness. That means it makes no assumptions (infinite population,
# well mixed, specific model to fit and find an asymptote). It is observed richness, corrected
# with the number of singletons (species only seen once) and doubletons (species seen twice).
# There is actually some fancy statistics behind using singletons and doubletons as correctors.
chao1 = function(abundance_table){
  if (is.vector(abundance_table)){
    abundance_table = matrix(abundance_table,1,length(abundance_table))
  }
  count = colSums(abundance_table,na.rm=TRUE)
  f1 = sum(count == 1)
  f2 = sum(count == 2)
  S = sum(count > 0)
  return(S + f1*(f1-1)/(2*(f2+1)))
}  

### Chao2 ###
# Similar to chao1 but uses incidence (pressence/absence) data instead of abundance.
# Running this on a single sample will return observed richness.
# instead of using doubleton and singleton species based on abundance, it uses species
# that were only found in one or two samples.
# This function will convert abundance data to incidence. If you have incidence data already,
# it can either be in 1s (pressence) and 0s (absence) or TRUE and FALSE.
chao2 = function(abundance_table){
  abundance_table[abundance_table > 0] = 1
  if (is.vector(abundance_table)){
    abundance_table = matrix(abundance_table,1,length(abundance_table))
  }
  count = colSums(abundance_table,na.rm=TRUE)
  q1 = sum(count == 1)
  q2 = sum(count == 2)
  S = sum(count > 0)
  m = nrow(abundance_table)
  return(S + (m-1)/m * q1*(q1-1)/(2*(q2+1)))
}


### Run all diversity and richness functions ###
# Specify which columns are the species counts if metadata is also included in the
# table. Optionally, specify a column to filter by (e.g., to run the diversity indices
# for each lake in a column name Lake use filter='Lake').
run_all_diversity = function(abundance_table,species_col=TRUE,filter_by=""){
  if (filter_by != ""){
    filter_unique = unique(abundance_table[[filter_by]])
    num_unique = length(filter_unique)
  }else{
    filter_unique = TRUE
    num_unique = 1
    which_row = TRUE
  }
  
  richness_col = rep(0,num_unique)
  shannons_col = rep(0,num_unique)
  simpsons_col = rep(0,num_unique)
  inv_simpsons_col = rep(0,num_unique)
  chao1_col = rep(0,num_unique)
  chao2_col = rep(0,num_unique)
  N_col = rep(0,num_unique)
  
  for (u in 1:num_unique){
    if (filter_by != ""){
      which_row = abundance_table[[filter_by]]==filter_unique[u]
    }
    richness_col[u] = richness(abundance_table[which_row,species_col])
    shannons_col[u] = shannons(abundance_table[which_row,species_col])
    simpsons_col[u] = simpsons(abundance_table[which_row,species_col])
    inv_simpsons_col[u] = 1/simpsons(abundance_table[which_row,species_col])
    chao1_col[u] = chao1(abundance_table[which_row,species_col])
    chao2_col[u] = chao2(abundance_table[which_row,species_col])
    N_col[u] = sum(which_row)
  }
  if (filter_by != ""){
    output = data_frame(N=N_col,
                        filter=filter_unique,
                        richness=richness_col,
                        shannons=shannons_col,
                        simpsons=simpsons_col,
                        inv_simpsons=inv_simpsons_col,
                        chao1=chao1_col,
                        chao2=chao2_col)
    names(output)[2] = filter_by
  }else{
    output = data_frame(N=N_col,
                        richness=richness_col,
                        shannons=shannons_col,
                        simpsons=simpsons_col,
                        inv_simpsons=inv_simpsons_col,
                        chao1=chao1_col,
                        chao2=chao2_col)
  }
  return(output)
}

### Rarefaction function ###
# Takes an abundance table and calculates a rarefaction curve for a specified function 'fun'.
# The default calculates a richness rarefaction, but any function that takes an abundance table works.
# Change the input 'species_col' to indicate which columns are the species counts if the table also
#     includes metadata.
# Change the input 'rand_combinations' to a non-zero integer to use only that number of random samples
#     instead of all possible. e.g., rand_combinations=100 would take the average of richness (or other fun)
#     for 100 different combinations of random samples (these can add up very quickly!). If the number
#     of samples in your data set is very smaple (~10 or so), you can set rand_combinations=Inf to use all
#     possible combinations.
# If you want to run the rarefaction multiple times filtered for each unique value in a specific column,
#     set the input 'filter_col' equal to the name of that column.
# Turn 'verbose' to FALSE if you don't want to see the output as it works.
rarefaction = function(abundance_table,
                       fun=richness,
                       species_col=T,
                       rand_combinations=100,
                       filter_col="",
                       verbose=TRUE){
  if (filter_col != ""){
    filter_unique = unique(abundance_table[[filter_col]])
    filter_num = which(names(abundance_table)==filter_col)
    num_unique = length(filter_unique)
  }else{
    num_unique = 1
  }
  
  n_sample = nrow(abundance_table[,species_col])    # the number of samples in the dataset
  output = rep(0,n_sample*num_unique)
  output_95 = rep(0,n_sample*num_unique)
  
  # loop over each unique element in the filter column if specified
  for (u in 1:num_unique){
    # filter the abundance table for the first element
    if (filter_col != ""){
      filter_abundance_table = filter_at(abundance_table,filter_num,all_vars(.==filter_unique[u]))
    }else{
      filter_abundance_table = abundance_table
    }
    filter_abundance_table = filter_abundance_table[,species_col]   # keep only the species count data
    
    # loop from 1 to the number of samples
    for (n in 1:n_sample){
      n_combinations = choose(n_sample,n) # the total number of combinations of n samples
      
      # check if the number of combinations is larger than specified
      if (rand_combinations<n_combinations){
        n_combinations = rand_combinations
        sample_order = matrix(0,n,rand_combinations)
        i = 1
        while (i <= rand_combinations){ # loop over the number of combinations
          sample_order[,i] = sample(n_sample,n) # random order of samples
          for (j in 1:i){ # check to make sure the sample has not already been chosen
            if (j != i){
              if (all(sample_order[,i]-sample_order[,j] == 0)){
                i = i-1
                break
              }
            }
          }
          i = 1+i
        }
      }else{
        sample_order = combn(n_sample,n) # find all possible combinations of n samples
      }
      
      temp_out = rep(0,n_combinations)
      for (i in 1:n_combinations){
        temp_sample = filter_abundance_table[sample_order[,i],] # use only chosen samples
        temp_out[i] = fun(temp_sample) # calculate function
      }
      output[n + (u-1)*n_sample] = mean(temp_out)
      output_95[n + (u-1)*n_sample] = 1.96*sd(temp_out)/sqrt(i) # calculate 95 confidence interval of the mean
      if (verbose) print(paste(n+(u-1)*n_sample,'/',num_unique*n_sample,', Number of Samples:',n_combinations,', Value:',output[n],sep=""))
    }
  }
  if (filter_col != ""){
    return(tibble(n_samples=rep(1:n_sample,times=num_unique),
                  value=output,
                  conf_95=output_95,
                  filter=rep(filter_unique,each=n_sample)))
  }else{
    return(tibble(n_samples=1:n_sample,
                  value=output,
                  conf_95=output_95))
  }
}


#################
# MDS Functions #
#################
# Function to run the ellipse calculation for each of column specified by group
# Uasge: calculate_ellipse(seagrass,"region")
# Will output a data frame with three columns,
#    one called "filter" for the column you are filtering by
#    and one each for MDS1 and MDS2
calculate_ellipse = function(df,filter_by){
  # Function to calculate the 95% ci ellipses for each sample
  veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) 
  {
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))
  }
  
  df_ellipse <- data.frame()
  for(g in unique(df[[filter_by]])){
    df_ellipse <- rbind(df_ellipse,cbind(as.data.frame(with(df [df[filter_by]==g,],
                                                            veganCovEllipse(cov.wt(cbind(MDS1,MDS2),
                                                                                   wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2))))),filter=g))
  }
  return(df_ellipse)
}

# Function to test stress levels at different
# Usage: checkMDSdim(seagrass)
# Optional arguments include
#    iter: the number of times to run MDS for each dimension
#    dimension: a vector of dimensions to run
#    distance, try, and trymax as inputs of metaMDS
checkMDSdim = function(data,iter=1,dimensions=1:5,distance="bray",try=20,trymax=20){
  n = length(dimensions*iter)
  stress = rep(0,n)
  c = 0
  for (d in 1:n){
    temp = rep(0,iter)
    for (i in 1:iter){
      c = c + 1
      stress[c] = metaMDS(data,k=d,distance=distance,try=try,trymax=trymax)$stress
    }
  }
  output = data.frame(dimension = rep(dimensions,each=iter),
                      stress = stress)
  return(output)
}

combine_pvalues = function(p){
  return(1-pchisq(-2*sum(log(p),na.rm=T),2*sum(!is.na(p))))
}