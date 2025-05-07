source("NAM_VI.R") # Load the NAM VI function
# Load the necessary library for parallel execution
if (!require(parallel)) install.packages("parallel", dependencies = TRUE); suppressPackageStartupMessages(library(parallel))
# Define the number of cores to use
numberOfCores <- detectCores() - 1  # Use one less than the total cores

# Run the function in parallel
Run_NAM_VI_MultiRun_Parallel <- function(numReplicate, X, Y, H, L, maxIter, epsilon){
  out <- mclapply(1:numReplicate, function(i) {
    Run_NAM_VI(X = X, Y = Y,
               H = H, L = L,
               maxIter = maxIter, epsilon = epsilon)}
    , mc.cores = numberOfCores)
  
  
  ELBO = 0 # Define ELBO values at convergence to store
  for(i in 1:length(out)){
    ELBO[i] = out[[i]]$ELBO_val_conv
  }
  
  # Choose the best output based on maximum ELBO value at 
  best_outputs = out[[which.max(ELBO)]]
  return(best_outputs)
}