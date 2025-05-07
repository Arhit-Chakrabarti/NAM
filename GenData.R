# Function to Generate Data with equal sample sizes
Generate_Data <- function(T.true, alpha0_x, J, q, lambda_x, L.true, p, alpha0_y, n, lambda_y, seed.x, seed.y, plot = TRUE){
  # Generate the data
  if (!require(extraDistr)) install.packages("extraDistr", dependencies = TRUE); suppressPackageStartupMessages(library(extraDistr))
  if (!require(MASS)) install.packages("MASS", dependencies = TRUE); suppressPackageStartupMessages(library(MASS))
  # True weights
  alpha_x.true = rgamma(n = 1, shape = alpha0_x, scale = 1)
  pi_star.true = as.numeric(rdirichlet(n = 1, alpha = rep(alpha_x.true/T.true, T.true))) # True beta
  
  set.seed(seed.x)
  
  s.true = sample(1:T.true, size = J, prob = pi_star.true, replace = TRUE)
  lambda_x.inv = 1/lambda_x
  
  mu_x = rep(0, q)
  sigma_x.true = stats::rWishart(n = T.true, df = q+5, Sigma = diag(1, q))
  mean_x.true = sapply(1:T.true, function(j){MASS::mvrnorm(n = 1, mu = mu_x, Sigma = lambda_x.inv * solve(sigma_x.true[, ,j]))}) # True mean
  
  X = sapply(1:J, function(j){
    MASS::mvrnorm(n = 1, mu = mean_x.true[ , s.true[j]], Sigma = solve(sigma_x.true[, ,s.true[j]]))
  }) # Covariate data
  
  set.seed(seed.y)
  alpha_y.true = rgamma(n = 1, shape = alpha0_y, scale = 1)
  
  S_t <- list()
  
  for(t in 1:T.true){
    S_t[[t]] <- which(s.true == t)
  }
  
  pi_true <- list()
  
  for(t in 1:length(S_t)){
    pi_true[[t]] <- rdirichlet(n = 1, alpha = rep(alpha_y.true/L.true, L.true))
  }
  
  pi_true_matrix <- NULL
  
  for(t in 1:length(S_t)){
    pi_true_matrix <- rbind(pi_true_matrix, pi_true[[t]])
  }
  
  z_true <- replicate(J, list(0))
  
  for(t in 1:length(pi_true)){
    for(j in S_t[[t]]){
      for(i in 1:n[j]){
        z_true[[j]][i] = extraDistr::rcat(n = 1, prob = pi_true_matrix[t, ])
      }
    }
  }
  
  mu_y = rep(0, p)
  lambda_y.inv = 1/lambda_y
  
  sigma_y.true = stats::rWishart(n = L.true, df = p+5, Sigma = diag(1, p))
  mean_y.true = sapply(1:L.true, function(j){MASS::mvrnorm(n = 1, mu = mu_y, Sigma = lambda_y.inv * solve(sigma_y.true[ , , j]))}) # True mean
  
  Y = lapply(1:J, function(j) matrix(0, nrow = p, ncol = n[j])) 
  
  for(j in 1:J){
    for(i in 1:n[j]){
      Y[[j]][, i] <- MASS::mvrnorm(n = 1, mu = mean_y.true[  ,z_true[[j]][i]], Sigma = solve(sigma_y.true[, , z_true[[j]][i]]))
    }
  }
  
  
  Y_matrix = NULL
  
  for(j in 1:J){
    Y_matrix = cbind(Y_matrix, Y[[j]])
  } 
  
  Group = NULL
  Population = NULL
  
  for(j in 1:J){
    group <- paste0("Group ", rep(s.true[j], n[j]))
    Group = c(Group, group)
    population = paste0("Individual ", rep(j, n[j]))
    Population = c(Population, population)
  }
  
  if (!require(stringr)) install.packages("stringr", dependencies = TRUE); suppressPackageStartupMessages(library(stringr))
  
  Y_data.frame = data.frame(Y1 = Y_matrix[1, ],
                            Y2 = Y_matrix[2, ],
                            Group = factor(Group),
                            Population = factor(Population, 
                                                levels = unique(str_sort(Population, numeric = TRUE))),
                            Cluster = factor(as.numeric(unlist(z_true))))
  
  if(plot == TRUE){
    if (!require(tidyverse)) install.packages("tidyverse", dependencies = TRUE); suppressPackageStartupMessages(library(tidyverse))
    if (!require(pals)) install.packages("pals", dependencies = TRUE); suppressPackageStartupMessages(library(pals))
    if (!require(latex2exp)) install.packages("latex2exp", dependencies = TRUE); suppressPackageStartupMessages(library(latex2exp))
    
    myvalues = unname(c(kelly(n = 22),
                        alphabet2(n = (9)))[-1])
    
    names(myvalues) = 1:30
    
    plot_x <- data.frame(X1 = X[1,], X2 = X[2,], Group = factor(s.true)) %>% ggplot(aes(x = X1, y = X2, col = Group)) + geom_point(size = 3) + scale_color_manual(values = myvalues) + 
      theme_minimal() +  
      theme(
        # LABLES APPEARANCE
        # panel.grid.major = element_blank(), 
        # panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        plot.title = element_text(hjust = 0.5, size=20, face= "bold", colour= "black" ),
        plot.subtitle = element_text(hjust = 0.5, size=16, face= "bold", colour= "black" ),
        axis.title.x = element_text(size=16, face="bold", colour = "black"),    
        axis.title.y = element_text(size=16, face="bold", colour = "black"),    
        axis.text.x = element_text(size=16, face="bold", colour = "black"), 
        axis.text.y = element_text(size=16, face="bold", colour = "black"),
        strip.text.x = element_text(size = 14, face="bold", colour = "black" ),
        strip.text.y = element_text(size = 14, face="bold", colour = "black"),
        axis.line.x = element_line(color="black", size = 0.3),
        axis.line.y = element_line(color="black", size = 0.3),
        panel.border = element_rect(colour = "black", fill=NA, size=0.3), 
        legend.title=element_text(size=16),
        legend.text=element_text(size=14)
      ) + guides(color = guide_legend(title = "True GC")) + labs(x = TeX("$X_1$", bold = TRUE),
                                                                 y = TeX("$X_2$", bold = TRUE))
    
    myvalues_y = myvalues
    names(myvalues_y) = paste0("Group ", 1:length(myvalues_y))
    
    plot_y <- Y_data.frame %>% ggplot(aes(x = Y1, y = Y2, col = Group, shape = Cluster)) + geom_point(size = 2)  + scale_color_manual(values = myvalues_y) + facet_wrap(~Population) + theme_minimal() +  
      theme(
        # LABLES APPEARANCE
        # panel.grid.major = element_blank(), 
        # panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        plot.title = element_text(hjust = 0.5, size=20, face= "bold", colour= "black" ),
        plot.subtitle = element_text(hjust = 0.5, size=16, face= "bold", colour= "black" ),
        axis.title.x = element_text(size=16, face="bold", colour = "black"),    
        axis.title.y = element_text(size=16, face="bold", colour = "black"),    
        axis.text.x = element_text(size=16, face="bold", colour = "black"), 
        axis.text.y = element_text(size=16, face="bold", colour = "black"),
        strip.text.x = element_text(size = 14, face="bold", colour = "black" ),
        strip.text.y = element_text(size = 14, face="bold", colour = "black"),
        axis.line.x = element_line(color="black", size = 0.3),
        axis.line.y = element_line(color="black", size = 0.3),
        panel.border = element_rect(colour = "black", fill=NA, size=0.3), 
        legend.title=element_text(size=16),
        legend.text=element_text(size=14)
      ) + guides(color = guide_legend(title = "True GC"),
                 shape = guide_legend(title = "True OC")) + labs(x = TeX("$Y_1$", bold = TRUE),
                                                                 y = TeX("$Y_2$", bold = TRUE)) + 
      scale_shape_manual(values = c(1, 2, 3, 4, 5, 8, 11, 13))
    
    return(list(J = J, n = n, 
                pi_star.true = pi_star.true, pi_true = pi_true, pi_true_matrix = pi_true_matrix, 
                alpha_x.true = alpha_x.true, alpha_y.true = alpha_y.true,
                
                s.true = s.true, 
                z_true = z_true,
                mean_x.true = mean_x.true, sigma_x.true = sigma_x.true,
                mean_y.true = mean_y.true, sigma_y.true = sigma_y.true,
                X = X, Y = Y,
                plot_x = plot_x, plot_y = plot_y))
    
  }else{
    return(list(J = J, n = n, 
                pi_star.true = pi_star.true, pi_true = pi_true, pi_true_matrix = pi_true_matrix, 
                alpha_x.true = alpha_x.true, alpha_y.true = alpha_y.true,
                
                s.true = s.true, 
                z_true = z_true,
                mean_x.true = mean_x.true, sigma_x.true = sigma_x.true,
                mean_y.true = mean_y.true, sigma_y.true = sigma_y.true,
                X = X, Y = Y))
  }
}