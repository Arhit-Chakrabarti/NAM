Estimated_Plot_Function <- function(Data, Out_NAM_VI_MultiRun_Parallel){
  
  if (!require(tidyverse)) install.packages("tidyverse", dependencies = TRUE); suppressPackageStartupMessages(library(tidyverse))
  if (!require(pals)) install.packages("pals", dependencies = TRUE); suppressPackageStartupMessages(library(pals))
  if (!require(latex2exp)) install.packages("latex2exp", dependencies = TRUE); suppressPackageStartupMessages(library(latex2exp))
  
  X = Data$X
  Y = Data$Y
  GC = Out_NAM_VI_MultiRun_Parallel$est.s
  OC = Out_NAM_VI_MultiRun_Parallel$z.est
  ARI_GC = aricode::ARI(Data$s.true, GC)
  ARI_OC = 0
  for(j in 1:length(OC)){
    ARI_OC[j] = aricode::ARI(Data$z_true[[j]], OC[[j]])
  }
  
  
  
  myvalues = unname(c(kelly(n = 22),
                      alphabet2(n = (9)))[-1])
  
  names(myvalues) = 1:30
  
  plot_x.est <- data.frame(X1 = X[1,], X2 = X[2,], Group = factor(GC)) %>% ggplot(aes(x = X1, y = X2, col = Group)) + geom_point(size = 3) + scale_color_manual(values = myvalues) + 
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
    ) + guides(color = guide_legend(title = "Cluster")) + labs(x = TeX("$X_1$", bold = TRUE),
                                                               y = TeX("$X_2$", bold = TRUE),
                                                               title = "Estimated GC",
                                                               subtitle = paste0("ARI = ", round(ARI_GC, 3)))
  
  
  Y_matrix = NULL
  J = length(Y)
  for(j in 1:J){
    Y_matrix = cbind(Y_matrix, Y[[j]])
  } 
  
  Group = NULL
  Population = NULL
  Population_ARI = NULL
  n = unlist(lapply(Y, ncol))
  
  for(j in 1:J){
    group <- paste0("Group ", rep(GC[j], n[j]))
    Group = c(Group, group)
    population = paste0("Individual ", rep(j, n[j]))
    Population = c(Population, population)
    Population = c(Population, population)
    
    population_ARI = paste0("Individual ", rep(j, n[j]), "\nARI = ", round(ARI_OC[j], 3))
    Population_ARI = c(Population_ARI, population_ARI)
    
  }
  
  if (!require(stringr)) install.packages("stringr", dependencies = TRUE); suppressPackageStartupMessages(library(stringr))
  
  Y_data.frame = data.frame(Y1 = Y_matrix[1, ],
                            Y2 = Y_matrix[2, ],
                            Group = factor(Group),
                            Population = factor(Population, 
                                                levels = unique(str_sort(Population, numeric = TRUE))),
                            Population_ARI = factor(Population_ARI,
                                                    levels = unique(str_sort(Population_ARI, numeric = TRUE))),
                            Cluster = factor(as.numeric(unlist(OC))))
  
  myvalues_y = myvalues
  
  names(myvalues_y) = paste0("Group ", 1:length(myvalues_y))
  
  plot_y.est <- Y_data.frame %>% ggplot(aes(x = Y1, y = Y2, col = Group, shape = Cluster)) + geom_point(size = 2)  + scale_color_manual(values = myvalues_y) + facet_wrap(~Population_ARI) + theme_minimal() +  
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
    ) + guides(color = guide_legend(title = "Estimated GC"),
               shape = guide_legend(title = "Estimated OC")) + labs(x = TeX("$Y_1$", bold = TRUE),
                                                                    y = TeX("$Y_2$", bold = TRUE)) + 
    scale_shape_manual(values = c(1, 2, 3, 4, 5, 8, 11, 13, 15, 17, 19, 21, 23, 25))
  
  return(list(plot_x.est = plot_x.est,
              plot_y.est = plot_y.est))
}