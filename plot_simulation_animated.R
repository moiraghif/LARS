## RUN THIS CODE AFTER simulation.R

library(tidyverse)
library(gganimate)


#Function to create db
create_db <- function(x) {
  db <- data.frame(cbind(1:length(x$err_log), x$err_log, x$log, x$variance_log, x$bias_log))
  names(db) <- c("Iter", "Errors", "Bias² + Variance", "Variance", "Bias²")
  return(db)
}

#db <- create_db(lars_simulation_fixed)
db <- create_db(lars_simulation_random)

#Arrow plot for Variance
if(is.na(db$Variance[1])) db$Variance[1] <- 0
arrow <- c()
for(k in 2:monte_carlo) {
  if(db$Variance[k]<db$Variance[k-1])
      arrow <- c(arrow,"\u2193")
  else if (db$Variance[k]==db$Variance[k-1])
    arrow <- c(arrow,"\u2192")
  else
    arrow <- c(arrow,"\u2191")
}

#Arrow plot for Bias
arrow_bias <- c()
for(k in 2:monte_carlo) {
  if(db$"Bias²"[k]<db$"Bias²"[k-1])
      arrow_bias <- c(arrow_bias,"\u2193")
  else if (db$"Bias²"[k]==db$"Bias²"[k-1])
      arrow_bias <- c(arrow_bias,"\u2192")
  else
      arrow_bias <- c(arrow_bias,"\u2191")
}

##Convergence animated Plot
z <- seq(1,nrow(db))
cols <- c("Bias² + Variance"="#00AFBB","Estimates Fluctuation"="#FC4E07","Expected Value"="#C4961A")
conv <- ggplot() +
  geom_hline(yintercept = mean(db$Errors),linetype="dashed")+
  #Theoretical Error
  geom_point(data = db, aes(x = db$Iter, y = db$`Bias² + Variance`, colour="Bias² + Variance"), size=6) +
  geom_line(data = db, aes(x = db$Iter, y = db$`Bias² + Variance`, colour="Bias² + Variance"), size=0.1) +
  #Calculated Error
  geom_point(data = db, aes(x = db$Iter, y = db$Errors, colour="Estimates Fluctuation"), size=6) +
  geom_line(data = db, aes(x = db$Iter, y = db$Errors, colour="Estimates Fluctuation"), size=0.1) +
  labs(title = "Convergence to the Expected MSE (Fixed-X Setting)")+
  ylim(0,5000)+
  xlim(0,length(z))+
  xlab("Iteration")+
  ylab("Value") +
  #Final points (Averages)
  geom_point(aes(x=length(z),y=mean(db$Errors),colour="Expected Value"), size=6)+
  geom_text(data=data.frame(z=z),
            mapping = aes(x = 0.8*length(z), y = 4700,
                          label = paste0("ITERATION N.  ",z,"\n","E(Bias²):       ",as.integer(db$"Bias²")," ",arrow_bias,"\n",
                                         "E(Variance):    ",as.integer(db$Variance)," ",arrow)))+
  theme(plot.title = element_text(size = 20, face = "bold"),
        legend.title=element_text(size=15, face = "bold"), 
        legend.text=element_text(size=15),
        plot.caption = element_text(face = "italic", hjust = 0.5, size = 15))+
  theme_bw()+
  labs(colour="MSE Values")+
  transition_reveal(z)+
  ease_aes("linear")+
  enter_appear()
plot1 <- animate(conv, fps=10, end_pause = 10)
plot1
#Save
anim_save("/Users/riccardocervero/Desktop/Plot1.gif",plot1)

##E(Bias) and E(Variance)
z <- seq(1,nrow(db))
cols <- c("E(Bias²)"="#00AFBB","E(Variance)"="#C4961A")
BiasVariance <- ggplot() +
  #Bias
  geom_line(data = db, aes(x = db$Iter, y = db$`Bias²`, colour="E(Bias²)"), size=0.5) +
  #Variance
  geom_line(data = db, aes(x = db$Iter, y = db$Variance, colour="E(Variance)"), size=0.5) +
  labs(title = "Bias² and Variance Trend")+
  ylim(0,1.01*max(db))+
  xlim(0,length(z))+
  xlab("Iteration")+
  ylab("Value") +
  geom_text(data=data.frame(z=z),
            mapping = aes(x = 0.8*length(z), y = 1.01*max(db),
                          label = paste0("ITERATION N.  ",z)))+
  theme(plot.title = element_text(size = 25, face = "bold"),
        legend.title=element_text(size=20, face = "bold"), 
        legend.text=element_text(size=20),
        plot.caption = element_text(face = "italic", size = 15))+
  theme_bw()+
  labs(colour="Value")+
  transition_reveal(z)+
  ease_aes("linear")+
  enter_appear()
plot2 <- animate(BiasVariance, fps=10, end_pause = 10)
plot2
#Save
anim_save("/Users/riccardocervero/Desktop/Plot2.gif",plot2)
