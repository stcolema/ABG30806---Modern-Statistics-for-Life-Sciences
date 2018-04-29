#!/usr/bin/env Rscript

library(tidyverse)

# Multiple plot function from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots <- length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots / cols)),
                     ncol = cols, nrow = ceiling(numPlots / cols)
    )
  }
  
  if (numPlots == 1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(
        layout.pos.row = matchidx$row,
        layout.pos.col = matchidx$col
      ))
    }
  }
}


insect_dose <- tribble(
  ~ dose, ~ num_insects, ~ num_died,
  #-----| ------------| ----------|
  1.6907,           59,          6,
  1.7242,           60,         13,
  1.7552,           62,         18,
  1.7842,           56,         28,
  1.8113,           63,         52,
  1.8369,           59,         53,
  1.8610,           62,         61,
  1.8839,           60,         60
)

# Add proportion variable
insect_dose <- insect_dose %>%
  mutate(prop_dead = num_died / num_insects)

# Linear regression
prop_lm <- lm(prop_dead ~ dose, data = insect_dose)

# Ugh, does not look good
par(mfrow = c(2, 2))
plot(prop_lm)
par(mfrow = c(1, 1))

# Create a plot of fitted values and errors
predicted <- prop_lm$fitted.values
p_lm <- qplot(data = insect_dose, x = dose, y = prop_dead)+
  geom_line(y = predicted, size = 1, colour = "blue")+
  geom_segment(aes(xend = dose, yend = predicted, color = "error"))+
  labs(title ="LM regression",
       x = "Dose",
       y = "Proportion of sample dead",
       color ="series")

# Compare to Local Polynomial Regression
# (overfitting, but should give better idea of ideal curve)
loess_mod <- loess(prop_dead ~ dose, data = insect_dose)

p_loess <- ggplot(data = insect_dose, aes(x = dose, y = prop_dead)) + 
  geom_point() + 
  geom_smooth(se=FALSE, method = "loess") +
  geom_segment(aes(xend = dose, yend = fitted(loess_mod),  color = "error")) +
  labs(title ="LOESS regression",
       x = "Dose",
       y = "Proportion of sample dead",
       color ="series")

# Compare plots; LOESS reveals sigmoidal relationship
# Logistic regression or Cloglog or Probit then
multiplot(p_lm, p_loess, cols = 1)


prop_log_reg <- glm(prop_dead ~ dose, family = binomial, data = insect_dose)
summary(prop_log_reg)


plot_fitted_errors <- function(model, data, x, y){
  predictions <- fitted(model)
  p <- qplot(data = data, x = x, y = y)+
    geom_line(y = predictions, size = 1, colour = "blue")+
    geom_segment(aes(xend = x, yend = predictions, color = "error"))
    # labs(title ="LM regression",
    #      x = "Dose",
    #      y = "Proportion of sample dead",
    #      color ="series")
  return(p)
}

plot_fitted_errors(prop_log_reg, insect_dose, "dose", "prop_dead")

