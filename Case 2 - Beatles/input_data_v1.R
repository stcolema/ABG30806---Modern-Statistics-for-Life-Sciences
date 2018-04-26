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

prop_lm <- lm(prop_dead ~ dose, data = insect_dose)

par(mfrow = c(2, 2))
plot(prop_lm)
par(mfrow = c(1, 1))

# int <- prop_lm$coefficients[[1]]
# slope <- prop_lm$coefficients[[2]]

# lm_coeff <- summary(prop_lm)$coefficients
# 
# int <- coeff[1,1]
# slope <- coeff[2,1]
# int_se <- coeff[1,2]
# slope_se <- coeff[2,2]

# p1 <- ggplot(data = insect_dose, mapping = aes(x = dose, y = prop_dead)) + 
#   geom_point() +
#   labs(
#     # title = "Initial visualisation",
#     x = "Dose",
#     y = "Proportion of sample dead"
#     # fill = "This is the fill",
#     # caption = "This is a caption"
#   )
# 
# 
# p1 +
#   geom_abline(intercept = int, slope = slope, colour = "red", size = 1.2) + 
#   labs(title = "Model fitted line")
# 
# p1 + 
#   geom_smooth(method = "loess", se = FALSE) + 
#   labs(title = "LOESS fitted line")

predicted <- prop_lm$fitted.values
p_lm <- qplot(data = insect_dose, x = dose, y = prop_dead)+
  geom_line(y = predicted, size = 1, colour = "blue")+
  geom_segment(aes(xend = dose, yend = predicted, color = "error"))+
  labs(title ="LM regression errors",
       x = "Dose",
       y = "Proportion of sample dead",
       color ="series")

loess_mod <- loess(prop_dead ~ dose, data = insect_dose)

p_loess <- ggplot(data = insect_dose, aes(x = dose, y = prop_dead)) + 
  geom_point() + 
  geom_smooth(se=FALSE, method = "loess") +
  geom_segment(aes(xend = dose, yend = fitted(loess_mod),  color = "error")) +
  labs(title ="LOESS regression errors",
       x = "Dose",
       y = "Proportion of sample dead",
       color ="series")

multiplot(p_lm, p_loess, cols = 1)
