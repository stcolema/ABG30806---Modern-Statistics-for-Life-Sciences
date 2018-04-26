#!/usr/bin/env Rscript

# install.packages(c("tidyverse", "gdata", "lmtest", "car", "GGally"), dep=T)

# Look into using olsrr for diagnosing models (DFBETAs and Cook's Distance), see
# https://cran.r-project.org/web/packages/olsrr/vignettes/influence_measures.html

# Also use plot(compareFits(coef(model_1), coef(model_2))) for comparing the weights
# attributed to variables in different models

# Use plot( augPred(model), aspect = "xy", grid = T ) to plot predicted
# versus actual values

# Call the requisite packackages
require(tidyverse)
require(gdata)
library(olsrr)
library(GGally)
library(emmeans)
require(lmtest)
require(car)
library(nlme)
library(leaps)
library(faraway)
library(sjmisc)
library(ggfortify)

# Functions -----------------------------------------------------------------

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


# Main program -------------------------------------------

# Set to working directory
my_wd <- paste0(
  "C:/Users/steph/Desktop/Bioinformatics/ABG30806",
  " - Modern Statistics for Life Sciences/Case 1 ",
  "(Chris)"
)

setwd(my_wd)

# Read in the xls file containing the data
datafile_name <- "TomatoesChris.xls"
tomato_data <- read.xls(datafile_name, sheet = 1, header = TRUE)

# Set the significance threshold for tests
sig_threshold <- 0.05

# Inspect the data (check for missing values or irregularities)
summary(tomato_data)

# Set the PlantNr to the row name and remove as a descriptive variale
rownames(tomato_data) <- tomato_data$PlantNr
tomato_data <- select(tomato_data, -PlantNr)

# Exploratory data analysis -----------------------------------------------

# Inspect the edited data frame
summary(tomato_data)

# Check if the design is balanced
table(tomato_data$type.)

# Pairwise plots
ggpairs(tomato_data)

# Seperate the measurement data into a new data frame (purely numeric)
measurement_data <- tomato_data %>%
  dplyr::select(-type.)

# Inspect the data for each unique type
levels(tomato_data$type.)

# Investigate relationship between measurements based on type
# Scatter plots
ggplot(tomato_data) +
  geom_point(mapping = aes(x = diam_stem_mm, y = leafarea_cm2, colour = type.))

ggplot(tomato_data) +
  geom_point(mapping = aes(x = leafarea_cm2, y = firmness_log_mm, colour = type.))

ggplot(tomato_data) +
  geom_point(mapping = aes(x = firmness_log_mm, y = leafarea_cm2, colour = type.))

# Plot the densities of the measurements based on type
p1 <- ggplot(tomato_data) +
  geom_density(mapping = aes(x = diam_stem_mm, colour = type.))

p2 <- ggplot(tomato_data) +
  geom_density(mapping = aes(x = leafarea_cm2, colour = type.))

p3 <- ggplot(tomato_data) +
  geom_density(mapping = aes(x = firmness_log_mm, colour = type.))

multiplot(p1, p2, p3, cols = 2)

# We know data is unbalanced and that the variables are normally distributed (ish)

# Investigating across types

# Inspect the data for each unique type
levels(tomato_data$type.)

# Seperate data frames for each type
b_tomato <- tomato_data %>%
  dplyr::filter(type. == "b") %>%
  dplyr::select(-type.)

c_tomato <- tomato_data %>%
  dplyr::filter(type. == "c") %>%
  dplyr::select(-type.)

r_tomato <- tomato_data %>%
  dplyr::filter(type. == "r") %>%
  dplyr::select(-type.)

# See that the relationships between measurements vary significantly
# across the types (note the negative correlation between stem diameter
# and leaf area in cherry tomatoes, and the complete lack of relationship
# between firmness and the other variables in cherry tomatoes).
cor(b_tomato)
cor(c_tomato)
cor(r_tomato)

ggpairs(c_tomato)
ggpairs(r_tomato)
ggpairs(b_tomato)

# Interaction plots and type specific model comparison -------------------------------

# Draw histogram of each variable
# Decide on bins
# Create factors for interaction plot
# lm across type, check overlap in intervals

# Create a df for playing with categorical variables
cat_df <- tomato_data

# Investigate distribution of each measurement variable
hist(tomato_data$diam_stem_mm, xlab = "Stem diameter in mm", main = "", breaks = 5)
hist(tomato_data$leafarea_cm2, xlab = "Leaf area in cm2", main = "", breaks = 5)
hist(tomato_data$firmness_log_mm, xlab = "Log(firmness in mm)", main = "", breaks = 5)

# Create categorical variables for measurements
# (6 is pretty arbitrary - probably too large)
cat_df$stem_cat <- as.ordered(ntile(tomato_data$diam_stem_mm, 6))
cat_df$leaf_area_cat <- as.ordered(ntile(tomato_data$leafarea_cm2, 6))
cat_df$firmness_cat <- as.ordered(ntile(tomato_data$firmness_log_mm, 6))

# Interaciton plots
interaction.plot(cat_df$leaf_area_cat, cat_df$type., cat_df$diam_stem_mm)
interaction.plot(cat_df$leaf_area_cat, cat_df$type., cat_df$firmness_log_mm)
interaction.plot(cat_df$stem_cat, cat_df$type., cat_df$firmness_log_mm)

# Try linear models across type
stem_lm <- nlme::lmList(diam_stem_mm ~ firmness_log_mm * leafarea_cm2 | type., tomato_data)
leafarea_lm <- nlme::lmList(leafarea_cm2 ~ firmness_log_mm * diam_stem_mm | type., tomato_data)

# Lets look at the overlap in models
plot(intervals(stem_lm))
plot(intervals(leafarea_lm))

# Notice cherry's far larger error bars (due to unbalanced nature of experiment)


# Are there differences for the measurement types -------------------------

# Some boxplots of the variables
box_p1 <- ggplot(tomato_data) +
  geom_boxplot(mapping = aes(x = type., y = diam_stem_mm))

box_p2 <- ggplot(tomato_data) +
  geom_boxplot(mapping = aes(x = type., y = leafarea_cm2))

box_p3 <- ggplot(tomato_data) +
  geom_boxplot(mapping = aes(x = type., y = firmness_log_mm))

multiplot(box_p1, box_p2, box_p3, cols = 2)
# Note how type c tends to have less / no overlap with b and r

## We showed here that there are differences ##
mod <- lm(diam_stem_mm ~ type., data = tomato_data)
anova(mod)
summary(lm(diam_stem_mm ~ type., data = tomato_data))

par(mfrow = c(2, 2))
plot(lm(diam_stem_mm ~ type., data = tomato_data))
summary(lm(firmness_log_mm ~ type., data = tomato_data))
plot(lm(firmness_log_mm ~ type., data = tomato_data))
summary(lm(leafarea_cm2 ~ type., data = tomato_data))
plot(lm(leafarea_cm2 ~ type., data = tomato_data))
par(mfrow = c(1, 1))


# Relationship between measurement variables ----------------------------

# First check relationship between stem diameter and other variables
stem_leafarea_lm <- lm(diam_stem_mm ~ leafarea_cm2, data = tomato_data)
stem_firmness_lm <- lm(diam_stem_mm ~ firmness_log_mm, data = tomato_data)

# Inspect models
ols_plot_cooksd_bar(stem_leafarea_lm)
ols_plot_resid_stud(stem_leafarea_lm)

ols_plot_cooksd_bar(stem_firmness_lm)
ols_plot_resid_stud(stem_firmness_lm)

# F-test
anova(stem_leafarea_lm)
summary(stem_leafarea_lm)

anova(stem_firmness_lm)
summary(stem_firmness_lm)

# Similarly for leaf area and firmness
# leafarea_stem_lm <- lm(leafarea_cm2 ~ diam_stem_mm, data = tomato_data)
leafarea_firmness_lm <- lm(leafarea_cm2 ~ firmness_log_mm, data = tomato_data)

# F-test
# anova(leafarea_stem_lm) # result is symmetric
anova(leafarea_firmness_lm)

ols_plot_cooksd_bar(leafarea_firmness_lm)
ols_plot_resid_stud(leafarea_firmness_lm)

# All three measurement variables are related


# Relationship between firmness and leaf area dependent upon type? ---------

# build the model with different intercepts but common slope
leafarea_firmness_lm <- lm(leafarea_cm2 ~ firmness_log_mm + type.,
  data = tomato_data
)

Anova(leafarea_firmness_lm, type = "II")
emm <- emmeans(leafarea_firmness_lm, ~ type.)
cld(emm, alpha = 0.05, adjust = "tukey")
plot(emm)

# Model with different slopes for each type (for comparison)
leafarea_firmness_lm_slopes <- lm(leafarea_cm2 ~ firmness_log_mm * type.,
  data = tomato_data
)

Anova(leafarea_firmness_lm_slopes, type = "II")

# Plot the predicted lines for the model
c2 <- coef(leafarea_firmness_lm)
plot(leafarea_cm2 ~ firmness_log_mm,
  col = c(1, 2, 3)[type.],
  data = tomato_data,
  xlab = "Firmness (Log) [mm]",
  ylab = "Leaf Area [cm2]"
)
coef_b <- c(c2[1], c2[2]) # ic and slope of AL
coef_c <- c(c2[1] + c2[3], c2[2]) # ic and slope of GL
coef_r <- c(c2[1] + c2[4], c2[2]) # ic and slope of GL

abline(coef_b, col = 1)
abline(coef_c, col = 2)
abline(coef_r, col = 3)
title("Model Predictions by Type")
legend("topright", title = "Type", c("b", "c", "r"), fill = c(1, 2, 3))

# Inspect model residuals and data leverage
ols_plot_cooksd_bar(leafarea_firmness_lm)
ols_plot_resid_stud(leafarea_firmness_lm)

# QQ plot to test normality
autoplot(leafarea_firmness_lm, which = 2)
# plot(leafarea_firmness_lm)[[2]]
