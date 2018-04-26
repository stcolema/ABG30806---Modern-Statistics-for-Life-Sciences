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

# Functions -----------------------------------------------------------------

# Function to return a data frame of the results of the shapiro wilkes test
# applied to the columns of the input data frame
find_variable_normality <- function(numerical_data,
                                    sig_threshold = 0.05,
                                    multiple_testing_method = "BH") {

  ## Perform the Shapiro-Wilk test of normality
  lshap <- lapply(numerical_data, shapiro.test)

  # Put p-values of normality test into data frame
  df <- data.frame("Condition" = character(), "P.value" = double())

  for (i in 1:length(lshap)) {
    entry <- data.frame("Condition" = names(lshap)[i], "P.value" = lshap[[i]]$p.value)
    df <- rbind(df, entry)
  }

  # Use BH method to correct p-values for multiple testing
  df$P.Adjust <- p.adjust(df$P.value, multiple_testing_method)
  df <- df[order(df$P.Adjust), ]

  # Any entries with p adjusted values below 0.05 (pre-defined threshold)
  # are considered significantly normal
  df$Significant <- 0
  df$Significant[df$P.Adjust < sig_threshold] <- 1

  return(df)
}

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

# Check normality of measurement variables
df_norm <- find_variable_normality(measurement_data, sig_threshold = sig_threshold)

# Seperate out non-normal conditions
non.normal.conditions <- df_norm[df_norm$Significant == 1, "Condition"]
conditions <- colnames(measurement_data)
conditions.to.investigate <- conditions[!conditions %in% non.normal.conditions]

qqPlot(measurement_data$diam_stem_mm)
qqPlot(measurement_data$leafarea_cm2)
qqPlot(measurement_data$firmness_log_mm)

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
anova(stem_firmness_lm)

# Similarly for leaf area and firmness
# leafarea_stem_lm <- lm(leafarea_cm2 ~ diam_stem_mm, data = tomato_data)
leafarea_firmness_lm <- lm(leafarea_cm2 ~ firmness_log_mm, data = tomato_data)

# F-test
# anova(leafarea_stem_lm) # result is symmetric
anova(leafarea_firmness_lm)

ols_plot_cooksd_bar(leafarea_firmness_lm)
ols_plot_resid_stud(leafarea_firmness_lm)

# All three measurement variables are related

# Model without type (including interaction)
stem_model_no_type <- glm(diam_stem_mm ~ leafarea_cm2 * firmness_log_mm,
  family = gaussian,
  data = tomato_data
)

# Model with type
stem_model_type <- glm(diam_stem_mm ~ leafarea_cm2 * firmness_log_mm * type.,
  family = gaussian,
  data = tomato_data
)

# Check significance of predictors using type II ANOVA
Anova(stem_model_no_type, type = "II")
Anova(stem_model_type, type = "II")

# Inspect models
summary(stem_model_no_type)
summary(stem_model_type)

# Check if deviance is in expected range
pchisq(stem_model_no_type$deviance, stem_model_no_type$df.residual, lower.tail = T)
pchisq(stem_model_type$deviance, stem_model_type$df.residual, lower.tail = T)

# Leaf area only important when type is uknown

# Investigate relationships between firmness and the other variables

# Model without type (including interaction)
firmness_model_no_type <- glm(firmness_log_mm ~ leafarea_cm2 * diam_stem_mm,
  family = gaussian,
  data = tomato_data
)

# Model with type
firmness_model_type <- glm(firmness_log_mm ~ leafarea_cm2 * diam_stem_mm * type.,
  family = gaussian,
  data = tomato_data
)

# Check significance of predictors using type II ANOVA
Anova(firmness_model_no_type, type = "II")
Anova(firmness_model_type, type = "II")

# Inspect models
summary(firmness_model_no_type)
summary(firmness_model_type)

# Check if deviance is in expected range
# CHECK WITH CHRIS re lower.tail
pchisq(firmness_model_no_type$deviance, firmness_model_no_type$df.residual, lower.tail = T)
pchisq(firmness_model_type$deviance, firmness_model_type$df.residual, lower.tail = T)

# Both variables are considered significant even accounting for type; however the
# deviance is out of whack

# Relationship between firmness and leaf area dependent upon type? ---------

# build the saturated model (different slopes and intercepts)
model_sat <- glm(leafarea_cm2 ~ firmness_log_mm * type.,
  family = gaussian,
  data = tomato_data
)

# build the model with different intercepts but common slope
model_int <- lm(leafarea_cm2 ~ firmness_log_mm + type.,
  data = tomato_data
)

Anova(model_int, type = "II")
emm <- emmeans(model_int, ~ type.)
cld(emm, alpha = 0.05, adjust = "tukey")
plot(emm)

# model with different intercepts and slope of 0
model_int_slope_0 <- glm(leafarea_cm2 ~ type.,
  family = gaussian,
  data = tomato_data
)

# model with common intercept and slope
model_int_slope <- glm(leafarea_cm2 ~ type.,
  family = gaussian,
  data = tomato_data
)

# Compare models using the likelihood ratio test
lrtest(model_sat, model_int) # insignificant difference
lrtest(model_int, model_int_slope_0) # significant difference
lrtest(model_int, model_int_slope) # significant difference

# Can claim that leaf area relationship with firmness does
# depend on the type (different intercepts for each type)
# but the the realtive impact of a change in firmness upon
# leaf area is common to all types (common slope)

# The optimal model for describing the leaf area includes the
# type and firmness as features but not their interaction
final_model <- model_int
summary(final_model) # the deviance is awful. really bad model. no relationship?


# Compare with model including stem diameter
leaf_area_model_diam <- glm(leafarea_cm2 ~ firmness_log_mm + type. + diam_stem_mm,
  family = gaussian,
  data = tomato_data
)

summary(leaf_area_model_diam)
Anova(leaf_area_model_diam, type = "II")
lrtest(final_model, leaf_area_model_diam)

# Between ANOVA and likelihood-ratio test diameter is not relvant (little commonality)




# Jarno --------------------------------------------------

# Model with all terms
stem_anova_wo_int <- lm(diam_stem_mm ~ firmness_log_mm + leafarea_cm2 + type.,
  data = tomato_data
)

stem_anova_w_int <- lm(diam_stem_mm ~ firmness_log_mm + leafarea_cm2 + type. + leafarea_cm2 * firmness_log_mm,
  data = tomato_data
)

firmness_anova_w_int <- lm(firmness_log_mm ~ diam_stem_mm + leafarea_cm2 + type.,
  data = tomato_data
)

firmness_anova_w_int_2 <- lm(firmness_log_mm ~ type. + diam_stem_mm + leafarea_cm2 + diam_stem_mm * leafarea_cm2,
  data = tomato_data
)

firmness_anova_wo_type <- lm(firmness_log_mm ~ diam_stem_mm + leafarea_cm2,
  data = tomato_data
)

anova(firmness_anova_w_int, firmness_anova_w_int_2)

plot(firmness_anova_w_int_2)

Anova(firmness_anova_wo_type, type = "II")
plot(firmness_anova_wo_type)
summary(firmness_anova_w_int)

# Check for groups
emm <- emmeans(firmness_anova_w_int, ~ type.)
cld(emm, alpha = 0.05, adjust = "tukey")

### ekf ###
x <- model.matrix(firmness_anova_w_int_2)[, -1]
y <- tomato_data$diam_stem_mm
adjr <- leaps(x, y, method = "adjr2")

maxadjr(adjr, 8)

anova(stem_anova_wo_int, stem_anova_w_int)
Anova(stem_anova_wo_int, type = "II")

# Firmness * Area
anov_tom_typ_no <- lm(diam_stem_mm ~ type. + leafarea_cm2 + firmness_log_mm,
  data = tomato_data
)

emm <- emmeans(anov_tom_typ_no, ~ type.)
cld(emm, alpha = 0.05, adjust = "tukey")
