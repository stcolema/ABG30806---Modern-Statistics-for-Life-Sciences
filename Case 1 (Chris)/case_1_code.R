#!/usr/bin/env Rscript

# install.packages(c("tidyverse", "gdata", "lmtest", "car"), dep=T)

# Call the requisite packackages
require(tidyverse)
require(gdata)
require(lmtest)
require(car)

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

# Inspect the edited data frame
summary(tomato_data)

# Are there differences for the measurement types -----------

# Seperate the measurement data into a new data frame (purely numeric)
measurement_data <- tomato_data %>%
  dplyr::select(-type.)


# t_test_data <- data.frame(t_stem_diam_leaf_area = t.test(measurement_data$diam_stem_mm,
#                                             measurement_data$leafarea_cm2)$p.value,
#              t_stem_diam_firmness = t.test(measurement_data$diam_stem_mm,
#                                             measurement_data$firmness_log_mm)$p.value,
#              t_leaf_area_firmness = t.test(measurement_data$leafarea_cm2,
#                                             measurement_data$firmness_log_mm)$p.value)
#
#
#

# Inspect the data for each unique type
levels(tomato_data$type.)

# Seperate data frames for each type
b_tomato <- tomato_data %>%
  dplyr::filter(type. == "b")

c_tomato <- tomato_data %>%
  dplyr::filter(type. == "c")

r_tomato <- tomato_data %>%
  dplyr::filter(type. == "r")

# Investigate relationship between measurements based on type
ggplot(tomato_data) +
  geom_point(mapping = aes(x = diam_stem_mm, y = leafarea_cm2, colour = type.))

ggplot(tomato_data) +
  geom_point(mapping = aes(x = diam_stem_mm, y = firmness_log_mm, colour = type.))

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

df_norm <- find_variable_normality(measurement_data, sig_threshold = sig_threshold)

# Seperate out non-normal conditions
non.normal.conditions <- df_norm[df_norm$Significant == 1, "Condition"]
conditions <- colnames(measurement_data)
conditions.to.investigate <- conditions[!conditions %in% non.normal.conditions]


# Investigate the normality of the different conditions
# In the following plots the data should form a straight line (the line formed by the qqline
# command)
qqnorm(measurement_data[, 1])
qqline(measurement_data[, 1])

## Have a look at the densities
plot(density(measurement_data[, 1]))
plot(density(measurement_data[, 2]))
plot(density(measurement_data[, 3]))

# They're grand. Shapiro pulled it off

# Check means
t.test(measurement_data[, 1], tomato_data[, 2])
t.test(tomato_data[, 1], tomato_data[, 3])
t.test(tomato_data[, 2], tomato_data[, 3])

names(measurement_data)

diam_lm <- glm(diam_stem_mm ~ leafarea_cm2 * firmness_log_mm,
  family = gaussian,
  data = measurement_data
)

Anova(diam_lm, type = "II")

par(mfrow = c(2, 2))
plot(diam_lm)
par(mfrow = c(1, 1))

# remove_outliers <- function(x, na.rm = TRUE, ...) {
#   qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
#   H <- 1.5 * IQR(x, na.rm = na.rm)
#   y <- x
#   y[x < (qnt[1] - H)] <- NA
#   y[x > (qnt[2] + H)] <- NA
#   y
# }


# Check variance in measurements using F-test
var.test(measurement_data[, 1], measurement_data[, 2], alternative = "two.sided")
var.test(measurement_data[, 3], measurement_data[, 2], alternative = "two.sided")
var.test(measurement_data[, 3], measurement_data[, 1], alternative = "two.sided")

# Check correlation
cor(measurement_data)

# Build some models
lm(diam_stem_mm ~ type., tomato_data)
lm(firmness_log_mm ~ type., tomato_data)
lm(leafarea_cm2 ~ type., tomato_data)

plot(measurement_data)

# Check impact of scaling
scaled_measurement <- as.data.frame(scale(measurement_data))
summary(scaled_measurement)

plot(scaled_measurement)
prcomp(scaled_measurement)


# Compare relationships ------------------------------

names(measurement_data)

measurement_data <- measurement_data %>%
  mutate(log(diam_stem_mm), log(leafarea_cm2))

stem_diam_model_full <- glm(diam_stem_mm ~ log(leafarea_cm2) * firmness_log_mm,
  measurement_data,
  family = gaussian
)

stem_diam_model_1 <- glm(diam_stem_mm ~ log(leafarea_cm2) + firmness_log_mm,
  measurement_data,
  family = gaussian
)

stem_diam_model_2 <- glm(diam_stem_mm ~ log(leafarea_cm2),
  measurement_data,
  family = gaussian
)

stem_diam_model_3 <- glm(diam_stem_mm ~ firmness_log_mm,
  measurement_data,
  family = gaussian
)

stem_diam_model_4 <- glm(diam_stem_mm ~ 1,
  measurement_data,
  family = gaussian
)

# Likelihood ratio test for differences between models
# Want model most similar to full balanced by lack of complexity
lrtest(stem_diam_model_full, stem_diam_model_1) # no significant difference, move to simpler model
lrtest(stem_diam_model_1, stem_diam_model_2) # significant difference, retain complex model
lrtest(stem_diam_model_1, stem_diam_model_3) # ditto
lrtest(stem_diam_model_1, stem_diam_model_4) # very much so

# Analyse the value of each variable using the Anova function from car
# note: Type II Anova
Anova(stem_diam_model_1, type = "II")
summary(stem_diam_model_1)

plot(stem_diam_model_1)


# Relationship between firmness and leaf area depends on type? -------

leafarea_model <- glm(leafarea_cm2 ~ firmness_log_mm * type.,
  tomato_data,
  family = gaussian
)

summary(leafarea_model)

leafarea_model_b <- glm(leafarea_cm2 ~ firmness_log_mm,
  b_tomato,
  family = gaussian
)

leafarea_model_c <- glm(leafarea_cm2 ~ firmness_log_mm,
  c_tomato,
  family = gaussian
)

leafarea_model_r <- glm(leafarea_cm2 ~ firmness_log_mm,
  r_tomato,
  family = gaussian
)

summary(leafarea_model_r)
summary(leafarea_model_c)
summary(leafarea_model_b)

ggplot(data = tomato_data, mapping = aes(x = firmness_log_mm, y = leafarea_cm2)) +
  geom_point(mapping = aes(colour = type.)) +
  geom_smooth(se = TRUE)

ggplot(data = tomato_data, mapping = aes(x = firmness_log_mm, y = leafarea_cm2, colour = type.)) +
  geom_point() +
  geom_smooth(se = FALSE)

# Compare some models -------------------------------------------------

# build the saturated model (different slopes and intercepts)
model_sat <- glm(leafarea_cm2 ~ firmness_log_mm * type.,
  family = gaussian,
  data = tomato_data
)

# build the model with different intercepts but common slope
model_int <- glm(leafarea_cm2 ~ firmness_log_mm + type.,
  family = gaussian,
  data = tomato_data
)

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
