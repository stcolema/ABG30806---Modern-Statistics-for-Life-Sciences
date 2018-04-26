#!/usr/bin/env Rscript

# Code for density plots and part 3 of the case study showing the 
# relaitonship between firmness and leaf area does depend upon type.

# Please remember to set the working directory to the location of TomatoesChris.xls

require(lmtest)
require(tidyverse)
require(gdata)

# Read in the xls file containing the data
datafile_name <- "TomatoesChris.xls"
tomato_data <- read.xls(datafile_name, sheet = 1, header = TRUE)

# Set the PlantNr to the row name and remove as a descriptive variale
rownames(tomato_data) <- tomato_data$PlantNr
tomato_data <- select(tomato_data, -PlantNr)

# Density plots -------------------------------------------------------

# Plot the densities of the measurements based on type
ggplot(tomato_data) +
  geom_density(mapping = aes(x = diam_stem_mm, colour = type.))

ggplot(tomato_data) +
  geom_density(mapping = aes(x = leafarea_cm2, colour = type.))

ggplot(tomato_data) +
  geom_density(mapping = aes(x = firmness_log_mm, colour = type.))

# Compare some models -------------------------------------------------

# build the saturated model (different slopes and intercepts)
model_sat <- glm(leafarea_cm2 ~ firmness_log_mm * type.,
                 family = gaussian,
                 data = tomato_data)

# build the model with different intercepts but common slope
model_int <- glm(leafarea_cm2 ~ firmness_log_mm + type.,
                 family = gaussian,
                 data = tomato_data)

# model with different intercepts and slope of 0
model_int_slope_0 <- glm(leafarea_cm2 ~ type.,
                         family = gaussian,
                         data = tomato_data)

# model with common intercept and slope
model_int_slope <- glm(leafarea_cm2 ~ type.,
                       family = gaussian,
                       data = tomato_data)

# Compare models using the likelihood ratio test
lrtest(model_sat, model_int) # insignificant difference
lrtest(model_int, model_int_slope_0) # significant difference
lrtest(model_int, model_int_slope) # significant difference

# Can claim that leaf area relationship with firmness does
# depend on the type (different intercepts for each type)
# but the the relative impact of a change in firmness upon 
# leaf area is common to all types (common slope)
