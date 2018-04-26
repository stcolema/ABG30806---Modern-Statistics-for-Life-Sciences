#!/usr/bin/env Rscript

rm(list = ls())
library(car)
library(ggplot2)
library(reshape2)
library(emmeans)
library(car)

# Set to working directory
my_wd <- paste0(
  "C:/Users/steph/Desktop/Bioinformatics/ABG30806",
  " - Modern Statistics for Life Sciences/Case 1 ",
  "(Chris)"
)

setwd(my_wd)

## Read in datafile ##
# df_tomato <- read.csv("TomatoesChris.csv", sep = ";")

# Read in the xls file containing the data (Stephen version)
require(gdata)
datafile_name <- "TomatoesChris.xls"
df_tomato <- read.xls(datafile_name, sheet = 1, header = TRUE)

# Function ----------------------------------------------------------

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

# Data exploration ------------------------------------------------------


## Convert the type column ##
df_tomato$type. <- as.factor(df_tomato$type.)
rownames(df_tomato) <- df_tomato$PlantNr
df_tomato <- select(df_tomato, -PlantNr)

# Check if the design is balanced
table(df_tomato$type.)
summary(df_tomato)

# Comment this please
sem <- function(y) {
  sd(y) / sqrt(length(y))
}

apply(df_tomato[, 2:4], 2, sd)
apply(df_tomato[, 2:4], 2, var)
apply(df_tomato[, 2:4], 2, sem)

dim(df_tomato)

## hist of quantive variable ##
par(mfrow = c(2, 2))
hist(df_tomato$diam_stem_mm, xlab = "Stem diameter in mm", main = "", breaks = 12)
hist(df_tomato$leafarea_cm2, xlab = "Leaf area in cm2", main = "", breaks = 12)
hist(df_tomato$firmness_log_mm, xlab = "Log(firmness in mm)", main = "", breaks = 12)
plot(df_tomato$type.)
par(mfrow = c(1, 1))
#### Are there differences between the three measurements and types of tomatoes? ####

## qqPlot to check normal distribution
qqPlot(df_tomato$diam_stem_mm)
qqPlot(df_tomato$leafarea_cm2)
qqPlot(df_tomato$firmness_log_mm)

## Check if there are any relationships
plot(df_tomato)

## All boxplots in one ##
dat_m <- melt(df_tomato,
  id.vars = "type.",
  measure.vars = c("diam_stem_mm", "firmness_log_mm", "leafarea_cm2")
)

p <- ggplot(dat_m) +
  geom_boxplot(aes(x = type., y = value, color = variable))

p1 <- ggplot(df_tomato) +
  geom_boxplot(mapping = aes(x = type., y = diam_stem_mm))

p2 <- ggplot(df_tomato) +
  geom_boxplot(mapping = aes(x = type., y = leafarea_cm2))

p3 <- ggplot(df_tomato) +
  geom_boxplot(mapping = aes(x = type., y = firmness_log_mm))

multiplot(p1, p2, p3, cols = 2) # Note how type c tends to have less / no overlap with b and r

## 3times boxplot ##
par(mfrow = c(2, 2))
plot(df_tomato$type., df_tomato$leafarea_cm2, xlab = "Type", ylab = "Leafarea cm2")
plot(df_tomato$type., df_tomato$diam_stem_mm, xlab = "Type", ylab = "Diameter stem mm")
plot(df_tomato$type., df_tomato$firmness_log_mm, xlab = "Type", ylab = "Log(firmness in mm)")
par(mfrow = c(1, 1))

## We showed here that there are differences ##
# So is this the correct way of concluding significance?? #
par(mfrow = c(2, 2))
mod <- lm(diam_stem_mm ~ type., data = df_tomato)
anova(mod)
summary(mod)
plot(mod)

mod <- lm(diam_stem_mm ~ type., data = df_tomato)
anova(mod)
summary(mod)
plot(mod)

mod2 <- lm(firmness_log_mm ~ type., data = df_tomato)
anova(mod2)
summary(mod2)
plot(mod2)

mod3 <- lm(leafarea_cm2 ~ type., data = df_tomato)
anova(mod3)
summary(mod3)
plot(mod3)

par(mfrow = c(1, 1))

#### Are the measurements related to each other? ####
## Measure type difference ##

# Model with all terms
df_tomato_anova_wo_int <- lm(diam_stem_mm ~ firmness_log_mm + leafarea_cm2 + type., data = df_tomato)
df_tomato_anova_w_int <- lm(diam_stem_mm ~ firmness_log_mm + leafarea_cm2 + type. + leafarea_cm2 * firmness_log_mm, data = df_tomato)
df_tomato_anova_w_int_2 <- lm(firmness_log_mm ~ diam_stem_mm + leafarea_cm2 + type., data = df_tomato)
df_tomato_anova_w_int_3 <- lm(firmness_log_mm ~ type. + diam_stem_mm + leafarea_cm2 + diam_stem_mm * leafarea_cm2, data = df_tomato)
df_tomato_anova_w_int_2_wo_type <- lm(firmness_log_mm ~ diam_stem_mm + leafarea_cm2, data = df_tomato)

anova(df_tomato_anova_w_int_3, df_tomato_anova_w_int_2)

plot(df_tomato_anova_w_int_2)

Anova(df_tomato_anova_w_int_2_wo_type, type = "II")
plot(df_tomato_anova_w_int_2_wo_type)
summary(df_tomato_anova_w_int_2)
emm <- emmeans(df_tomato_anova_w_int_2, ~ type.)
cld(emm, alpha = 0.05, adjust = "tukey")

### ekf ###
library(leaps)
x <- model.matrix(df_tomato_anova_w_int_3)[, -1]
y <- df_tomato$diam_stem_mm
adjr <- leaps(x, y, method = "adjr2")
library(faraway)
maxadjr(adjr, 8)

anova(df_tomato_anova_wo_int, df_tomato_anova_w_int)
Anova(df_tomato_anova_wo_int, type = "II")

# Firmness * Area
anov_tom_typ_no <- lm(diam_stem_mm ~ type. + leafarea_cm2 + firmness_log_mm, data = df_tomato)
emm <- emmeans(anov_tom_typ_no, ~ type.)
cld(emm, alpha = 0.05, adjust = "tukey")
