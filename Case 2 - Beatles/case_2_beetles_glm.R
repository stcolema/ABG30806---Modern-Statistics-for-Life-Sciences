#!/usr/bin/env Rscript

library(tidyverse)
library(magrittr) # install.packages("magrittr", dep = T)
library(ggfortify) # install.packages("ggfortify", dep = T)
library(reshape2) # install.packages("reshape2", dep = T)
require(lmtest) # install.packages("lmtest", dep = T)

# === Functions ================================================================

# Multiple plot function from
# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
# Allows multiple ggplot2 plots on the same image
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

# Create a plot of predicted line vs actual points including error bars
plot_fitted_errors <- function(model, data, x_var, y_var,
                               title = paste(y_var, "vs", x_var),
                               x_name = x_var,
                               y_name = y_var) {
  predictions <- fitted(model)
  p <- ggplot(data = data, mapping = aes_string(x = x_var, y = y_var)) +
    geom_point() +
    geom_line(y = predictions, size = 1, colour = "blue") +
    geom_segment(aes(
      xend = data[[x_var]],
      yend = predictions,
      color = "error"
    )) +
    labs(
      title = title,
      x = x_name,
      y = y_name,
      color = "series"
    )
  return(p)
}

# === Main code ================================================================

# === Data input ===============================================================

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
insect_dose %<>%
  mutate(prop_dead = num_died / num_insects)

# === Linear regression ========================================================
# Build model
prop_lm <- lm(prop_dead ~ dose, data = insect_dose)

# Inspect model
# Ugh, does not look good
par(mfrow = c(2, 2))
plot(prop_lm)
par(mfrow = c(1, 1))

ggplot(data = insect_dose, mapping = aes(x = dose, y = prop_dead)) +
  geom_point() +
  labs(
    title = "Proportion dead vs dose",
    x = "Dose",
    y = "Proportion dead"
  )

# Plot fiited values versus actual
p_lm <- plot_fitted_errors(prop_lm, insect_dose, "dose", "prop_dead",
  title = "LM regression",
  x_name = "Dose",
  y_name = "Proportion of sample dead"
)

p_lm

# Compare to Local Polynomial Regression
# (overfitting, but should give better idea of ideal curve)
loess_mod <- loess(prop_dead ~ dose, data = insect_dose)

p_loess <- plot_fitted_errors(loess_mod, insect_dose, "dose", "prop_dead",
  title = "LOESS regression",
  x_name = "Dose",
  y_name = "Proportion of sample dead"
)

# Compare LM fit with LOESS fit in single graph
ggplot(data = insect_dose, mapping = aes(x = dose, y = prop_dead)) +
  geom_point() +
  geom_line(y = fitted(prop_lm), size = 1, colour = "blue") +
  geom_line(y = fitted(loess_mod), size = 1, colour = "red") +
  geom_abline(intercept = 0.59, slope = 0) +
  geom_vline(xintercept = 1.7912) +
  labs(
    title = "LM regression vs ideal curve",
    x = "Dose",
    y = "Proportion dead"
  )

# Inspect fitted values
predictions <- fitted(prop_lm)
residuals <- rstudent(prop_lm) # studentised residuals
lm_comp_df <- data.frame(
  actual = insect_dose$prop_dead,
  fitted = predictions,
  studentised_residuals = residuals,
  errors = insect_dose$prop_dead - predictions,
  index = 1:8
)

lm_comp_df$errors / lm_comp_df$actual

# Plot residuals
ggplot(data = lm_comp_df, aes(x = index, y = studentised_residuals)) +
  geom_point(colour = "blue") +
  geom_abline(intercept = 0, slope = 0, colour = "red") +
  geom_segment(xend = lm_comp_df$index, yend = 0) +
  ylim(-3, 3) +
  labs(
    title = "Studentised residuals for LM regression",
    x = "Index",
    y = "Studentised residuals"
  )

# Diagnostic plots
autoplot(prop_lm, which = 1:6, label.size = 3)

# Compare plots; LOESS reveals sigmoidal relationship
# Logistic regression or Cloglog or Probit then
multiplot(p_lm, p_loess, cols = 1)

# === Logistic Regression ======================================================

# Build a logistic regression model
prop_log_reg <- glm(prop_dead ~ dose, family = binomial, data = insect_dose)
summary(prop_log_reg)


p_log_reg <- plot_fitted_errors(prop_log_reg, insect_dose, "dose", "prop_dead",
  title = "Logistic regression",
  x_name = "Dose",
  y_name = "Proportion of sample dead"
)

# Compare plot with LOESS regression and LM
multiplot(p_log_reg, p_loess) # Note heavier tails in LOESS
multiplot(p_log_reg, p_lm)

par(mfrow = c(3, 2))
plot(prop_log_reg, which = 1:6)
par(mfrow = c(1, 1))
plot(prop_log_reg, which = 1)

# Inspect residuals and deviance
prop_log_reg$residuals
prop_log_reg$deviance

# Put fitted values in table with actual data
comparison_table <- data.frame(
  Actual = insect_dose$prop_dead,
  L.R.predicted = fitted(prop_log_reg)
)

# === Cloglog Regression =======================================================

# Try a different link function (with heavier tails)
cloglog_model <- glm(prop_dead ~ dose,
  family = binomial(link = "cloglog"),
  data = insect_dose
)

p_cloglog <- plot_fitted_errors(cloglog_model, insect_dose, "dose", "prop_dead",
  title = "Cloglog regression",
  x_name = "Dose",
  y_name = "Proportion of sample dead"
)

# compare with logistic regression
multiplot(p_log_reg, p_cloglog)

# === Comparison binomial models ===============================================

# Deviance comparison
anova(prop_log_reg, cloglog_model, test = "Chisq")
summary(prop_log_reg)
summary(cloglog_model)
# Cloglog has less residual deviance and they have common df

# add to table of predictions
comparison_table$Cloglog.predicted <- fitted(cloglog_model)

# add error variables for each model
comparison_table %<>%
  mutate(
    L.R.errors = Actual - L.R.predicted,
    Cloglog.errors = Actual - Cloglog.predicted,
    Index = rownames(.)
  ) %>%
  as.tibble()

# Compare summary of Logistic regression errors with Cloglog errors
summary(comparison_table)

# Compare errors
ggplot(data = comparison_table) +
  geom_density(aes(x = L.R.errors), fill = "grey") +
  geom_density(aes(x = Cloglog.errors), fill = "red")

# Convert the comparison table to long data for plotting purposes
smaller_table <- comparison_table %>%
  dplyr::select(Actual, L.R.errors, Cloglog.errors) %>%
  melt(id = c("Actual")) %>%
  rename(Model = variable, Error = value) %>%
  mutate(Model = as.character(Model)) %>%
  as.tibble() %>%
  mutate(
    Model = replace(Model, Model == "L.R.errors", "Logistic"),
    Model = replace(Model, Model == "Cloglog.errors", "Cloglog")
  )

summary(smaller_table)

# Plot the densities of the errors for the binomial models
ggplot(data = smaller_table) +
  geom_density(aes(x = Error, colour = Model)) +
  geom_vline(xintercept = 0, lty = 2) +
  xlim(-0.15, 0.15) +
  labs(
    title = "Distribution of errors for binomial models",
    x = "Error",
    y = "Density"
  )

# === Multivariate model =======================================================

insect_dose %<>% mutate(dose_quad = dose * dose)

# Build multivariate model
prop_log_reg_multi <- glm(prop_dead ~ dose + dose_quad,
  data = insect_dose,
  family = binomial(link = "logit")
)

# Likelihood ratio test
lrtest(prop_log_reg, prop_log_reg_multi) # p = 0.51 => retain simpler model
anova(prop_log_reg, prop_log_reg_multi, test = "Chisq")

# Compare non-nested models (use AIC)
AIC(cloglog_model, prop_log_reg, prop_log_reg_multi) # Cloglog wins

# Update the comparison table
multi_comp_table <- data.frame(
  Actual = insect_dose$prop_dead,
  Error = insect_dose$prop_dead - fitted(prop_log_reg_multi),
  Model = "Logistic multivariate"
)

smaller_table %<>% bind_rows(multi_comp_table)

# Plot the densities of the errors for the binomial models
ggplot(data = smaller_table) +
  geom_density(aes(x = Error, colour = Model)) +
  geom_vline(xintercept = 0, lty = 2) +
  xlim(-0.15, 0.15) +
  labs(
    title = "Distribution of errors for binomial models",
    x = "Error",
    y = "Density"
  )
