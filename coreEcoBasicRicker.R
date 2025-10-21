# This script fits a basic deterministic Ricker model by linear regression to
# data from a density experiment with Tribolium castaneum.

# This approach is suitable for an introductory ecology computer lab or a lab
# conducting a live Tribolium experiment. The simplified approach estimates the
# parameters for the deterministic model, ignoring the clearly complex
# stochastic processes.

library(dplyr)
library(ggplot2)
library(lme4)

# Read in existing data from Melbourne & Hastings 2008 or replace with your own

tribdata <- read.csv("data/ricker_data.csv")
tribdata <- tribdata |>
    rename(Nt = At,
           Ntp1 = Atp1)
head(tribdata)

# Basic base R plot

with(tribdata, plot(Nt, Ntp1))

# Basic ggplot

tribdata |>
    ggplot(aes(x = Nt, y = Ntp1, col=factor(batch))) +
    geom_point()

# Calculate r (nb several -Inf due to extinctions)

tribdata$r <- log(tribdata$Ntp1 / tribdata$Nt)


# Plot r vs Nt

# We see it is nice and linear and batch has little effect. We also see the much
# greater variance at small Nt and a bunch of -Infs where small initial
# populations went extinct. These issues are dealt with in the more complex
# stochastic models.

tribdata |>
    ggplot(aes(x = Nt, y = r, col=factor(batch))) +
    geom_point()

# Fit r_0 and alpha by linear mixed model to account for batch, excluding
# extinctions. We see the batch variance is estimated to be 0.

tribdata$fbatch <- factor(tribdata$batch)
fit <- lmer(r ~ Nt + (1|fbatch), data = filter(tribdata, r != -Inf))
fit

# Fit by ordinary linear regression

fit <- lm(r ~ Nt, data = filter(tribdata, r != -Inf))
fit
r_0 = coef(fit)[1]
alpha = -coef(fit)[2]
tribdata |>
    ggplot(aes(x = Nt, y = r)) +
    geom_point() + 
    geom_abline(slope = -alpha, intercept = r_0) +
    labs(title = "Fitted model")

# Predict the deterministic dynamics
# We see the population reaches carrying capacity at about
# generation 7 if started with N_0 = 20.

R <- exp(r_0)
t <- 0:15
N <- t * NA
N[1] <- 20
for (i in 1:max(t)) {
   N[i+1] <- R * N[i] * exp(-alpha * N[i])
}
rickersim <- data.frame(t, N)
rickersim |>
    ggplot(aes(x=t, y=N)) +
    geom_line() +
    geom_point() +
    ylim(0, 250) +
    labs(title = "Predicted dynamics")


