# Core Ecology class population ecology modeling 

# The purpose of this script is to explore the population ecology models in
# Melbourne & Hastings (2008). This script first fits a basic deterministic 
# Ricker model by linear regression to data from a density experiment with 
# Tribolium castaneum. This initial simplified approach estimates the 
# parameters for the deterministic model, ignoring the clearly complex
# stochastic processes. Next, the script explores some of the stochastic 
# complexity in the paper by modeling demographic stochasticity and 
# environmental stochasticity. This script is mostly 12_basic_ricker_fit.R and
# 05_plotprodfunc.R by Brett Melbourne with some additions by Claire Winfrey.
# Brett's original code : https://github.com/melbourne-lab/stochastic-ricker)

############################################
# SET UP 
############################################
library(dplyr)
library(ggplot2)
library(lme4)

# Read in existing data from Melbourne & Hastings 2008
tribdata <- read.csv("data/ricker_data.csv")
View(read.delim("data/ricker_data_about.txt")) #familiarize yourself with the data
head(tribdata)

# Re-name columns to match paper
tribdata <- tribdata |>
    rename(Nt = At,
           Ntp1 = Atp1)
head(tribdata)

############################################
# I. FITTING BASIC DETERMINISTIC RICKER MODEL
############################################

# Basic base R plot showing relationship between population size at gen t and
# population size next generation
with(tribdata, plot(Nt, Ntp1))

# Basic ggplot 
tribdata |>
    ggplot(aes(x = Nt, y = Ntp1, col=factor(batch))) +
    geom_point()

# Calculate r (nb several -Inf due to extinctions). r is the growth rate for each of the patches from time t to time t+1
tribdata$r <- log(tribdata$Ntp1 / tribdata$Nt)
# Check out these values. 
tribdata$r  

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
?lmer #look at man page to understand code syntax below
fit <- lmer(r ~ Nt + (1|fbatch), data = filter(tribdata, r != -Inf))
fit

# Fit by ordinary linear regression (can remove batch given results above)
fit <- lm(r ~ Nt, data = filter(tribdata, r != -Inf))
fit
r_0 = coef(fit)[1]
alpha = -coef(fit)[2] #refer to paper for a reminder of what alpha is!
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
N[1] <- 20 #set intitial population size
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

############################################
# II. EXPLORING DISTRIBUTIONS
############################################
# Why do we use different distributions to model different biological processes?
# Refer to Fig. 2 in  Melbourne & Hastings (2008) for an overview of the models.
# We'll explore part of the stochastic complexity in Melbourne & Hastings (2008),
# specifically demographic stochasticity and stochasticity driven by the random
# process of sex determination. The former is modeled using a Poisson and the
# latter Bernoulli. To understand why, below is a refresher on these
# distributions, which are commonly used in ecological modeling and biology 
# in general!

# 1. Poisson-- Probability of a given number of events occurring in given 
# interval (usually, but not necessarily, time). Single parameter is lambda, 
# known mean rate of events per interval. 
# In the paper: Demographic stochasticity is modeled with a poisson, with lambda 
# in birth rates within individuals

# i. Simulate 100 beetle parents using the experimentally-determined (from our 
# model with the beetle data above) birthrate as lambda
R #look at number we derived
?rpois #man page for the function
set.seed(19) #ensure reproducibility by setting seed before (pseudo) random
# rpois simulation
rpois(n= 100, lambda = R)

# ii. Given the mean rate of R, what is the probability that a beetle parent 
# has 10 offspring over the course of its life? 3?
?dpois #man page for the function
dpois(x= 10, lambda = R)
# Only 3 offspring?
dpois(x= 3, lambda = R)
# iii. Plot example
# Define range of number of births
births <- 0:10 #what is probability of differing amounts of offspring (0-10)
# given lambda?
# Plot probability mass function (total area under the curve = 1)
plot(births, dpois(births, lambda=R), type='h')

# 2. Bernoulli -- Probability of achieving a success where 2 outcomes are 
# possible. Your typical 'coin toss'. Has one parameter, the probability of
# success on each trial.
# As the paper did, we will model stochastic sex determination with a Bernoulli
# distribution (assuming a 0.5 chance of being female as the 'success').
?rbinom #note that a Bernoulli trial is a special case of the broader binomial,
# with size = 1

# i. Simulate various number of trials (trials = beetle eggs), randomly assigning either
# 0 = male, 1 = female
set.seed(19) #ensure reproducibility
sexRatio1000 <- rbinom(n = 1000, size = 1, prob = 0.5) #1000 beetles
table(sexRatio1000) #how many are in each category?
set.seed(19) #ensure reproducibility
sexRatio100 <- rbinom(n = 100, size = 1, prob = 0.5) #100 beetles
table(sexRatio100) #how many are in each category?
set.seed(19) #ensure reproducibility
sexRatio10 <- rbinom(n = 10, size = 1, prob = 0.5) #10 beetles
table(sexRatio10) #how many are in each category?
# What happens to the evenness of the sex ratio with changes in population size?

# ii. An additional parameter that the paper models with a Bernoulli is whether
# or not a red flour beetle is cannibalized (although it's more complicated in
# the paper than this). 
# Experiment with different probabilities of cannibalization (assume that your
# prob is probability of not being eaten, with 1 = alive and 0 = dead), by 
# changing the code above
objectNameHere <- rbinom(n = , size = 1, prob = ) 
table(objectNameHere) #how many are in each category?
objectNameHere<- rbinom(n = , size = 1, prob = ) 
table(objectNameHere) #how many are in each category?

############################################
# III. MODEL DEMOGRAPHIC STOCHASTICITY AND SEX
############################################
# For this part of the script, we will delve deeper into two of the more simple
# stochastic models from the paper, specifically those modeling demographic 
# stochasticity on the level of individual beetle (Poisson model) and 
# demographic stochasticity plus sex determination (as Bernoulli, combined model
# is Poisson Binomial). See Fig. 1 for a schematic of the model types. 
# we will reproduce a few of the plots in 
# Melbourne & Hasting (2008), Supplementary Figure 1.

# Define custom functions that simulate data using these distributions 
# (same functions used in the paper!)
Ricker <- function(Nt, R, alpha){
  Nt * R * exp(-1 * alpha * Nt)
}

RickerStBS <- function(Nt, R, alpha){
  # Births and density independent survival (R = births*(1-mortality); mortality
  # is binomial, so compound distribution is Poisson with mean R).
  births <- rpois( length(Nt), Nt * R )
  # Density dependent survival
  survivors <- rbinom( length(Nt), births, exp( -1 * alpha * Nt ) )
  return(survivors)
}

Ricker_poisbinom.var <- function(Nt,R,alpha) {
  poisvar <- Ricker(Nt,R,alpha)             #Poisson variance
  sexvar <- Ricker(Nt,R,alpha)^2 / Nt       #sex variance
  return( poisvar + sexvar )
}

RickerStSexBS <- function(Nt, R, alpha, p=0.5){
  females <- rbinom(length(Nt),Nt,p)
  # Births and density independent survival (R = births*(1-mortality) and
  # density dependent survival ( e^-alpha*Nt ); mortality is
  # binomial, so compound distribution is Poisson with mean Re^-aNt).
  survivors <- rpois( length(Nt), females * (1/p) * R * exp( -alpha * Nt ) )
  return(survivors)
}

# Parameters for simulating data
R <- 5
alpha <- 0.05
kD <- 0.5 #Shape parameter of gamma for birth rate
kE <- kD/alpha #k for env variation: Var same as NBdem-Ricker @ stat pt
Nt <- c(2,6,12,seq(20,160,by=10)) 
reps <- 100

# Set up for plots below:
xlim <- c(0,160) #x-axis limits
ylim <- c(0,100)
xtks <- seq(from = 0 , to = 160 , by = 50 ) #x ticks
ytks <- seq(from = 0 , to = 100 , by = 20 )
xpl <- 150 #x position of panel label
ypl <- 90
reps <- 100

# Poisson model
plot(1,1,xlim=xlim,ylim=ylim,type="n",axes=FALSE)
axis(1, at = xtks, labels = FALSE )
axis(2, at = ytks, las=1 )
box()
for (i in 1:reps) {
  Ntp1 <- RickerStBS(Nt,R,alpha)
  points(jitter(Nt),Ntp1,col="grey75")
}
lines(0:1100,Ricker(0:1100,R,alpha),col="grey50")
x <- Nt
y <- Ricker(x,R,alpha)
std <- sqrt(y)
segments(x,y-std,x,y+std,col="black")
text(xpl,ypl,"Poisson",pos=2,offset=0)

# Poisson-binomial
plot(1,1,xlim=xlim,ylim=ylim,type="n",axes=FALSE)
axis(1, at = xtks )
axis(2, at = ytks, las=1 )
box()
for (i in 1:reps) {
  Ntp1 <- RickerStSexBS(Nt,R,alpha)
  points(jitter(Nt),Ntp1,col="grey75")
}
lines(0:1100,Ricker(0:1100,R,alpha),col="grey50")
x <- Nt
y <- Ricker(x,R,alpha)
std <- sqrt(Ricker_poisbinom.var(x,R,alpha))
segments(x,y-std,x,y+std,col="black")
text(xpl,ypl,"Poisson-binomial",pos=2,offset=0)
