# Core Ecology class population ecology modeling 
# October 2025

# The purpose of this script is to explore the population ecology models in
# Melbourne & Hastings (2008). This script first fits a basic deterministic 
# Ricker model by linear regression to data from a density experiment with 
# Tribolium castaneum. This initial simplified approach estimates the 
# parameters for the deterministic model, ignoring the clearly complex
# stochastic processes. Next, the script explores some of the stochastic 
# complexity in the paper by modeling demographic stochasticity and 
# environmental stochasticity. This script is mostly 12_basic_ricker_fit.R and
# 05_plotprodfunc.R by Brett Melbourne with edits and additions by Claire
# Winfrey.
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
# is Poisson Binomial). Look over Fig. 1 for a schematic of the model types and to
# understand the biology the models represent. 

# 1. Define (and think through!) custom functions from the paper that simulate
# data using these distributions (same functions used in the paper!). Talk with
# your small groups through the code for these functions (and ask for help if
# you need it). 

# It's okay if some of how the code works isn't totally clear. The goal is more
# to demystify these models by breaking them down into their component variables 
# and distributions, while thinking about how these align with the biology processes
# that they are modeling. First, look over the functions below and try to 
# answer the questions accompanying them among yourselves. Then, make the plots
# and re-discuss.

# i. This first custom function matches what we did above! Does it have
# any stochasticity? How do you know?
# Make sure you remember what Nt, R, and alpha represent (see above or the
# paper)
Ricker <- function(Nt, R, alpha){
  Nt * R * exp(-1 * alpha * Nt)
}

# ii. Here we are scaling up and incorporating some stochasticity. This function
# is used in the Poisson model. How is the stochasticity added in and what does it 
# model biologically (hint, look at the respective model in Fig. 1!). 
# Remember that Bernoulli in the paper is a type of Binomial distribution.
RickerStBS <- function(Nt, R, alpha){
  # Births and density independent survival (R = births*(1-mortality); mortality
  # is binomial, so compound distribution is Poisson with mean R).
  births <- rpois( length(Nt), Nt * R )
  # Density dependent survival
  survivors <- rbinom( length(Nt), births, exp( -1 * alpha * Nt ) )
  return(survivors)
}

# iii. Here is the first function specifically needed for the Poisson Binomial.
# How does it incorporate stochasticity? Also, pay attention to how the binomial
# simulation is then plugged into the Poisson simulation. Why does it make 
# sense biologically to model the binomial first?
RickerStSexBS <- function(Nt, R, alpha, p=0.5){
  females <- rbinom(length(Nt),Nt,p)
  # Births and density independent survival (R = births*(1-mortality) and
  # density dependent survival ( e^-alpha*Nt ); mortality is
  # binomial, so compound distribution is Poisson with mean Re^-aNt).
  survivors <- rpois( length(Nt), females * (1/p) * R * exp( -alpha * Nt ) )
  return(survivors)
}

# iv. This the second of two functions specifically needed for the Poisson + 
# binomial (again, see Fig. 1). 
# Does this function incorporate any stochasticity? (When you get to the 
# plotting stage, look at the code and think about why or why not)
# Notice how it returns the sum of the poisvar and sexvar. Can you figure
# out why it does this (again, may need to plot first!)
Ricker_poisbinom.var <- function(Nt,R,alpha) {
  poisvar <- Ricker(Nt,R,alpha)             #Poisson variance
  sexvar <- Ricker(Nt,R,alpha)^2 / Nt       #sex variance
  return( poisvar + sexvar )
}

# 2. Now we'll reproduce a few of the plots in Melbourne & Hasting (2008),
# Supplementary Figure 1.
# i. Parameters for simulating data
R <- 5
alpha <- 0.05
Nt <- c(2,6,12,seq(20,160,by=10)) 
reps <- 100

# ii. Set up for plots below:
xlim <- c(0,160) #x-axis limits
ylim <- c(0,100)
xtks <- seq(from = 0 , to = 160 , by = 50 ) #x ticks
ytks <- seq(from = 0 , to = 100 , by = 20 )
xpl <- 150 #x position of panel label
ypl <- 90
reps <- 100

# iii. Poisson model
# Create a new window to view the plots (i.e. not in the Plots panel to the right)
quartz() #CHANGE TO windows() if on a Windows OS computer!! 
par(mfrow = c(3, 1))
plot(1, 1,
     xlim = xlim, ylim = ylim, type = "n", axes = FALSE,
     xlab = expression(N[t]), ylab = expression(N[t+1]))
axis(1, at = xtks, labels = xtks)
axis(2, at = ytks, las = 1, labels = ytks)
box() #draw a box around plot
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

# iv. Poisson-binomial
plot(1, 1,
     xlim = xlim, ylim = ylim, type = "n", axes = FALSE,
     xlab = expression(N[t]), ylab = expression(N[t+1]))
axis(1, at = xtks, labels = xtks)
axis(2, at = ytks, las = 1, labels = ytks)
box() #draw a box around plot
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

# v. Questions after plotting
# 1) In both models, how was Ntp1 simulated stochastically?
# 2) Compare the two plots. What differences do you see? How do the 
# simulations change when more stochasticity is added in?
# 3) Why do you think the (non-stochastic) Ricker function was used
# to add the standard deviations to the plots?

# vi. Finally, we will plot simulate and plot the model that was the best fit
# for the red flour beetle data, the Negative binomial-binomial gamma. This
# incorporates two additional sources of stochasticity using an additional
# distribution, the Gamma, which is similar to an exponential distribution. 
# However, gammas are more difficult to simulate and a little more complex,
# so we will not go through them in depth as we did with Poisson and
# Bernoulli/binomial. The point of this section is, instead, to compare the 
# models that you hopefully understand well now with the final model, thinking
# about how adding more stochasticity reflects the biology.

# First, we'll define some additional parameters needed for the simulation
kD <- 0.5 #Shape parameter of gamma for birth rate
kE <- kD/alpha #k for env variation: Var same as NBdem-Ricker @ stat pt

# Next, we need more custom functions that were made for the paper
# 1. RickerStSexBS_DEhB. How many sources of stochasticity does this function bring in and what biological
# processes occurring for the beetles (or Ricker's original fish as in Fig.1)
# do they reflect?
RickerStSexBS_DEhB <- function(Nt, R, alpha, kD, kE, p=0.5) {
  # Heterogeneity in birth rate/DI mortality between times or locations
  Rtx <- rgamma(length(Nt),shape=kE,scale=R/kE)
  females <- rbinom(length(Nt),Nt,p)
  # Heterogeneity in individual birth rate plus density independent survival
  # (R = births*(1-mortality)*(1/p); mortality is binomial, so compound distribution
  # is negative binomial with mean (1/p)*R*Nt)
  births <- rnbinom( length(Nt), size = kD * females + (females==0),
                     mu = (1/p) * females * Rtx ) # Add 1 to size when females=0 to avoid NaN's.
  # Density dependent survival
  survivors <- rbinom( length(Nt), births, exp( -1 * alpha * Nt ) )
  return(survivors)
}

# 2. Function for adding in std later on. Note how it returns a combination
# of the 4 types of stochasticity/4 models!
Ricker_nbinombinomgamma.var <- function(Nt,R,alpha,kD,kE) {
  poisvar <- Ricker(Nt,R,alpha)             #Poisson variance
  sexvar <- Ricker(Nt,R,alpha)^2 / Nt       #sex variance
  dhvar <- Ricker(Nt,R,alpha)^2 / (kD*Nt)   #raw dhvar for no-sex model
  evar <- Ricker(Nt,R,alpha)^2 / kE         #env variance
  return( poisvar + sexvar + 2*dhvar + evar )
}

# Finally, plot this!
# Negative binomial-binomial gamma
# Note that this will probably plot the plot in the Plots panel to the right.
# Starting back on line 257 (quartz) and running the code quickly 
# through this plot should print all three models together in the separate
# window.
plot(1, 1,
     xlim = xlim, ylim = ylim, type = "n", axes = FALSE,
     xlab = expression(N[t]), ylab = expression(N[t+1]))
axis(1, at = xtks )
axis(2, at = ytks, labels = FALSE )
box()
for (i in 1:reps) {
  Ntp1 <- RickerStSexBS_DEhB(Nt, R, alpha, kD, kE)
  points(jitter(Nt),Ntp1,col="grey75")
}
lines(0:1100,Ricker(0:1100,R,alpha),col="grey50")
x <- Nt
y <- Ricker(x,R,alpha)
std <- sqrt(Ricker_nbinombinomgamma.var(x,R,alpha,kD,kE))
segments(x,y-std,x,y+std,col="black")
text(xpl,ypl,"NB-binomial-gamma",pos=2,offset=0)

### Some final questions to consider with all three plots.
# 1. What additional stochastic elements did the NBBg model add and, crucially,
# which biological processes do they reflect?
# 2. How do the plots differ from one another? How does adding more
# stochastic elements change the plots?
# 3. Consider your 3 plots alongside Fig. 3 in the paper. How do your plots
# relate to the results in Fig. 3? Why is stochasticity important to consider
# when thinking about extinction risk?
# 4. Finally, think about the 4 types of stochasticity considered in the paper, 
# which are rooted in the biology of red flour beetles (and that work for
# Ricker's original fish too!). 
## i. Would similar models be a good starting place for your own system? Why 
# or why not?
## ii. Can you imagine any systems where you would need to add additional types
# of stochasticity? What might be some good model types?


