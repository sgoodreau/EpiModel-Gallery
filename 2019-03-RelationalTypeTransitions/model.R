
## Modeling an SI epidemic on a single network with two relational types marked by 
## an edge attribute
##
## Relational types differ by dissolution probability, coital frequency
##
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Author: Steven M. Goodreau (University of Washington)
## Date: March 2019

# Load EpiModel
suppressMessages(library(EpiModel))

# Standard Gallery unit test lines
rm(list = ls())
eval(parse(text = print(commandArgs(TRUE)[1])))


# Network model estimation ------------------------------------------------

## For now I am storing the list of cohabs as a separate network.  
## This is an inelegant way of doing this for sure - but it gets
##  around (1) needing to code up a new change statistics that 
## works on network's native edge attribute functionality 
## (edgecov does not, since it is a misnomer), and 
## (2) my untested assumption that tergmLite strips out edge covariates

n <- 1000
meandeg_non <- 0.3
meandeg_chb <- 0.5
diss_prob_non <- 0.01
diss_prob_chb <- 0.001
c2n_prob <- 0.03
concurrency <- 0.09
  
# Initialize the network
allrels <- network.initialize(n = n, directed = FALSE)
chbs <- network.initialize(n = n, directed = FALSE)

# Define and parametrize formation model
formation <- ~edges + concurrent

edges_non <- meandeg_non*n
edges_chb <- meandeg_chb*n
edges_tot <- edges_non + edges_chb
target.stats <- c(edges_tot, concurrency*n)

# Parameterize the dissolution model (same for Models 1 and 2)
coef.diss_main <- dissolution_coefs(
                    dissolution = ~offset(edges), duration = 50)
coef.diss

# Fit the models
est.mod1 <- netest(nw, formation.mod1, target.stats.mod1, coef.diss)

# Model diagnostics
dx.mod1 <- netdx(est.mod1, nsims = 10, ncores = 19, nsteps = 100,
                 set.control.ergm = control.simulate.ergm(MCMC.burnin = 1e5))
plot(dx.mod1, plots.joined = FALSE)






# Epidemic model simulation -----------------------------------------------

# Parameterizing an SIS epidemic
# Note: strain 1 is highly infectious but short-lived;
#       strain 2 has much lower infection but longer duration
param <- param.net(inf.prob = 0.5, inf.prob.st2 = 0.01,
                   rec.rate = 0.05, rec.rate.st2 = 0.005,
                   pct.st2 = 0.5,
                   act.rate = 2)

# Initial conditions
init <- init.net(i.num = 50)

# Read in the module functions
source("module-fx.R", echo = TRUE)

# Control settings
control <- control.net(nsims = 4,
                       ncores = 4,
                       nsteps = 1000,
                       infection.FUN = infection.2strains,
                       recovery.FUN = recov.2strains)

# Run the network model simulations with netsim
sim.mod1 <- netsim(est.mod1, param, init, control)
sim.mod2 <- netsim(est.mod2, param, init, control)


# Plotting results -----------------------------------------------

## In model 1, strain 1 dominates strain 2.
## In model 2, strain 1 goes extinct while strain 2 persists.
## The only difference between the two models was that one
## enforced a monogamy rule and the other did not.

par(mfrow = c(1, 2))
plot(sim.mod1, y = c("i.num.st1", "i.num.st2"),
     sim.lines = TRUE,
     sim.col = c("magenta", "lightgreen"),
     mean.line = TRUE, mean.lwd = 4,
     mean.col = c("purple", "green"),
     qnts = FALSE)

plot(sim.mod2, y = c("i.num.st1", "i.num.st2"),
     sim.lines = TRUE,
     sim.col = c("magenta", "lightgreen"),
     mean.line = TRUE, mean.lwd = 4,
     mean.col = c("purple", "green"),
     qnts = FALSE)


## Probing further  --------------------------------------------------------

# At what level of concurrency does the cross-over point occur?

if (interactive() == TRUE) { # won't run on Travis CI testing

# Check  how many nodes had concurrent ties on average in model 1
dx.mod1a <- netdx(est.mod1, nsims = 10, ncores = 10, nsteps = 100,
                  set.control.ergm = control.simulate.ergm(MCMC.burnin = 1e5),
                  nwstats.formula = ~edges+concurrent)
dx.mod1a             # Roughly 120

# Define model, and set up a vector of concurrency levels to explore
formation.mod3 <- ~edges + concurrent
concties.mod3 <- seq(0,120,10)
est.mod3 <- list()
sim.mod3 <- list()

# Run models
# Warning: this loop can take 30+ minutes to run
for (i in 1:length(concties.mod3)) {
  target.stats.mod3 <- c(300, concties.mod3[i])
  est.mod3[[i]] <- netest(nw, formation.mod3, target.stats.mod3, coef.diss)
  sim.mod3[[i]] <- netsim(est.mod3[[i]], param, init, control)
}

# Process output
i.num.st1.final.mod3 <- sapply(1:13, function(x) rowMeans(sim.mod3[[x]]$epi$i.num.st1)[1000])
i.num.st2.final.mod3 <- sapply(1:13, function(x) rowMeans(sim.mod3[[x]]$epi$i.num.st2)[1000])

# Plot results
par(mfrow = c(1, 1))
plot(concties.mod3, i.num.st1.final.mod3, type = "b", col = "purple",
     xlab = "exp. # of persons with concurrent ties",
     ylab = "prevalence of strains at time step 1000", ylim = c(0, 500))
points(concties.mod3, i.num.st2.final.mod3, type = "b", col = "green")
legend(0, 500, legend = c("strain 1", "strain 2"), col = c("purple", "green"), lwd = 2)

}
