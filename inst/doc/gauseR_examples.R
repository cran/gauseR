## -----------------------------------------------------------------------------
# load package
require(gauseR)

# load data
data(gause_1934_book_f22)

logistic_data<-gause_1934_book_f22[gause_1934_book_f22$Treatment=="Pa",]

plot(Volume_Species2~Day, logistic_data)

## -----------------------------------------------------------------------------
# calculate time-lagged abundance
lagged_data<-get_lag(x=logistic_data$Volume_Species2, time = logistic_data$Day)

# calculate per-capita growth rate
lagged_data$dNNdt <- percap_growth(x = lagged_data$x, laggedx = lagged_data$laggedx, dt = lagged_data$dt)

# plot relationship
plot(dNNdt~laggedx, lagged_data,
     xlab="Lagged Abundance (N)",
     ylab="Per-capita Growth Rate (dN/Ndt)",
     xlim=c(0, 250), ylim=c(0, 1))
abline(h=0, v=0, lty=3)

# fit model to relationship
mod<-lm(dNNdt~laggedx, lagged_data)
abline(mod, lwd=2, col=2)

# label parameters
arrows(25, 0.6, 1, 0.8, length = 0.1, lwd=1.5)
text(25, 0.6, "y-intercept: r", pos=1)

arrows(200, 0.4, 232, 0.01, length = 0.1, lwd=1.5)
text(200, 0.4, "x-intercept: K", pos=3)

arrows(80, predict(mod, newdata=data.frame(laggedx=80)),
         130, predict(mod, newdata=data.frame(laggedx=80)),
         length = 0.1, lwd=1.5)
arrows(130, predict(mod, newdata=data.frame(laggedx=80)),
       130, predict(mod, newdata=data.frame(laggedx=130)),
       length = 0.1, lwd=1.5)
text(130, predict(mod, newdata=data.frame(laggedx=80)), "slope: aii", pos=3)

## -----------------------------------------------------------------------------
plot(Volume_Species2~Day, logistic_data, xlab="Day", ylab="P. aurelia Std. Volume")

timelst<-seq(0, 25, by=0.1) # sequence of time values to plot
prediction<-get_logistic(time = timelst, N0 = 0.6, r = 0.8, K=230)

lines(timelst, prediction, lwd=2)

## -----------------------------------------------------------------------------
prediction_short<-get_logistic(time = logistic_data$Day, N0 = 0.6, r = 0.8, K=230)

test_goodness_of_fit(observed = logistic_data$Volume_Species2, predicted = prediction_short)

## -----------------------------------------------------------------------------
# load data from competition experiment
data(gause_1934_book_f32)

# keep all data - no separate treatments exist for this experiment
predatorpreydata<-gause_1934_book_f32

# get time-lagged observations for each species
prey_lagged<-get_lag(x = predatorpreydata$Individuals_Prey, time = predatorpreydata$Day)
predator_lagged<-get_lag(x = predatorpreydata$Individuals_Predator, time = predatorpreydata$Day)

# calculate per-capita growth rates
prey_dNNdt<-percap_growth(x = prey_lagged$x, laggedx = prey_lagged$laggedx, dt = prey_lagged$dt)
predator_dNNdt<-percap_growth(x = predator_lagged$x,
      laggedx = predator_lagged$laggedx, dt = predator_lagged$dt)

# fit linear models to dNNdt, based on average
# abundances between current and lagged time steps
prey_mod_dat<-data.frame(prey_dNNdt=prey_dNNdt, prey=prey_lagged$laggedx,
      predator=predator_lagged$laggedx)
mod_prey<-lm(prey_dNNdt~prey+predator, data=prey_mod_dat)

predator_mod_dat<-data.frame(predator_dNNdt=predator_dNNdt,
      predator=predator_lagged$laggedx, prey=prey_lagged$laggedx)
mod_predator<-lm(predator_dNNdt~predator+prey, data=predator_mod_dat)

# model summaries
summary(mod_prey)
summary(mod_predator)

# extract parameters
# growth rates
r1 <- unname(coef(mod_prey)["(Intercept)"])
r2 <- unname(coef(mod_predator)["(Intercept)"])

# self-limitation
a11 <- unname(coef(mod_prey)["prey"])
a22 <- unname(coef(mod_predator)["predator"])

# effect of Pa on Pc
a12 <- unname(coef(mod_prey)["predator"])
# effect of Pc on Pa
a21 <- unname(coef(mod_predator)["prey"])

# run ODE:
# make parameter vector:
parms <- c(r1, r2, a11, a12, a21, a22)
initialN <- c(4, 0.1)
out <- deSolve::ode(y=initialN, times=seq(1, 17, length=100), func=lv_interaction, parms=parms)
matplot(out[,1], out[,-1], type="l",
        xlab="time", ylab="N", col=c("black","red"), lty=c(1,3), lwd=2, ylim=c(0, 60))
legend("topright", c("Pc", "Dn"), col=c(1,2), lwd=2, lty=c(1,3))

# now, plot in points from data
points(predatorpreydata$Day, predatorpreydata$Individuals_Predator , col=2)
points(predatorpreydata$Day, predatorpreydata$Individuals_Prey, col=1)

## -----------------------------------------------------------------------------
# Data for the optimizer:
# Must have a column with time steps labeled 'time', and
# columns for each species in the community.
opt_data<-data.frame(time=predatorpreydata$Day, Prey=predatorpreydata$Individuals_Prey,
    Predator=predatorpreydata$Individuals_Predator)

# Save the signs of the parameters -
# optimizer works in log space, so these
# must be specified separately
parm_signs<-sign(parms)

# parameter vector for optimizer -
# must be a vector with, first, the
# starting abundances in log space,
# and second, the parameter values,
# again in log space
pars<-c(log(initialN), log(abs(parms)))

# run optimizer
optout<-optim(par = pars, fn = lv_optim, hessian = TRUE,
             opt_data=opt_data, parm_signs=parm_signs)

# extract parameter vector:
parms <- exp(optout$par[-c(1:2)])*parm_signs
initialN <- exp(optout$par[1:2])

out <- deSolve::ode(y=initialN, times=seq(1, 17, length=100), func=lv_interaction, parms=parms)
matplot(out[,1], out[,-1], type="l",
        xlab="time", ylab="N", col=c("black","red"), lty=c(1,3), lwd=2, ylim=c(0, 60))
legend("topright", c("Pc", "Dn"), col=c(1,2), lwd=2, lty=c(1,3))

# now, plot in points from data
points(predatorpreydata$Day, predatorpreydata$Individuals_Predator , col=2)
points(predatorpreydata$Day, predatorpreydata$Individuals_Prey, col=1)

## -----------------------------------------------------------------------------
#load competition data
data("gause_1934_science_f02_03")

#subset out data from species grown in mixture
mixturedat<-gause_1934_science_f02_03[gause_1934_science_f02_03$Treatment=="Mixture",]

#extract time and species data
time<-mixturedat$Day
species<-data.frame(mixturedat$Volume_Species1, mixturedat$Volume_Species2)
colnames(species)<-c("P_caudatum", "P_aurelia")

#run wrapper
gause_out<-gause_wrapper(time=time, species=species)

## -----------------------------------------------------------------------------
#load competition data
data("gause_1934_science_f02_03")

#subset out data from species grown in mixture
mixturedat<-gause_1934_science_f02_03[gause_1934_science_f02_03$Treatment=="Mixture",]

#extract time and species data
time<-mixturedat$Day
species<-data.frame(mixturedat$Volume_Species1, mixturedat$Volume_Species2)
colnames(species)<-c("P_caudatum", "P_aurelia")

#run wrapper
gause_out<-gause_wrapper(time=time, species=species)

# number of species
N<-ncol(gause_out$rawdata)-1
# parameters
pars_full<-c(gause_out$parameter_intervals$mu)
# data.frame for optimization
fittigdata<-data.frame(y=unlist(gause_out$rawdata[,-1]),
                       time=gause_out$rawdata$time,
                       N=N)

yest<-ode_prediction(pars_full, time=fittigdata$time, N=fittigdata$N)
plot(fittigdata$y, yest, xlab="observation", ylab="prediction")
abline(a=0, b=1, lty=2)

#example of how to apply function, using nls()
mod<-nls(y~ode_prediction(pars_full, time, N),
           start = list(pars_full=pars_full),
           data=fittigdata)
summary(mod)

