#'---
#' title: Habitat selection of otters (an SSF analysis)
#' author: "S. Muff, J. Signer, J. Fieberg"
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   html_document: 
#'     toc: true
#'---
#'
#+ include = FALSE
knitr::opts_chunk$set(warning = FALSE, message = FALSE, cache = TRUE)
#'  
#' **Purpose**: This code replicates the analysis presented in Muff, Signer, Fieberg (2019) Section 4.2 "Habitat selection of otters: an SSF analysis".
#'

#' ## Load libraries and read in data
#+ warning=FALSE, message=FALSE
library(survival)
library(TwoStepCLogit)
library(INLA)
library(glmmTMB)
options(width=150)
dat <-  read.csv("d_otter.csv")
str(dat)

#' NAT1, REST1 and STAU1 are the three factor levels of the factor variable habitat type, encoded as dummy variables, where
#'
#' - NAT1: natural habitat (reference category)
#' - REST1: residual water
#' - STAU1: a reservoir
#'
#' Further, the two continuous variables in the model are:
#' 
#' - Sohlbrei: the river width
#' - Breaks_Dis: step length
#' 
#' Finally, `Loc` is the binary response variable that indicates if a habitat point was used (1) or available (0).
#'
#' ### Some data manipulation:

#' Add numerical variable for animals:
dat$ANIMAL_ID <- as.numeric(as.factor(dat$NA_ANIMAL))

#' Stratum ID is given as "NA_ID" in the data; 
#' It is easier to have sequential enumeration, so let's generate a new stratum-ID variable str_ID:
d.map <- data.frame(NA_ID=unique(dat$NA_ID),str_ID=1:length(unique(dat$NA_ID)))
dat$str_ID <- d.map[match(dat$NA_ID,d.map$NA_ID),"str_ID"]
dat <- dat[order(dat$str_ID),]

#' Scale and center the two continuous variables river width (Sohlenbrei) and step length (Breaks_Dis)
dat$Sohlenbrei <- scale(dat$Sohlenbrei)
dat$Breaks_Dis <- scale(dat$Breaks_Dis)


#' ## Fixed effects models 
 
#' ### clogit
r.clogit <- clogit(Loc ~ STAU1 + REST1 + Sohlenbrei + Breaks_Dis  
                   +   strata(str_ID), data=dat) 

summary(r.clogit)$coef


#' ### INLA

#' Set mean and precision for the priors of slope coefficients
mean.beta <- 0
prec.beta <- 1e-4  

#' The model formula for INLA, where we set the stratum-specific intercept variance to $10^6$ (or rather: the precision to $10^{-6}$), is given as follows:
formula.fixed <-  Loc ~  STAU1 + REST1 + Sohlenbrei +  Breaks_Dis +
  f(str_ID,model="iid",hyper=list(theta=list(initial=log(1e-6),fixed=T))) 

#' Then fit the Poisson model
r.inla.fixed <- inla(formula.fixed, family ="Poisson", data=dat,
                       control.fixed = list(
                         mean = mean.beta,
                         prec = list(default = prec.beta)
                       )
                )

#' The summary for the posterior distribution of the fixed effects:
r.inla.fixed$summary.fixed

#' ### glmmTMB
 
#' Note that we now have to manually fix the variance of the intercept first. 
#' 
#' Start by setting up the model, but do not yet fit it:
TMBStruc.fix = glmmTMB(Loc ~ STAU1 + REST1 + Sohlenbrei +
                         Breaks_Dis +  (1|str_ID), 
                       family=poisson, data=dat, doFit=FALSE) 

#' Then fix the standard deviation of the first random term, which is the `(1|str_ID)` component  in the above model equation:
TMBStruc.fix$parameters$theta[1] = log(1e3) 

#' We need to tell `glmmTMB` not to change the variance by setting it to `NA`:
TMBStruc.fix$mapArg = list(theta=factor(c(NA)))

#' Then fit the model and look at the results:
glmm.TMB.fixed = glmmTMB:::fitTMB(TMBStruc.fix) 
summary(glmm.TMB.fixed)
#'
#'
#' Note also that there is a one-step way to carry out the regression with newer versions of glmmTMB:

glmm.TMB.fixed = glmmTMB(Loc ~ STAU1 + REST1 + Sohlenbrei +
                           Breaks_Dis +  (1|str_ID), 
                         family=poisson, data=dat, 
                         map=list(theta=factor(c(NA))),
                         start=list(theta=c(log(1e3))))

#' ## Random effects models 


#' ### 2StepCLogit 

#' The two-step procedure with independent random effect (D="UN(1)"):
r.Twostep <-  Ts.estim(formula = Loc ~  STAU1 + REST1 + Sohlenbrei   + strata(NA_ID) + Breaks_Dis +
                     cluster(NA_ANIMAL) ,data = dat, random = ~ STAU1 + REST1 + Sohlenbrei ,
                   all.m.1=F, D="UN(1)") 

#' Slope estiamtes and standard errors
r.Twostep$beta
r.Twostep$se

#' Variance estimates
r.Twostep$D

#' ### INLA  
#'
#' To fit the model with random slopes in INLA, we need to generate new (but identical) variables of individual ID (ID cannot be used multiple times in the model formula):
dat$ANIMAL_ID1 <- dat$ANIMAL_ID
dat$ANIMAL_ID2 <- dat$ANIMAL_ID
dat$ANIMAL_ID3 <- dat$ANIMAL_ID

#' Set the model formula as for the fixed-effects model, but now add three random slope terms, namely for river width and for the two levels of the categorical variable (STAU1 and REST1) which are not the reference categories. The priors for precision of the three random slopes are PC(3,0.05), while the intercept variance is again fixed:
formula.random <- Loc ~  -1 + STAU1 + REST1 + Sohlenbrei +  
  Breaks_Dis +
  f(str_ID,model="iid",hyper=list(theta = list(initial=log(1e-6),fixed=T))) +
  f(ANIMAL_ID1,Sohlenbrei,values=1:9,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(ANIMAL_ID2,STAU1,values=1:9,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(ANIMAL_ID3,REST1,values=1:9,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) 

#' Fit the model
r.inla.random <- inla(formula.random, family ="Poisson", data=dat, 
                      control.fixed = list(
                        mean = mean.beta,
                        prec = list(default = prec.beta)
                      )
)

#' The summary for the posterior distribution of the fixed effects:
r.inla.random$summary.fixed 

#' Since variances are parameterized and treated as precisions, the summary of the respective posterior distributions is given for the precisions:
r.inla.random$summary.hyperpar
 

#' Source R functions for calculating posterior means 
#' and medians of the precisions.
source("inla_emarginal.R")
source("inla_mmarginal.R")
inla_emarginal(r.inla.random)
inla_mmarginal(r.inla.random)


#' ### glmmTMB

#' Now the same model using glmmTMB(). Again start to set up the model without fitting it:
TMBStruc = glmmTMB(Loc ~ -1 + STAU1 + REST1 + Sohlenbrei +  
                     Breaks_Dis +  (1|str_ID) + 
                     (0 + STAU1 | ANIMAL_ID) + 
                     (0 + REST1 | ANIMAL_ID)   + 
                     (0 + Sohlenbrei | ANIMAL_ID), 
                   family=poisson, data=dat, doFit=FALSE) 

#' Set the value of the standard deviation of the first random effect (here (1|str_ID)):
TMBStruc$parameters$theta[1] = log(1e3) 

#' Tell glmmTMB not to change the first standard deviation, all other values are freely estimated (and are different from each other)
TMBStruc$mapArg = list(theta=factor(c(NA,1:3)))

#' Fit the model and look at the summary:
glmm.TMB.random <- glmmTMB:::fitTMB(TMBStruc)
summary(glmm.TMB.random)

#' 95\% CIs for fixed and random effects (standard deviations) are obtained via the confint() function:
confint(glmm.TMB.random)

#' Note: It it has been a problem in a previous version of glmmTMB that the confint() function only showed a table of the length equal to the number of parameters estimated. As the variance for str_ID was not estimated but fixed to 10^6, but was still listed, the last variance component (here the one for (0 + Sohlenbrei | ANIMAL_ID)) was not shown. If you face this problem, you can solve the issue by moving the component (1|str_ID) to the last position in the formula, and then use 
TMBStruc$parameters$theta[4] = log(1e3) 
TMBStruc$mapArg = list(theta=factor(c(1:3, NA)))
#' in the above code.
#' 
#' 
#' Note also that there is a one-step way to carry out the regression with newer versions of glmmTMB:

glmm.TMB.random = glmmTMB(Loc ~ -1 + STAU1 + REST1 + Sohlenbrei +  
                            Breaks_Dis +  (1|str_ID) + 
                            (0 + STAU1 | ANIMAL_ID) + 
                            (0 + REST1 | ANIMAL_ID)   + 
                            (0 + Sohlenbrei | ANIMAL_ID),
                          family=poisson, data=dat,
                          map=list(theta=factor(c(NA,1:3))),
                          start=list(theta=c(log(1e3),0,0,0))
)
 
#' ## Illustration of convergence problem for 2StepCLogit
knitr::opts_chunk$set(error = TRUE)

#' Remove all strata of animal 6 where at least one point (used or available) lied in residual water (REST1) habitat; keep only strata that exclusively contain natural (NAT1) and reservoir (STAU1) points
#' 
#' strata to be removed:
strata.remove <- unique(dat[which(dat$ANIMAL_ID=="6" & dat$REST1==1),"NA_ID"])
#' Keep only strata that are not in this list; store in a reduced dataset:
dat.red <- dat[!(dat$NA_ID%in%strata.remove),]

#' Number of total data points (used or available) that are removed (from the 41670 in total):
nrow(dat) - nrow(dat.red)

#' Try to start the two-step procedure, but it does not run:
# Ts.estim(formula = Loc ~  STAU1 + REST1 + Sohlenbrei   + strata(NA_ID) + Breaks_Dis +
#   cluster(NA_ANIMAL) ,data = dat.red, random = ~ STAU1 + REST1 + Sohlenbrei ,
# all.m.1=F, D="UN(1)")

#' On the other hand, the same dataset is no problem for the Poisson model: 
TMBStruc = glmmTMB(Loc ~ -1 + STAU1 + REST1 + Sohlenbrei +  
                     Breaks_Dis +  (1|str_ID) + 
                     (0 + STAU1 | ANIMAL_ID) + 
                     (0 + REST1 | ANIMAL_ID)   + 
                     (0 + Sohlenbrei | ANIMAL_ID), 
                   family=poisson, data=dat.red, doFit=FALSE) 

TMBStruc$parameters$theta[1] = log(1e3) 
TMBStruc$mapArg = list(theta=factor(c(NA,1:3)))
glmm.TMB.random <- glmmTMB:::fitTMB(TMBStruc)

summary(glmm.TMB.random)


#' ## Session Info

devtools::session_info()