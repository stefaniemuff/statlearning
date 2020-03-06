library(latticeExtra)
library(INLA)

# load image including the mcmc output
path <- "/home/steffi/Subversion/measurementError/Framingham/"
savepath <- "/home/steffi/IFSPM/Conferences/BB2014/scientific/graphics/"
#path <- ""
load(paste(path,"inla_result_fram.RData",sep=""))

# read the MCMC samples
betaSamples <- data.frame(read.table(paste(path, "betaSamples_nP.txt", sep="")))
alphaSamples <- data.frame(read.table(paste(path, "alphaSamples_nP.txt", sep="")))
varUSamples <- as.matrix(read.table(paste(path, "varUSamples_nP.txt", sep="")))
varXSamples <- as.matrix(read.table(paste(path, "varXSamples_nP.txt", sep="")))

prior.beta = c(0, 0.01)
# naive estimates derived using INLA with comparable prior settings as in ME adjusted analysis
data<-read.table(file=paste(path,"frambook.txt",sep=""))
#data<-read.table("/home/steffi/IFSPM/Carroll_data/framingham/frambook.txt")
names(data)<-c("obs","age","sbp21","sbp22","sbp31","sbp32","smoke","chol2","chol3","Y")

data_indices<-which(data$age>=45 & data$chol3>=200 & data$chol3<=300)

data<-data[data_indices,]

n<-nrow(data)

z<-data$smoke
W1<-log(data$sbp31-50)
W2<-log(data$sbp32-50)
T1<-log(data$sbp21-50)
T2<-log(data$sbp22-50)

y<-data$Y

dd<-data.frame(cbind(z,W1,W2,T1,T2,y))

w1<-apply(dd[,2:3],1,mean)
w2<-apply(dd[,4:5],1,mean)
## using w1 and w2 explicitly
# w1 <- w1 -mean(w1)
# w2 <- w2 -mean(w2)
# 
# Ntrials = rep(1,2*n)
# w <- c(w1, w2)
# Z = c(z, z)
# 
# Y = matrix(NA, 2*n, 1)
# Y[1:n, 1] = y
# Y[n+(1:n), 1] = y
# 
# data_fram <- list(
#   Y=Y,
#   w=w,
#   Z=Z,
#   Ntrials=Ntrials)
# 
# formula = Y ~  w + Z
# 
# r_naive = inla(formula, Ntrials=Ntrials, data = data_fram,
#          family = "binomial",
#          control.fixed = list(
#            mean.intercept = prior.beta[1], 
#            prec.intercept = prior.beta[2],
#            mean = prior.beta[1],
#            prec = prior.beta[2]),
#          verbose=T)
# 
# coef_naive_b0 <- r_naive$summary.fixed["(Intercept)", c(1,3,5)]         
# coef_naive_bw <- r_naive$summary.fixed["w", c(1,3,5)]         
# coef_naive_bz <- r_naive$summary.fixed["Z", c(1,3,5)]         

## using the mean of w1 and w2
w <- 0.5*(w1+w2) 
w <- w - mean(w)


Ntrials = rep(1,n)

data_fram <- list(
  y=y,
  w=w,
  z=z,
  Ntrials=Ntrials)

formula = y ~  w + z

r_naive = inla(formula, Ntrials=Ntrials, data = data_fram,
         family = "binomial",
         control.fixed = list(
           mean.intercept = prior.beta[1], 
           prec.intercept = prior.beta[2],
           mean = prior.beta[1],
           prec = prior.beta[2]),
         verbose=F, control.inla=list(strategy="laplace"))

coef_naive_b0 <- r_naive$summary.fixed["(Intercept)", c(1,3,5)]         
coef_naive_bw <- r_naive$summary.fixed["w", c(1,3,5)]         
coef_naive_bz <- r_naive$summary.fixed["z", c(1,3,5)]       
    
    
         
# Ntrials = rep(1,n)
# regr<-glm(cbind(y,Ntrials-y)~w+z,family="binomial")
# coef_naive<-summary(regr)$coef[,1]
# CI_naive_b0<- confint(regr)[1,] 
# CI_naive_bx<- confint(regr)[2,] 
# CI_naive_bz<- confint(regr)[3,] 

# INLA
coef_inla_b0<- r$summary.fixed["mu",c(1,3,5)]
coef_inla_bx<-r$summary.hyperpar["Beta for idx.xx",c(1,3,5)]
coef_inla_bz<- r$summary.fixed["Z",c(1,3,5)]

# MCMC
coef_mcmc_b0<-mean(betaSamples[,1])
coef_mcmc_bx<-mean(betaSamples[,2])
coef_mcmc_bz<-mean(betaSamples[,3])

save(coef_naive_bw, coef_inla_bx, file=paste(path,"Fram_estimates.Rdata",sep=""))

CI_mcmc_b0<-quantile(betaSamples[,1], probs=c(0.025, 0.975))
CI_mcmc_bx<-quantile(betaSamples[,2], probs=c(0.025, 0.975))
CI_mcmc_bz<-quantile(betaSamples[,3], probs=c(0.025, 0.975))

# Bootstrap
beta<-data.frame(read.table(paste(path, "bootstrap_results.txt", sep="")))
betax<-beta[,1]
betaz<-beta[,2]

mybreaks <- function(data){
  nbins <- nclass.scott(data)
  h <- diff(pretty(data, nbins))[1L]
  x0 = -h/1000
  
  first <- floor((min(data) - x0)/h)
  last <- ceiling((max(data) - x0)/h)
  breaks <- x0 + h * c(first:last)
  return(breaks)
}


bx<-cbind(c(coef_naive_bw[[1]],1.76,coef_inla_bx[[1]]),
c(coef_naive_bw[[2]],0.70,coef_inla_bx[[2]]),
c(coef_naive_bw[[3]],2.82,coef_inla_bx[[3]]))
bz<-cbind(c(coef_naive_bz[[1]],0.38,coef_inla_bz[[1]]),
c(coef_naive_bz[[2]],-0.23,coef_inla_bz[[2]]),
c(coef_naive_bz[[3]],0.99,coef_inla_bz[[3]]))



data <- data.frame(rbind(bx, bz))
data[,"param"] <- rep(c("betax", "betaz"), each=3)
data[,"method"] <- rep(c("NAIVE", "ML", "INLA"), 2)
names(data) <- c("est", "low", "up", "param", "method")
data$method <- factor(data$method, levels=c( "INLA",  "ML", "NAIVE"))
data$param <- factor(data$param, levels=c("betax", "betaz"))

xlim2 <- c(0.6,3.2) 
xlim3 <- c(-0.3,1.1) 


data <- data
save(data, file=paste(path,"Fram_CI.Rdata",sep=""))

library(latticeExtra)
pdf(paste(savepath,"pe_fram.pdf",sep=""), width=11, height=3.5)
segplot(method~low+up|param, data=data,
        scales=list(x=list(limits=list(xlim2,xlim3),
            relation="free", cex=1.2),y=list(draw=TRUE,# labels=c("SIMEX", "INLA", "MCMC", "Naive"), 
            cex=1.2)),
        col=1, col.symbol=1,
        ylab="", xlab="",
        draw.bands = FALSE, centers = est,
        horizontal=T,
        strip = strip.custom(factor.levels = c(expression(beta[x]),
                                               expression(beta[z])),
                             bg="grey90"),
        segments.fun = panel.arrows, ends = "both", 
        angle = 90, length = 1, unit = "mm",
        panel=function(x,y, z, subscripts, centers, ...){
#        print(data[subscripts,])
        panel.grid(h=-1,v=-1)
        panel.segplot(x,y, z,subscripts,
                        centers=centers,
                        ...)
        dataSub <- data[subscripts,]
        panel.abline(v=subset(dataSub, method=="NAIVE")$est, col=1,lty=2)
        },
        par.strip.text=list(cex=1.5, lines=1.7)
)
dev.off()
