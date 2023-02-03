# ---
# title: "Problem 3 of recommended exercise 5"
# author: "S. Muff"
# date: "`r format(Sys.time(), '%d %B, %Y')`"
# output:
#  html_document:
#    toc: yes
# ---

#+ include = FALSE
knitr::opts_chunk$set(warning = FALSE, message = FALSE, cache = TRUE)


#' ## Selection bias and the "wrong way to do CV".

 
  
#' ### Generate data

#' * Simulate high dimensional data ($p=5000$ predictors) from independent or correlated normal variables, but with few samples ($n=50$).

#' * Randomly assign class labels (here only 2). This means that the "truth"" is that the misclassification rate can not get very small. What is the expected misclassification rate (for this random set)?
#' 

#' ### Classification task:

#' * We choose a few ($d=25$) of the predictors (how? we just select those with the highest correlation to the outcome).
#' * Perform a classification rule (here: logistic empirical Bayes) on these predictors.
#' * Then we run CV ($k=5$) on either only the $d$ (=wrong way), or on all $c+d$ (=right way) predictors. 
#' * Report misclassification errors for both situations.

#' One possible version of this is presented in the R-code below. Go through the code and explain what is done in each step, then run the code and observe if the results are in agreement with what you expected. Make changes to the R-code if you want to test out different strategies.

#' 
#' We start by generating data for $n=50$ observations


library(boot)

#' GENERATE DATA; use a seed for reproducibility
set.seed(4268)
n=100 #number of observations
p=5000 #number of predictors
d=10 #top correlated predictors chosen
#' 
#' 
#' Generating predictor data
xs=matrix(rnorm(n*p,0,1),ncol=p,nrow=n) #simple way to to uncorrelated predictors
dim(xs) # n times p
xs[1:10,1:10]

#' Generate class labels independent of predictors - so if all classifies as class 1 we expect 50% errors in general
ys=c(rep(0,n/2),rep(1,n/2)) #now really 50% of each
table(ys)
 

#' ## *WRONG CV*:
#' 
#' Select the 25 most correlated predictors outside the CV.

#'
#' Calculate the correlations for all predictors and the response
corrs=apply(xs,2,cor,y=ys)
hist(corrs)
 

 
#' Select the d predictors with the highes correlation
selected=order(corrs^2,decreasing = TRUE)[1:d]  

#' Then store the highest correlated predictors into a reduced data frame, together with the response y:
data=data.frame(ys,xs[,selected])
 

#'Then run CV around the fitting of the classifier - use logistic regression and built in `cv.glm()` function
 
logfit=glm(ys~ . ,family="binomial",data=data)

cost <- function(r, pi = 0) mean(abs(r-pi) > 0.5)
kfold=10

cvres=cv.glm(data=data,glmfit=logfit,cost=cost,K=10)
cvres$delta
 

#'Observe a zero misclassification rate!
#'


#' ## *CORRECT CV*
#' 
#' Do not pre-select predictors outside the CV, but as part of the CV. We need to code this ourselves:

 
reorder=sample(1:n,replace=FALSE)
validclass=NULL
for (i in 1:kfold)
{
  neach=n/kfold
  trainids=setdiff(1:n,(((i-1)*neach+1):(i*neach)))
  traindata=data.frame(xs[reorder[trainids],],ys[reorder[trainids]])
  validdata=data.frame(xs[reorder[-trainids],],ys[reorder[-trainids]])
  colnames(traindata)=colnames(validdata)=c(paste("X",1:p),"y")
  foldcorrs= apply(traindata[,1:p],2,cor,y=traindata[,p+1]) 
  selected=order(foldcorrs^2,decreasing = TRUE)[1:d] #top d correlated selected
  data=traindata[,c(selected,p+1)]
  trainlogfit=glm(y~.,family="binomial",data=data)
  pred=plogis(predict.glm(trainlogfit,newdata=validdata[,selected]))
  validclass=c(validclass,ifelse(pred > 0.5, 1, 0))
}
table(ys[reorder],validclass)
1-sum(diag(table(ys[reorder],validclass)))/n
 


 