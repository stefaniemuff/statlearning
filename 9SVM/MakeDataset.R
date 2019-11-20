setwd("/Users/tjroksva/Documents/H2017/StatL/Module9/")

set.seed(9)
library(MASS)
linje=function(x){return(1.3-1.1*x)}
giklasse=function(coords){
  klasse=c()
  for(i in 1:dim(coords)[1]){
    y=coords[i,2];x=coords[i,1]
    if(y-1.2+1.1*x+0.08<0){
      klasse[i]=-1
    }
    if(y-1.2+1.1*x-0.08>0){
      klasse[i]=1
    }
  }
  return(klasse)
}


giklasse2=function(coords){
  klasse=c()
  for(i in 1:dim(coords)[1]){
    
    y=coords[i,2];x=coords[i,1]
    
    if((y-1.3+1.1*x<0) & (y-1+1.1*x>0)){
    
        klasse[i]=sample(c(-1,1),size=1)
  
          }else if(y-1.3+1.1*x<0){
      klasse[i]=-1
      
    }else{
      klasse[i]=1
    }
  }
  return(klasse)
}

n=100
coords=cbind(runif(n,0,1),runif(n,0,1))
klasse=giklasse(coords)
par(mfrow=c(1,3));par(pty="s")
plot(NA,xlab="x",ylab="y",xlim=c(0,1),ylim=c(0,1));title("Forest 1")
points(coords[klasse==-1,],col="lightcoral",pch=19,cex=0.9)
points(coords[klasse==1,],col="darkseagreen",pch=19,cex=0.9)

ind_na=which(is.na(klasse)==TRUE)
forest1=cbind(coords[-ind_na,],klasse[-ind_na])

write.table(forest1,file="forest1.txt",col.names=FALSE,row.names=FALSE)

#lines(seq(0,1,length.out=100),linje(seq(0,1,length.out=100)),lwd=2)
plot(NA,xlab="x",ylab="y",xlim=c(0,1),ylim=c(0,1));title("Forest 2")
coords=cbind(runif(n,0,1),runif(n,0,1))
klasse=giklasse2(coords)
if(length(klasse)!=dim(coords)[1]){klasse[length(klasse)+1]=NA}

points(coords[klasse==-1,],col="lightcoral",pch=19,cex=0.9)
points(coords[klasse==1,],col="darkseagreen",pch=19,cex=0.9)

forest2=cbind(coords,klasse)

write.table(forest2,file="forest2.txt",col.names=FALSE,row.names=FALSE)

plot(NA,xlab="x",ylab="y",xlim=c(0,1),ylim=c(0,1));title("Forest 3")
clust1=mvrnorm(40,mu=c(0.2,0.1),Sigma=diag(c(0.05,0.07)))
points(clust1,col="lightcoral",pch=19,cex=0.9)
clust2=mvrnorm(20,mu=c(0.5,0.5),Sigma=diag(c(0.02,0.02)))
points(clust2,col="darkseagreen",pch=19,cex=0.9)
clust3=mvrnorm(40,mu=c(0.9,0.8),Sigma=diag(c(0.02,0.03)))
points(clust3,col="lightcoral",pch=19,cex=0.9)

forest3=rbind(cbind(clust1,-1),cbind(clust2,1),cbind(clust3,-1))
if(length(which(forest3[,2]<0))>0){forest3=forest3[-which(forest3[,2]<0),]}
if(length(which(forest3[,1]<0))>0){traub3=forest3[-which(forest3[,1]<0),]}
if(length(which(forest3[,2]>1))>0){forest3=forest3[-which(forest3[,2]>1),]}
if(length(which(forest3[,1]>1))>0){forest3=forest3[-which(forest3[,1]>1),]}


write.table(forest3,file="forest3.txt",col.names=FALSE,row.names=FALSE)


##TRAIN SET##
forest1=read.table(file="forest1.txt")
forest2=read.table(file="forest2.txt")
forest3=read.table(file="forest3.txt")

par(mfrow=c(1,3));par(pty="s")
plot(NA,xlab="x",ylab="y",xlim=c(0,1),ylim=c(0,1));title("Forest 1")
points(forest1[forest1[,3]==-1,1:2],pch=19,col="lightcoral",cex=0.9)
points(forest1[forest1[,3]==1,1:2],pch=19,col="darkseagreen",cex=0.9)

plot(NA,xlab="x",ylab="y",xlim=c(0,1),ylim=c(0,1));title("Forest 2")
points(forest2[forest2[,3]==-1,1:2],pch=19,col="lightcoral",cex=0.9)
points(forest2[forest2[,3]==1,1:2],pch=19,col="darkseagreen",cex=0.9)

plot(NA,xlab="x",ylab="y",xlim=c(0,1),ylim=c(0,1));title("Forest 3")
points(forest3[forest3[,3]==-1,1:2],pch=19,col="lightcoral",cex=0.9)
points(forest3[forest3[,3]==1,1:2],pch=19,col="darkseagreen",cex=0.9)

###TEST SET###
set.seed(23)
n=20
coords=cbind(runif(n,0,1),runif(n,0,1))
klasse=giklasse(coords)
par(mfrow=c(1,3));par(pty="s")
plot(NA,xlab="x",ylab="y",xlim=c(0,1),ylim=c(0,1));title("Forest 1")
points(coords[klasse==-1,],col="lightcoral",pch=19,cex=0.9)
points(coords[klasse==1,],col="darkseagreen",pch=19,cex=0.9)

ind_na=which(is.na(klasse)==TRUE)
seeds1=cbind(coords[-ind_na,],klasse[-ind_na])
write.table(seeds1,file="seeds1.txt",col.names=FALSE,row.names=FALSE)

plot(NA,xlab="x",ylab="y",xlim=c(0,1),ylim=c(0,1));title("Forest 2")
coords=cbind(runif(n,0,1),runif(n,0,1))
klasse=giklasse2(coords)
if(length(klasse)!=dim(coords)[1]){klasse[length(klasse)+1]=NA}

points(coords[klasse==-1,],col="lightcoral",pch=19,cex=0.9)
points(coords[klasse==1,],col="darkseagreen",pch=19,cex=0.9)

seeds2=cbind(coords,klasse)
write.table(seeds2,file="seeds2.txt",col.names=FALSE,row.names=FALSE)

plot(NA,xlab="x",ylab="y",xlim=c(0,1),ylim=c(0,1));title("Forest 3")
clust1=mvrnorm(7,mu=c(0.2,0.1),Sigma=diag(c(0.05,0.07)))
points(clust1,col="lightcoral",pch=19,cex=0.9)
clust2=mvrnorm(5,mu=c(0.5,0.5),Sigma=diag(c(0.02,0.02)))
points(clust2,col="darkseagreen",pch=19,cex=0.9)
clust3=mvrnorm(7,mu=c(0.9,0.8),Sigma=diag(c(0.02,0.03)))
points(clust3,col="lightcoral",pch=19,cex=0.9)

seeds3=rbind(cbind(clust1,-1),cbind(clust2,1),cbind(clust3,-1))
if(length(which(seeds3[,2]<0))>0){seeds3=seeds3[-which(seeds3[,2]<0),]}
if(length(which(seeds3[,1]<0))>0){seeds3=seeds3[-which(seeds3[,1]<0),]}
if(length(which(seeds3[,2]>1))>0){seeds3=seeds3[-which(seeds3[,2]>1),]}
if(length(which(seeds3[,1]>1))>0){seeds3=seeds3[-which(seeds3[,1]>1),]}


write.table(seeds3,file="seeds3.txt",col.names=FALSE,row.names=FALSE)


#All together:
forest1=read.table(file="forest1.txt"); seeds1=read.table(file="seeds1.txt")
forest2=read.table(file="forest2.txt");seeds2=read.table(file="seeds2.txt")
forest3=read.table(file="forest3.txt");seeds3=read.table(file="seeds3.txt")

par(mfrow=c(1,3));par(pty="s")
plot(NA,xlab="x",ylab="y",xlim=c(0,1),ylim=c(0,1));title("Forest 1")
points(forest1[forest1[,3]==-1,1:2],pch=19,col="lightcoral",cex=0.9)
points(forest1[forest1[,3]==1,1:2],pch=19,col="darkseagreen",cex=0.9)
points(seeds1[seeds1[,3]==-1,1:2],pch=21,col="black",bg="lightcoral",cex=1.2)
points(seeds1[seeds1[,3]==1,1:2],pch=21,col="black",bg="darkseagreen",cex=1.2)

plot(NA,xlab="x",ylab="y",xlim=c(0,1),ylim=c(0,1));title("Forest 2")
points(forest2[forest2[,3]==-1,1:2],pch=19,col="lightcoral",cex=0.9)
points(forest2[forest2[,3]==1,1:2],pch=19,col="darkseagreen",cex=0.9)
points(seeds2[seeds2[,3]==-1,1:2],pch=21,col="black",bg="lightcoral",cex=1.2)
points(seeds2[seeds2[,3]==1,1:2],pch=21,col="black",bg="darkseagreen",cex=1.2)

plot(NA,xlab="x",ylab="y",xlim=c(0,1),ylim=c(0,1));title("Forest 3")
points(forest3[forest3[,3]==-1,1:2],pch=19,col="lightcoral",cex=0.9)
points(forest3[forest3[,3]==1,1:2],pch=19,col="darkseagreen",cex=0.9)
points(seeds3[seeds3[,3]==-1,1:2],pch=21,bg="lightcoral",col="black",cex=1.2)
points(seeds3[seeds3[,3]==1,1:2],pch=21,bg="darkseagreen",col="black",cex=1.2)

