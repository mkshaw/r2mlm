## code to prepare `teachsat` dataset goes here

usethis::use_data(teachsat, overwrite = TRUE)

#load things
library(rockchalk)
library(tidyverse)
library(lme4)
library(nlme)
library(lmerTest)
library(MASS)
library(truncnorm)

#satisfaction <- gamma00 + gamma01*salary_m + gamma02*control_m + gamma03*s_t_ratio +
#  gamma10*salary_c + gamma20*control_c + u0j + u1j*salary_c + u2j*control_c + eij

#fixed component of intercept
gamma00 <- 10

#fixed component of school mean salary
gamma01 <- .5

#fixed component of school mean control
gamma02 <- .5

#fixed component of student-teacher ratio
gamma03 <- -.5

#fixed component of school-mean-centered salary
gamma10 <- 1

#fixed component of school-mean-centered control
gamma20 <- 4

##random effect (co)variances
tau00 <- 75
tau11 <- .5
tau22 <- 6
tau01 <- 0
tau02 <- 0
tau12 <- 0
Tau <- matrix(c(tau00,tau01,tau02,
                tau01,tau11,tau12,
                tau02,tau12,tau22),3,3)
sigma2 <- 100

#set sample size
clusters <- 300
clustersize <- 30

##make dataset
teachsat <- as.data.frame(matrix(NA,clusters*clustersize,12))
colnames(teachsat) <- c("schoolID","teacherID","satisfaction","control_c","salary_c","control_m","salary_m","s_t_ratio","u0j","u1j","u2j","eij")

teachsat[,"schoolID"] <- rep(seq(clusters),each=clustersize)
teachsat[,"teacherID"] <- rep(seq(clustersize),clusters)

#generate predictors
teachsat[,"control_c"] <- (rtruncnorm(clusters*clustersize,a=-3,b=3,mean=0,sd=3))
teachsat[,"control_m"] <- (rep(rtruncnorm(clusters,a=3,b=7,mean=5,sd=2),each=clustersize))
teachsat[,"salary_c"] <- (rtruncnorm(clusters*clustersize,a=-20,b=20,mean=0,sd=10))
teachsat[,"salary_m"] <- rep((rtruncnorm(clusters,a=50,b=100,70,15)),each=clustersize)

teachsat[,"s_t_ratio"] <- rep(sample(c(15:50),clusters,replace=T),each=clustersize)

#generate errors
teachsat[,"eij"] <- rnorm(clusters*clustersize,0,sqrt(sigma2))
randomeffects <- mvrnorm(clusters,c(0,0,0),Tau)
teachsat[,"u0j"] <- rep(randomeffects[,1],each=clustersize)
teachsat[,"u1j"] <- rep(randomeffects[,2],each=clustersize)
teachsat[,"u2j"] <- rep(randomeffects[,3],each=clustersize)

#group-mean-center level-1 vars
teachsat<-gmc(teachsat,c("control_c","salary_c"),"schoolID")
teachsat$control_c <-teachsat$control_c_dev
teachsat$salary_c <-teachsat$salary_c_dev

#remove unnecessary vars
teachsat<- subset(teachsat, select = -c(control_c_dev,salary_c_dev,control_c_mn,salary_c_mn))

#generate outcome
for(i in seq(clusters*clustersize)){
  teachsat[i,"satisfaction"] <- gamma00 + gamma01*teachsat[i,"salary_m"]+ gamma02*teachsat[i,"control_m"]+ gamma03*teachsat[i,"s_t_ratio"]+
    gamma10*teachsat[i,"salary_c"]+ gamma20*teachsat[i,"control_c"]+
    teachsat[i,"u0j"] + teachsat[i,"u1j"]*teachsat[i,"salary_c"]+ teachsat[i,"u2j"]*teachsat[i,"control_c"]+ teachsat[i,"eij"]
}

##rescale
teachsat[,"satisfaction"] <- ((teachsat[,"satisfaction"] -mean(teachsat[,"satisfaction"] ))/sd(teachsat[,"satisfaction"]))*1.5 +6

#make sure satisfcation is bound by 1 and 10
for(i in seq(nrow(teachsat))){
  if(teachsat[i,"satisfaction"] > 10) teachsat[i,"satisfaction"]  <- 10
  if(teachsat[i,"satisfaction"] < 1) teachsat[i,"satisfaction"]  <- 1
}

#remove unnecessary vars
teachsat<-teachsat[,-c(9:12)]

