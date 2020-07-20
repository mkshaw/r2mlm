## code to prepare `DATASET` dataset goes here

usethis::use_data(DATASET, overwrite = TRUE)

#load things
library(rockchalk)
library(tidyverse)
library(lme4)
library(nlme)
library(lmerTest)
library(MASS)

#fixed component of intercept
gamma00 <- 20

#fixed component of school-mean salary
gamma01 <- 1.3

#fixed component of school-mean control
gamma02 <- 4

#fixed slope of s_t_ratio
gamma03 <- -.4

#fixed component of school-mean-centered salary
gamma10 <- 1.5

#fixed component of school-mean-centered control
gamma20 <- 2.5

#variance components
tau00 <- 20
tau11 <- 1
tau22 <- 4
tau01 <- 0
tau02 <- 0
tau12 <- 0
Tau <- matrix(c(tau00,tau01,tau02,
                tau01,tau11,tau12,
                tau02,tau12,tau22),3,3)
sigma2 <- 40

#set sample size
clusters <- 300
clustersize <- 30

#create dataset
teachsat <- as.data.frame(matrix(NA,clusters*clustersize,12))
colnames(teachsat) <- c("schoolID","teacherID","satisfaction","control_c","salary_c","control_m","salary_m","s_t_ratio","u0j","u1j","u2j","eij")

#ID variables
teachsat[,"schoolID"] <- rep(seq(clusters),each=clustersize)
teachsat[,"teacherID"] <- rep(seq(clustersize),clusters)

#generate predictors
teachsat[,"control_c"] <- rnorm(clusters*clustersize,0,1)
teachsat[,"control_m"] <- rep(rnorm(clusters,0,1),each=clustersize)
teachsat[,"salary_c"] <- rnorm(clusters*clustersize,0,2)
teachsat[,"salary_m"] <- rep(rnorm(clusters,50,2),each=clustersize)
teachsat[,"s_t_ratio"] <- rep(sample(c(10:30),clusters,replace=T),each=clustersize)

#generate errors
teachsat[,"eij"] <- rnorm(clusters*clustersize,0,sqrt(sigma2))
randomeffects <- MASS::mvrnorm(clusters,c(0,0,0),Tau)
teachsat[,"u0j"] <- rep(randomeffects[,1],each=clustersize)
teachsat[,"u1j"] <- rep(randomeffects[,2],each=clustersize)
teachsat[,"u2j"] <- rep(randomeffects[,3],each=clustersize)

#group mean center level-1 predictors
teachsat<- rockchalk::gmc(teachsat,c("control_c","salary_c"),"schoolID")
teachsat$control_c <-teachsat$control_c_dev
teachsat$salary_c <-teachsat$salary_c_dev
teachsat<- subset(teachsat, select = -c(control_c_dev,salary_c_dev,control_c_mn,salary_c_mn))

#grand mean center level-2 predictors
teachsat$salary_m <- teachsat$salary_m - mean(teachsat$salary_m)
teachsat$control_m <- teachsat$control_m - mean(teachsat$control_m)
teachsat$s_t_ratio <- teachsat$s_t_ratio - mean(teachsat$s_t_ratio)

#generate outcome
for(i in seq(clusters*clustersize)){
  teachsat[i,"satisfaction"] <- gamma00 + gamma01*teachsat[i,"salary_m"]+ gamma02*teachsat[i,"control_m"]+ gamma03*teachsat[i,"s_t_ratio"]+
    gamma10*teachsat[i,"salary_c"]+ gamma20*teachsat[i,"control_c"]+
    teachsat[i,"u0j"] + teachsat[i,"u1j"]*teachsat[i,"salary_c"]+ teachsat[i,"u2j"]*teachsat[i,"control_c"]+ teachsat[i,"eij"]
}

#remove error terms from dataset
teachsat<- subset(teachsat, select = -c(u0j,u1j,u2j,eij))

