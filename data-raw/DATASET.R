## code to prepare `DATASET` dataset goes here

usethis::use_data(DATASET, overwrite = TRUE)

simulate_data <- function() {

  #int
  gamma00 <- 10

  #mean salary
  gamma01 <- 1

  #mean control
  gamma02 <- 1

  #leadership
  gamma03 <- 1

  #centered salary
  gamma10 <- 1

  #centered control
  gamma20 <- 1


  tau00 <- 10
  tau11 <- 1
  tau22 <- 1
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
  #population decomposition

  teachsat <- as.data.frame(matrix(NA,clusters*clustersize,14))
  colnames(teachsat) <- c("schoolID","teacherID","satisfaction","control_c","salary_c","control_m","salary_m","leadership","u0j","u1j","u2j","u3j","eij","type")

  teachsat[,"schoolID"] <- rep(seq(clusters),each = clustersize)
  teachsat[,"teacherID"] <- rep(seq(clustersize),clusters)

  #generate predictors
  teachsat[,"control_c"] <- rnorm(clusters*clustersize,0,1)
  teachsat[,"control_m"] <- rep(rnorm(clusters,0,1),each = clustersize)
  teachsat[,"salary_c"] <- rnorm(clusters*clustersize,0,2)
  teachsat[,"salary_m"] <- rep(rnorm(clusters,50,2),each = clustersize)
  teachsat[,"leadership"] <- rep(rnorm(clusters,0,1),each = clustersize)



  #generate errors
  teachsat[,"eij"] <- rnorm(clusters*clustersize, 0, sqrt(sigma2))
  randomeffects <- MASS::mvrnorm(clusters, c(0,0,0), Tau)
  teachsat[,"u0j"] <- rep(randomeffects[,1], each = clustersize)
  teachsat[,"u1j"] <- rep(randomeffects[,2], each = clustersize)
  teachsat[,"u2j"] <- rep(randomeffects[,3], each = clustersize)

  teachsat <- rockchalk::gmc(teachsat, c("control_c","salary_c"), "schoolID")

  teachsat$control_c <- teachsat$control_c_dev
  teachsat$salary_c <- teachsat$salary_c_dev

  teachsat <- subset(teachsat, select = -c(control_c_dev,salary_c_dev,control_c_mn,salary_c_mn))

  #generate outcome
  for (i in seq(clusters*clustersize)) {
    teachsat[i,"satisfaction"] <- gamma00 + gamma01*teachsat[i,"salary_m"] + gamma02*teachsat[i,"control_m"] + gamma03*teachsat[i,"leadership"] +
      gamma10*teachsat[i,"salary_c"] + gamma20*teachsat[i,"control_c"] +
      teachsat[i,"u0j"] + teachsat[i,"u1j"]*teachsat[i,"salary_c"] + teachsat[i,"u2j"]*teachsat[i,"control_c"] + teachsat[i,"eij"]
  }

  teachsat$type <- rep(c("public","private","charter"), each = clusters*clustersize/3)

  return(teachsat)

}

teachsat = simulate_data() %>%
  dplyr::select(schoolID:leadership, type)

