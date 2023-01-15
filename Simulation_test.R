    rp <-1  
# 
    thin <- 2;STEPS <- 50#for test
#   
#   thin <- 100;STEPS <- 10^4
#   
  
  N <- 30;  Nnode <- 40;  G <-3  
  sim.data <- readRDS(paste("simdata_","Nnet",N,"_Nnode",Nnode,"_",rp,".rds",sep=''))
  packages <- c("scales",
                "foreach", "doParallel",
                "dplyr", "xtable",
                "pwr", "ergm",
                "coda", "flexclust", "ClusterR",
                "gtools", "sna","mvtnorm","mcclust","Bergm")
  sapply(packages, require, character.only = TRUE)
  source("mixture.ergm.functions.R")
  source("mixture.ergm.functions_SR2.R")
  source("mixture.ergm.functions_SR3.R")
  source("mixture.ergm.functions_SR4.R")
  source("mixture.ergm.functions_SR5.R")
  source("DPM_ERGM_functions2.R")
  source("DPM_ERGM_functions3.R")
  source("DPM_ERGM_functions4.R")
  
  alpha = 0.001; cluster_num <- G; p<-3; K<- G*2
  burn.in <- ceiling(STEPS*0.6)
  form_sim <- y ~ edges + gwesp(0.25, fixed=TRUE) + nodematch("nodeattr")
  sim.x <- sim.data$sim.x
  sim.k <- sim.data$sim.k

  output_Y3<-finite_full(sim.x=sim.x,sim.k=sim.k,N=N,Nnode=Nnode,alpha=alpha,K=K,STEPS=STEPS,burn.in = burn.in,thin = thin,post_process=TRUE)
  saveRDS(output_Y3,paste("Output3/","Y3_","Nnet",N,"_Nnode",Nnode,"_",rp,".rds",sep=''))


