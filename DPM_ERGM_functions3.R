logAPL <- function(theta, 
                   y, X, 
                   weights,
                   theta_MLE,
                   theta_MPLE,W, 
                   logC){
  theta_transf <- theta_MPLE + W %*% (theta - theta_MLE)
  xtheta_transf <- c(X %*% theta_transf)
  log.apl <- logC + sum(dbinom(weights * y, weights, 
                               expit(xtheta_transf), 
                               log = TRUE)) 
  return(log.apl)
}
class_prob_infi_APL <- function(obs.s = obs.s,
                            APL=APL,
                            theta = theta_sub_sub,
                            tau = tau_sub_sub){
  N <- length(obs.s)
  K <- ncol(tau)
  n<- nrow(tau)
  log_prob_array <- array(NA, dim=c(N, K, n))
  for (i in 1:N){
    for (j in 1:K){
      for (m in 1:n){
        #log_prob_array[i,j,m]<-log(tau[m,j])+pseudo.loglikelihood(obs.s[[i]],theta[m,,j])
        log_prob_array[i,j,m]<-log(tau[m,j])+logAPL(theta = theta[m,,j],
                                                    y= obs.s[[i]]$response,
                                                    X = obs.s[[i]]$predictor,
                                                    weights = obs.s[[i]]$weights,
                                                    theta_MLE = APL[[i]]$MLE,
                                                    theta_MPLE = APL[[i]]$MPLE,
                                                    W = APL[[i]]$W,
                                                    logC = APL[[i]]$logC)
      }
    }
  }
  prob_matrix <- matrix(NA, nrow = N, ncol = K)
  for(i in 1:N){
    for(j in 1:K){
      prob_matrix[i,j] <- mean(exp(log_prob_array[i,j,] - apply(matrix(log_prob_array[i,,],ncol=K),
                                                                1,
                                                                function(x) sna::logSum(x[which(x != -Inf)]) ))) 
    }}
  return(prob_matrix)
}

infinite_APL <- function(sim.x=sim.x,sim.k=sim.k,N=N,Nnode=Nnode,c=c,STEPS=STEPS,burn.in = burn.in,thin = thin,post_process=post_process){
  start_timeA <- Sys.time()
  LOWESTLL=-1e8
  ###1.Generate simulated dataset####
  sim.ego.terms<-c('edges','gwesp(0.25, fixed = TRUE)','nodematch("nodeattr")')
  Nterms <- length(sim.ego.terms)
  #attr<-sample((0:1),Nnode,replace=TRUE)##discrete attribute 0/1 #random 2
  ergmformula <- paste("~", paste(sim.ego.terms,collapse="+"),sep="")
  # sim.x<-list()
  # g.use <- network(Nnode,density=0.15,directed=FALSE)##name a network #random 3
  # set.vertex.attribute(g.use, "nodeattr", attr)
  # for (i in 1:N){
  #   sim.x[[i]]<-simulate(as.formula(paste("g.use",ergmformula)),coef=c(sim.theta[,sim.k[i]])) #no difference as only one network is generated once
  # } #random 4
  # 
  n.statistics <- matrix(nrow=N,ncol=Nterms)
  for (i in 1:N){
    n.statistics[i,] <- summary(as.formula(paste("sim.x[[i]]",ergmformula)))
  }
  obs.s<-list()#list of every ego-network with 4 objects
  APL <- list()
  for(i in 1:N){
    obs.s[[i]]<-ergmMPLE(as.formula(paste("sim.x[[i]]",ergmformula))) # pseudolikelihood statistics
    temp_ergmAPL_rlt <-ergmAPL(formula = as.formula(paste("sim.x[[i]]",ergmformula)),estimate="CD")
    W <- temp_ergmAPL_rlt$W
    logC <- temp_ergmAPL_rlt$logC
    MLE<- temp_ergmAPL_rlt$Theta_MLE
    MPLE <- temp_ergmAPL_rlt$Theta_PL
    APL[[i]] <- list(MLE=MLE,MPLE=MPLE,W=W,logC=logC)
  }
  #pairs(n.statistics,col=sim.k)
  #################2. Estimation###########################
  #set.seed(814)
  int.stat <- mean(c(max(n.statistics[,1]),min(n.statistics[,1])))
  g.init <- sim.x[[which.min(abs(n.statistics[,1]- int.stat))]]
  #mu_0 <- c(-3,rep(0,Nterms-1))#length of Nterms, prior
  #sigma_0 <- diag(16,nrow=Nterms) #dimension of Nterms, prior
  mu_0 <- c(-1,rep(0,Nterms-1))#prior mean in Yin
  sigma_0 <- diag(25,nrow=Nterms) #prior variance in Yin
  sigma_q <- diag(0.05^2,nrow=Nterms)#variance of proposal 
  ###pre-specified value, MCMH for mixture
  nan <- 10# number of auxiliary networks 
  nsteps1 <- 2
  #auxiliary_path1 <- array(dim=c(nan,Nterms,nsteps1+2))
  nsteps2 <- 5
  #auxiliary_path2 <- array(dim=c(nan,Nterms,nsteps2+2))
  z <- rep(1,N)
  theta <- matrix(0,nrow=Nterms,ncol=N)#initialization of theta
  theta[1,] <- -2 #set the first parameter -2
  v <- NULL; u <- NULL#latent variable 1
  xi <- exp(-1:-N)#deterministic decreasing sequence
  Ncy <- NULL
  acr <- rep(0,N)
  #save results
  Nc_all <- NULL
  w_all <- NULL
  z_all1 <- NULL
  z_all2 <- NULL
  theta2_all <- list()
  # ###Iteration####

  for (l in 1:STEPS){
    #1. Update u
    for (ia in 1:N){
      u[ia] <- runif(1,0,xi[z[ia]])
      Ncy[ia] <- max(which(xi>u[ia]))
    }
    Nc <- max(Ncy)#new number of groups
    Nc_all <- c(Nc_all,Nc)
    #NCR <- rep(1,Nc)
    
    
    ###1. Updata oc,m####
    oc <- list()
    m <- NULL
    for (j in 1:Nc)
    {
      oc[[j]] <- which(z==j)
      m[j] <- length(oc[[j]])
    }
    
    ##2. Update parameter theta ####
    NCR_path1 <- NULL
    for (k in 1:Nc){ # different clusters
      if (m[k]!=0){
        b <- oc[[k]]
        p1 <- NULL
        p2 <- NULL
        a <- matrix(rmvnorm(1,mean=theta[,k],sigma=sigma_q),ncol=1)
        for (ki in 1:m[k]){
          p1[ki] <- logAPL(theta = a,
                           y= obs.s[[b[ki]]]$response,
                           X = obs.s[[b[ki]]]$predictor,
                           weights = obs.s[[b[ki]]]$weights,
                           theta_MLE = APL[[b[ki]]]$MLE,
                           theta_MPLE = APL[[b[ki]]]$MPLE,
                           W = APL[[b[ki]]]$W,
                           logC = APL[[b[ki]]]$logC)
          p2[ki] <- logAPL(theta = theta[,k],
                           y= obs.s[[b[ki]]]$response,
                           X = obs.s[[b[ki]]]$predictor,
                           weights = obs.s[[b[ki]]]$weights,
                           theta_MLE = APL[[b[ki]]]$MLE,
                           theta_MPLE = APL[[b[ki]]]$MPLE,
                           W = APL[[b[ki]]]$W,
                           logC = APL[[b[ki]]]$logC)
          #p1[ki] <- pseudo.loglikelihood(obs.s[[b[ki]]],a)
          #p2[ki] <- pseudo.loglikelihood(obs.s[[b[ki]]],theta[,k])
        }
        p3 <- sum(p1)-sum(p2)
        pr1 <- dmvnorm(as.numeric(a),mean=mu_0,sigma=sigma_0)
        pr2 <- dmvnorm(as.numeric(theta[,k]),mean=mu_0,sigma=sigma_0)
        p_a <- exp(p3)*pr1/pr2
        if (p_a > runif(1)){
          acr[k] <- acr[k]+1
          theta[,k] <- a}
      }else{
        theta[, k] <- rmvnorm(1,mean=mu_0,sigma=sigma_0)
      }
    }
    
    theta2_all[[l]] <- theta
    
    ##3. Update v####
    for (j in 1:(Nc-1)){
      v[j] <- rbeta(1,(1+m[j]),(c+sum(m[(j+1):Nc])))
    }
    v[Nc] <- rbeta(1,(1+m[Nc]),c)
    
    
    ##4. Update w with new number of clusters####
    w <- rep(0,N)
    w[1]<-v[1]
    for (id in (2:Nc)){
      w[id]<-v[id]*prod(1-v[1:(id-1)])
    }
    w_all <- rbind(w_all,w)
    
    ###5.Update indicator Z ####
    #posterior probability for every sample in every group
    lambda <- matrix(nrow=N,ncol=Nc)#probability
    for (i in 1:N){
      for (j in 1:Nc){
        if (xi[j] > u[i]){
          #lambda[i,j] <-log(w[j])-log(xi[j])+pseudo.loglikelihood(obs.s[[i]],theta[,j])
          lambda[i,j] <-log(w[j])-log(xi[j])+logAPL(theta = theta[,j],
                                                    y= obs.s[[i]]$response,
                                                    X = obs.s[[i]]$predictor,
                                                    weights = obs.s[[i]]$weights,
                                                    theta_MLE = APL[[i]]$MLE,
                                                    theta_MPLE = APL[[i]]$MPLE,
                                                    W = APL[[i]]$W,
                                                    logC = APL[[i]]$logC)
          
        }else{
          lambda[i,j] <- LOWESTLL }
      }
    }
    #sample z
    z<-apply(lambda,1,which.max)
    
    #finish
    z_all1<-c(z_all1,length(unique(z)))#the number of clusters
    z_all2 <- rbind(z_all2,z)
  }

  
  K_nonzero.thin = z_all1[seq(burn.in+1, STEPS, by=thin)]
  tau.thin = w_all[seq(burn.in+1, STEPS,by=thin),]
  
  Z.thin = z_all2[seq(burn.in+1, STEPS,by=thin),]
  bi_id <- seq(burn.in+1, STEPS,by=thin)#burn in id
  theta.thin <- array(dim=c(length(bi_id),Nterms,N))
  for (id in 1:length(bi_id)){
    theta.thin[id,,] <- theta2_all[[bi_id[id]]]
  }
  
  
  
  K_est_m <- mean(K_nonzero.thin)
  RI <- apply(apply(Z.thin,1,comPart,sim.k),1,mean)
  VI <- mean(apply(Z.thin,1,vi.dist,sim.k))
  
  #the mode estimation for K
  Khat_table = table(factor(K_nonzero.thin,levels=1:max(K_nonzero.thin)))
  Khat <- as.numeric(which(Khat_table == max(Khat_table))[1])
  
  #new add in 23/10/2022
  #### K-centroid clustering
  # first identify those iterations with the number of non-zero components = Khat
  if (post_process==TRUE){
  sub_id <- which(K_nonzero.thin == Khat)
  # keep those after burn_in
  # sub_id <- sub_id[sub_id >= (burn.in+1)]
  M0 <- length(sub_id)
  # pull out the theta.thin's to form a matrix : p/Nterms columns, Khat * M0 rows
  theta_sub <- matrix(NA, nrow = M0*Khat, ncol = Nterms)
  tau_sub <- rep(NA, length=M0*Khat)
  # browser()
  # alpha_sub <- rep(NA, length=M0)
  for(jj in 1:M0){
    if(Nterms == 1){
      theta_sub[(1+(jj-1)*Khat):(jj*Khat),] <- theta.thin[sub_id[jj], ,unique(Z.thin[sub_id[jj],])]
    } else{
      theta_sub[(1+(jj-1)*Khat):(jj*Khat),] <- t(theta.thin[sub_id[jj], ,unique(Z.thin[sub_id[jj],])])
    }
    tau_sub[(1+(jj-1)*Khat):(jj*Khat)] <- tau.thin[sub_id[jj],unique(Z.thin[sub_id[jj],])]
    # alpha_sub[jj] <- alpha.thin[jj]
  }
  # k-centroid clustering
  # kcca_rlt <- kcca(Nclus, k=Khat, family=kccaFamily("kmedians"),
  #             control=list(initcent="kmeanspp"))
  # do clustering if and only if Khat > 1
  # browser()
  if(Khat > 1){
    #
    kcca_rlt <- Cluster_Medoids(theta_sub,
                                clusters = Khat,
                                distance_metric = "mahalanobis")
    #
    # kcca_rlt$clusters is a sequence of cluster membership indicators
    # check if the clustering result is valid for each iteration of parameters
    theta_sub_sub <-  array(NA, dim = c(M0, Nterms, Khat))
    tau_sub_sub <- matrix(NA, nrow = M0, ncol = Khat)
    alpha_sub_sub <- rep(NA, length = M0)
    M0_sub = 0 # number of iterations retained
    ###########
    for(jj in 1:M0){
      # if the cluster assignment is a permutation
      if(all(sort(kcca_rlt$clusters[(1+(jj-1)*Khat):(jj*Khat)]) == c(1:Khat)) ){
        if(Nterms == 1){
          theta_sub_sub[M0_sub+1,,kcca_rlt$clusters[(1+(jj-1)*Khat):(jj*Khat)]] <- (theta_sub[(1+(jj-1)*Khat):(jj*Khat),])
        } else{
          theta_sub_sub[M0_sub+1,,kcca_rlt$clusters[(1+(jj-1)*Khat):(jj*Khat)]] <- t(theta_sub[(1+(jj-1)*Khat):(jj*Khat),])
        }
        
        # alpha_sub_sub[M0_sub+1] <- alpha_sub[jj]
        tau_sub_sub[M0_sub+1,kcca_rlt$clusters[(1+(jj-1)*Khat):(jj*Khat)]] <- tau_sub[(1+(jj-1)*Khat):(jj*Khat)]
        # tau_sub_sub[jj,kcca_rlt$clusters[(1+(jj-1)*Khat):(jj*Khat)]] <- tau.thin[jj, ]
        # theta[sub_id[jj], ,unique(Z[sub_id[jj],])]
        M0_sub = M0_sub + 1
      }
    }
    
    # only keep the complete.cases
    theta_sub_sub = theta_sub_sub[1:M0_sub,,]
    # alpha_sub_sub = alpha_sub_sub[1:M0_sub]
    tau_sub_sub = tau_sub_sub[1:M0_sub,]
    #
    tau_sub_sub <- sweep(tau_sub_sub, 1, apply(tau_sub_sub,1,sum), "/")
    # normalize the tau's to make sure rowsum is 1
    #### the cluster with smaller value on edge parameter should have smaller cluster id number
    ####
    cluster_rank <- rank(kcca_rlt$medoids[,1])
    ####
    for(jj in 1:nrow(tau_sub_sub)){
      tau_sub_sub[jj,cluster_rank] <- tau_sub_sub[jj,]
      theta_sub_sub[jj,,cluster_rank] <- theta_sub_sub[jj,,]
    }
    #
  } else if(Khat == 1){
    ###########
    theta_sub_sub <-  array(NA, dim = c(M0, Nterms, Khat))
    tau_sub_sub <- matrix(NA, nrow = M0, ncol = Khat)
    alpha_sub_sub <- rep(NA, length = M0)
    M0_sub = 0 # number of iterations retained
    ###########
    for(jj in 1:M0){
      # when Khat = 1, no need to check permutation
      # if the cluster assignment is a permutation
      # if(all(sort(kcca_rlt$clusters[(1+(jj-1)*Khat):(jj*Khat)]) == c(1:Khat)) ){
      # theta_sub_sub[M0_sub+1,, kcca_rlt$clusters[(1+(jj-1)*Khat):(jj*Khat)] ] <- t(theta_sub[(1+(jj-1)*Khat):(jj*Khat),])
      # # alpha_sub_sub[M0_sub+1] <- alpha_sub[jj]
      # tau_sub_sub[M0_sub+1, kcca_rlt$clusters[(1+(jj-1)*Khat):(jj*Khat)] ] <- tau_sub[(1+(jj-1)*Khat):(jj*Khat)]
      # tau_sub_sub[jj,kcca_rlt$clusters[(1+(jj-1)*Khat):(jj*Khat)]] <- tau.thin[jj, ]
      # theta[sub_id[jj], ,unique(Z[sub_id[jj],])]
      if(Nterms == 1){
        theta_sub_sub[M0_sub+1,, 1:Khat ] <- theta_sub[(1+(jj-1)*Khat):(jj*Khat),]
      } else{
        theta_sub_sub[M0_sub+1,, 1:Khat ] <- t(theta_sub[(1+(jj-1)*Khat):(jj*Khat),])
      }
      
      # alpha_sub_sub[M0_sub+1] <- alpha_sub[jj]
      tau_sub_sub[M0_sub+1, 1:Khat ] <- tau_sub[(1+(jj-1)*Khat):(jj*Khat)]
      M0_sub = M0_sub + 1
      # }
    }
    tau_sub_sub <- sweep(tau_sub_sub, 1, apply(tau_sub_sub,1,sum), "/")
  }
  
  ############################# posterior class probability
  
  post_class_prob <- class_prob_infi_APL(obs.s = obs.s,
                                         APL=APL,
                                     tau = tau_sub_sub,
                                     theta = theta_sub_sub)
  ###########################
  post_class_label <- apply(post_class_prob,1,which.max)
  set.seed(NULL)
  } else {post_class_label <- z
  tau_sub_sub <- NULL
  theta_sub_sub <- NULL}
  
  
 
  
  RI2 <- comPart(post_class_label,sim.k)
  VI2 <- vi.dist(post_class_label,sim.k)
  end_timeA <- Sys.time()
  time_used <- end_timeA-start_timeA
  
  output1 <- rbind(c(K_est_m,Khat,RI,VI,time_used),c(K_est_m,Khat,RI2,VI2,time_used))
  output_label <- rbind(post_class_label,sim.k)
  output2 <- list(K_est=z_all1, z_est=z_all2,theta_est=theta2_all,w_est=w_all,w_est2 = tau_sub_sub,
                  theta_est2= theta_sub_sub,time_used=time_used)
  cat("SR","infinite_APL","K_est_m",K_est_m,"K_mode",Khat,"RI",RI,"VI",VI, "\n")
  cat("SR","infinite_APL","K_est_m",K_est_m,"K_mode",Khat,"RI",RI2,"VI",VI2, "\n")
  cat("time_used", time_used,"\n")
  list(output1=output1,output_label=output_label,output2=output2)
}


