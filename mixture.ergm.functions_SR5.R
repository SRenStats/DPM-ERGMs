finite_full <- function(sim.x=sim.x,sim.k=sim.k,N=N,Nnode=Nnode,alpha=alpha,K=K,STEPS=STEPS,burn.in = burn.in,thin = thin,post_process=post_process){
  samp.start.time <- Sys.time()
  order.constraints = TRUE
  sim.ego.terms<-c('edges','gwesp(0.25, fixed = TRUE)','nodematch("nodeattr")')
  Nterms <- length(sim.ego.terms)
  p <- Nterms
  ergmformula <- paste("~", paste(sim.ego.terms,collapse="+"),sep="")
  n.statistics <- matrix(nrow=N,ncol=Nterms)
  for (i in 1:N){
    n.statistics[i,] <- summary(as.formula(paste("sim.x[[i]]",ergmformula)))
  }
  MLE <- matrix(nrow=N,ncol=Nterms)
  for (i in 1:N){
    MLE[i,] <- ergmAPL(formula = as.formula(paste("sim.x[[i]]",ergmformula)),estimate = "CD" )$Theta_MLE
  }

  int.stat <- mean(c(max(n.statistics[,1]),min(n.statistics[,1])))
  g.init <- sim.x[[which.min(abs(n.statistics[,1]- int.stat))]]
  mu_0 <- c(-1,rep(0,Nterms-1))#prior mean in Yin
  sigma_0 <- diag(25,nrow=Nterms) #prior variance in Yin
  sigma_q <- diag(0.05^2,nrow=Nterms)#variance of proposal 
  ###pre-specified value, MCMH for mixture
  nan <- 10# number of auxiliary networks 
  nsteps1 <- 2
  #auxiliary_path1 <- array(dim=c(nan,Nterms,nsteps1+2))
  nsteps2 <- 5
  LOWESTLL=-1e8
  tau <- matrix(NA, nrow = 1+STEPS, ncol = K)
  theta <- array(NA, dim = c(1+STEPS, p, K))
  Z <- matrix(NA, nrow = 1+STEPS, ncol = N)
  acc.counts <- matrix(FALSE, nrow=1+STEPS, ncol = K)
  K_nonzero <- rep(NA, 1+STEPS) # keep track of the number of non-zero components
  
  #Z[1,] <- sample(x=1:K,size=N,replace=TRUE,prob=rep(alpha[1]/(alpha[1] * K), K))
  kmeans_ini <- kmeans(MLE,K)
  if (nrow(kmeans_ini$centers)==K){
    Z[1,] <-  kmeans_ini$cluster
    theta[1,,] <- t(kmeans_ini$centers)
  }else{
    Z[1,] <- sample(x=1:K,size=n,replace=TRUE,prob=rep(alpha[1]/(alpha[1] * K), K))
    for(k in 1:K){
      theta[1,,k] <- runif(n=p, min=-0.1,max=0.1)
    }
  }

  
  tau[1,] <- rdirichlet(n=1, alpha = rep(alpha[1], K)) 
  K_nonzero[1] <- length(table(Z[1,]))
  
  j = 2
  while(j <= 1+STEPS){
    # -- 1. update tau --#
    # update tau (from dirichlet), now we have Z[j,]
    tau[j,] <- gtools::rdirichlet(n=1, 
                                  alpha = (rep(alpha[1],K) + table(factor(Z[j-1,], levels = 1:K))) ) # K-dimensional vector
    
    # -- 2. update theta --#
    nK <- sample(1:K, K)
    NCR_path1 <- NULL
    for (k in nK){ # different clusters
      theta1 <- as.vector(rmvnorm(1,mean=theta[j-1,,k], sigma = sigma_q))
      pr <- dmvnorm( rbind(theta1, theta[j-1,,k]), mean=mu_0,sigma=sigma_0,log=TRUE)
      b <- which(Z[j-1,] == k) 
      if (length(b)>=1){
        if (length(b)>1){
          sy <- colSums(n.statistics[b,])
        }else{
          sy <- n.statistics[b,]
        }
        AC1 <- ( theta1-theta[j-1,,k]) %*% sy+pr[1]-pr[2]
        path <- sapply(seq(from = 0, to = 1 , length.out = nsteps1+2),function(u) cbind(theta1* (1-u) + theta[j-1,,k] * u))
        auxiliary_path1 <- lapply(2:(nsteps1+2),function(ji) simulate(as.formula(paste("g.init",ergmformula)),nsim=nan,coef=path[,ji],output='stats'))
        for (kj in 2:(nsteps1+2)){
          path_mid1 <- auxiliary_path1[[kj-1]]%*%(path[,kj-1]-path[,kj])
          NCR_path1[kj-1] <- log(mean(exp(path_mid1-max(path_mid1))))+max(path_mid1)
        }
        LNCR_g1<- sum(NCR_path1)
        AC2 <- length(b)*LNCR_g1
        AC <- AC1-AC2
      }else{
        AC <- pr[1]-pr[2]
      }
      if (AC >= log(runif(1)))
      { theta[j,,k] <- theta1
      acc.counts[j,k] <- TRUE
      }else {
        theta[j,,k] <- theta[j-1,,k]
      }
    }
    
    # 3. update the assignment label (total of N labels) Z
    LNCR_g<-rep(0,K) #log normalizing constant ratio among groups
    to <- c(-2,rep(0,Nterms-1))
    NCR_path2 <- matrix(1,nrow=K, ncol=nsteps2+1)
    for (jj in 1:K){
      from <- theta[j,,jj]
      path <-sapply(seq(from = 0 , to = 1, length.out = nsteps2+2),function(u) cbind(from * (1 - u)+to * u ))
      auxiliary_path2 <- lapply(2:(nsteps2+2),function(ji) simulate(as.formula(paste("g.init",ergmformula)),nsim=nan,coef=path[,ji],output='stats'))
      for (k in 2:(nsteps2+2)){
      path_mid <- auxiliary_path2[[k-1]]%*%(path[,k-1]-path[,k])
      NCR_path2[jj,k-1] <- log(mean(exp(path_mid-max(path_mid))))+max(path_mid)
      }
      LNCR_g[jj] <- sum(NCR_path2[jj,])
    }
    #posterior probability for every sample in every group
    lambda <- matrix(nrow=N,ncol=K)#probability
    for (i in 1:N){
      for (k in 1:K){
          lambda[i,k] <-log(tau[j,k])+theta[j,,k] %*% n.statistics[i,]- LNCR_g[k]
      }
    }
    Z[j,] <- apply(lambda,1,which.max)
    

    # I am not sure whether this is necessary, 10/19
    if(order.constraints){
      theta_order <- order(theta[j,1,], decreasing = FALSE)
      theta[j,,] <- theta[j,,theta_order]
      tau[j,] <- tau[j,theta_order]
      # also need to permute the latent cluster assignment...
      temp_Z <- data.frame(old_label = Z[j,])
      reference_Z <- data.frame(old_label=1:K, correct_label = theta_order)
      correct_Z <- dplyr::left_join(temp_Z, reference_Z, by = c("old_label" = "old_label"))
      Z[j,] <- correct_Z$correct_label
    } else{
      # randomly permute their labels
      theta_order <- sample(1:K,size=K,replace = FALSE)
      theta[j,,theta_order] <- theta[j,,]
      tau[j,theta_order] <- tau[j,]
      # also need to permute the latent cluster assignment...
      temp_Z <- data.frame(old_label = Z[j,])
      reference_Z <- data.frame(old_label=1:K, correct_label = theta_order)
      correct_Z <- dplyr::left_join(temp_Z, reference_Z, by = c("old_label" = "old_label"))
      Z[j,] <- correct_Z$correct_label
    }
    
    # update K_nonzero: number of components with at least one observation
    K_nonzero[j] <- length(table(Z[j,])) # only clusters with observations are counted
    # update j
    j <- j+1
  }
  
  
  K_nonzero.thin = K_nonzero[seq(burn.in+1, STEPS+1, by=thin)]
  tau.thin = tau[seq(burn.in+1, STEPS+1,by=thin),]
  Z.thin = Z[seq(burn.in+1, STEPS+1,by=thin),]
  theta.thin = theta[seq(burn.in+1, STEPS+1,by=thin),,]
  
  K_est_m <- mean(K_nonzero.thin)
  RI <- apply(apply(Z.thin,1,comPart,sim.k),1,mean)
  VI <- mean(apply(Z.thin,1,vi.dist,sim.k))
  
  Khat_table = table(factor(K_nonzero.thin, 
                            levels=1:K))
  Khat = as.numeric(which(Khat_table == max(Khat_table))[1])
  ####
  ####
  #### K-centroid clustering 
  # first identify those iterations with the number of non-zero components = Khat
  if(post_process){
    sub_id <- which(K_nonzero.thin == Khat)
    # keep those after burn_in
    # sub_id <- sub_id[sub_id >= (burn.in+1)]
    M0 <- length(sub_id)
    # pull out the theta.thin's to form a matrix : p columns, Khat * M0 rows
    theta_sub <- matrix(NA, nrow = M0*Khat, ncol = p)
    tau_sub <- rep(NA, length=M0*Khat)
    # browser()
    # alpha_sub <- rep(NA, length=M0)
    for(jj in 1:M0){
      if(p == 1){
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
      theta_sub_sub <-  array(NA, dim = c(M0, p, Khat))
      tau_sub_sub <- matrix(NA, nrow = M0, ncol = Khat)
      alpha_sub_sub <- rep(NA, length = M0)
      M0_sub = 0 # number of iterations retained 
      ###########
      for(jj in 1:M0){
        # if the cluster assignment is a permutation
        if(all(sort(kcca_rlt$clusters[(1+(jj-1)*Khat):(jj*Khat)]) == c(1:Khat)) ){
          if(p == 1){
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
      theta_sub_sub <-  array(NA, dim = c(M0, p, Khat))
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
        if(p == 1){
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
    # class_prob_PL <- function(mplesetup_list, eta, theta, samp.size,...)
    post_class_prob <- NA
    post_class_label <- rep(NA, length = N)
    #browser()
    post_class_prob <-  class_prob_infi_full(n.statistics = n.statistics,
                                             tau = tau_sub_sub,
                                             theta = theta_sub_sub,
                                             g.init=g.init,
                                             ergmformula=ergmformula,
                                             nan=nan,
                                             nsteps2=nsteps2)
    ###########################
    post_class_label <- apply(post_class_prob,1,which.max)
    set.seed(NULL)
  }else{post_class_label <- Z[j-1,]
  tau_sub_sub <- NULL
  theta_sub_sub <- NULL
  }
  
  RI2 <- comPart(post_class_label,sim.k)
  VI2 <- vi.dist(post_class_label,sim.k)
  samp.end.time <- Sys.time()
  time_used <-samp.end.time-samp.start.time
  
  output1 <- rbind(c(K_est_m,Khat,RI,VI,time_used),c(K_est_m,Khat,RI2,VI2,time_used))
  output2 <- list(K_est=K_nonzero,z_est=Z,theta_est=theta,w_est=tau, z_est2 = post_class_label,w_est2 = tau_sub_sub, 
                  theta_est2= theta_sub_sub,time_used=time_used)
  output_label <- rbind(post_class_label,sim.k)
  cat("finite_true","K_est_m",K_est_m,"K_mode",Khat,"RI",RI,"VI",VI, "\n")
  cat("finite_true","K_est_m",K_est_m,"K_mode",Khat,"RI",RI2,"VI",VI2, "\n")
  cat("time_used", time_used,"\n")
  list(output1=output1,output_label=output_label,output2=output2)
  
}
