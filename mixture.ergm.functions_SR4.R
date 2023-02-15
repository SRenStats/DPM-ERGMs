finite_APL <- function(dat=sim_networks,
                          sim.k=sim.k,
                          form=form_sim,
                          p=p,
                          K=cluster_num*2,
                          prior.mean=c(-1,0,0),
                          prior.sigma=diag(25,p),
                          alpha = alpha,
                          main.iters = main.iters,
                          sigma.epsilon = diag(rep(0.0025,p)),
                          thin = thin,
                          num.cores=1,
                          burn.in = burn.in,
                          order.constraints = order.constraints,
                          post.class.calc=post.class.calc,
                          size.offset = FALSE
                          ){
  samp.start.time <- Sys.time()
  n <- length(dat)
  netsize <- do.call("c",lapply(dat, function(x) network::network.size(x) ) )
  sy <- matrix(NA, nrow = n, ncol = p)
  mplesetup_list <- vector(mode="list", length = n)
  data.glm.initial_list <- vector(mode="list", length=n)
  
  obs.s<-list()#list of every ego-network with 4 objects
  APL <- list()
  for(i in 1:n){
    obs.s[[i]]<-ergmMPLE(as.formula(paste("sim.x[[i]]",ergmformula))) # pseudolikelihood statistics
    temp_ergmAPL_rlt <-ergmAPL(formula = as.formula(paste("sim.x[[i]]",ergmformula)),estimate = "CD" )
    W <- temp_ergmAPL_rlt$W
    logC <- temp_ergmAPL_rlt$logC
    MLE<- temp_ergmAPL_rlt$Theta_MLE
    MPLE <- temp_ergmAPL_rlt$Theta_PL
    APL[[i]] <- list(MLE=MLE,MPLE=MPLE,W=W,logC=logC)
  }
  
  tau <- matrix(NA, nrow = 1+main.iters, ncol = K)
  theta <- array(NA, dim = c(1+main.iters, p, K))
  Z <- matrix(NA, nrow = 1+main.iters, ncol = n)
  acc.counts <- matrix(FALSE, nrow=1+main.iters, ncol = K)
  K_nonzero <- rep(NA, 1+main.iters) # keep track of the number of non-zero components
  
  Z[1,] <- sample(x=1:K,size=n,replace=TRUE,prob=rep(alpha[1]/(alpha[1] * K), K))
  tau[1,] <- rdirichlet(n=1, alpha = rep(alpha[1], K)) 
  for(k in 1:K){
    theta[1,,k] <- runif(n=p, min=-0.1,max=0.1)
  }
  
  K_nonzero[1] <- length(table(Z[1,]))
  
  j = 2
  while(j <= 1+main.iters){
    # -- 1. update tau --#
    # update tau (from dirichlet), now we have Z[j,]
    tau[j,] <- gtools::rdirichlet(n=1, 
                                  alpha = (rep(alpha[1],K) + table(factor(Z[j-1,], levels = 1:K))) ) # K-dimensional vector
    
    # -- 2. update theta --#
    nK <- sample(1:K, K)
    for(k in nK){
      # proposing a new theta for k-th component
      theta1 <- rmvnorm(1,mean=theta[j-1,,k], sigma = sigma.epsilon)
      # log prior vector
      pr <- dmvnorm( rbind(theta1, theta[j-1,,k]), mean= prior.mean, sigma = prior.sigma, log=TRUE)
      k_id <- which(Z[j-1,] == k) # index of observed data whose corresponding assigned cluster is k 
      if(length(k_id) >= 1){
        ll_mat <- matrix(NA, nrow=2,ncol = length(k_id))
        for(l in 1:length(k_id)){
          #ll_mat[2,l] <- logPL(theta = theta[j-1,,k], y=mplesetup_list[[k_id[l]]]$response,
                               # X = mplesetup_list[[k_id[l]]]$predictor, weights = mplesetup_list[[k_id[l]]]$weights,
                               # size = mplesetup_list[[k_id[l]]]$size, size.offset = size.offset)

          ll_mat[2,l] <- logAPL(theta = theta[j-1,,k],
                 y= obs.s[[k_id[l]]]$response,
                 X = obs.s[[k_id[l]]]$predictor,
                 weights = obs.s[[k_id[l]]]$weights,
                 theta_MLE = APL[[k_id[l]]]$MLE,
                 theta_MPLE = APL[[k_id[l]]]$MPLE,
                 W = APL[[k_id[l]]]$W,
                 logC = APL[[k_id[l]]]$logC)
          # ll_mat[1,l] <- logPL(theta = c(theta1), y=mplesetup_list[[k_id[l]]]$response,
          #                      X = mplesetup_list[[k_id[l]]]$predictor, weights = mplesetup_list[[k_id[l]]]$weights,
          #                      size = mplesetup_list[[k_id[l]]]$size, size.offset = size.offset)
          ll_mat[1,l] <- logAPL(theta = c(theta1),
                                y= obs.s[[k_id[l]]]$response,
                                X = obs.s[[k_id[l]]]$predictor,
                                weights = obs.s[[k_id[l]]]$weights,
                                theta_MLE = APL[[k_id[l]]]$MLE,
                                theta_MPLE = APL[[k_id[l]]]$MPLE,
                                W = APL[[k_id[l]]]$W,
                                logC = APL[[k_id[l]]]$logC)
        }
        lr <- apply(ll_mat,1,sum) # row-sum
      }      else { lr <- c(0,0)}
      # log acceptance ratio, proposal is symmetric
      beta.theta <- (lr[1] - lr[2]) + (pr[1] - pr[2]) 
      if(is.nan(beta.theta)){
        beta.theta <- (pr[1] - pr[2]) 
      }
      if(beta.theta >= log(runif(1))){
        theta[j,,k] <- theta1
        acc.counts[j,k] <- TRUE
      }
      else {
        theta[j,,k] <- theta[j-1,,k]
      }
    }
    
    # 3. update the assignment label (total of n labels) Z
    for(i in 1:n){
      unnorm_p <- rep(NA, length = K)
      unnorm_logp <- rep(NA, length = K)
      for(k in 1:K){
        # calculating log pseudolikelihood
        # temp_logPL <- logPL(theta = c(theta[j,,k]), 
        #                     y=mplesetup_list[[i]]$response,
        #                     X = mplesetup_list[[i]]$predictor, 
        #                     weights = mplesetup_list[[i]]$weights,
        #                     size = mplesetup_list[[i]]$size, size.offset = size.offset)
        temp_logPL<-logAPL(theta =c(theta[j,,k]),
               y= obs.s[[i]]$response,
               X = obs.s[[i]]$predictor,
               weights = obs.s[[i]]$weights,
               theta_MLE = APL[[i]]$MLE,
               theta_MPLE = APL[[i]]$MPLE,
               W = APL[[i]]$W,
               logC = APL[[i]]$logC)
        unnorm_logp[k] <- log(tau[j,k]) + temp_logPL
      }
      lsum = sna::logSum(unnorm_logp[which(unnorm_logp!=-Inf)])
      unnorm_p <- exp(unnorm_logp-lsum)
      # Update Z...
      Z[j,i] <- sample(x=1:K, size=1, prob = unnorm_p)
      #Z[j,i] <- sample(x=1:K, size=1, prob = unnorm_logp)
      #Z[j,i] <-which.max(unnorm_p)#SR's way
      
    }
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
  
  
  K_nonzero.thin = K_nonzero[seq(burn.in+1, main.iters+1, by=thin)]
  tau.thin = tau[seq(burn.in+1, main.iters+1,by=thin),]
  Z.thin = Z[seq(burn.in+1, main.iters+1,by=thin),]
  theta.thin = theta[seq(burn.in+1, main.iters+1,by=thin),,]
  
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
  if(post.class.calc){
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
    post_class_label <- rep(NA, length = n)
    #browser()
    post_class_prob <- class_prob_infi_APL(obs.s = obs.s,
                                           APL=APL,
                                           tau = tau_sub_sub,
                                           theta = theta_sub_sub)
    ###########################
    post_class_prob <- matrix(post_class_prob, ncol = Khat)
    for(i in 1:n){
      post_class_label[i] <- which(post_class_prob[i,] == max(post_class_prob[i,]))
    }
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
  cat("Yin_APL","K_est_m",K_est_m,"K_mode",Khat,"RI",RI,"VI",VI, "\n")
  cat("Yin_APL","K_est_m",K_est_m,"K_mode",Khat,"RI",RI2,"VI",VI2, "\n")
  cat("time_used", time_used,"\n")
  list(output1=output1,output_label=output_label,output2=output2)
  
}
