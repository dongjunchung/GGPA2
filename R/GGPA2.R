
# main GGPA2 function

GGPA2 <- function( gwasPval, pgraph=NULL, annotMat=NULL, nBurnin=10000, nMain=40000, lbPval=1e-10, verbose=1 ) {

  
	# summarizing setting for graph-GPA

	mcmcSetting <- list()
	mcmcSetting$nBurnin <- nBurnin
	mcmcSetting$nMain <- nMain
	mcmcSetting$lbPval <- lbPval

	# convert p-value to a matrix, if needed, for consistency

	if ( !is.matrix(gwasPval) ) {
		gwasPval <- as.matrix(gwasPval)
	}

	# check correctness of data
	# gwasPval: M * K matrix (p-value, [ 0, 1 ] )

	if ( any( gwasPval < 0 | gwasPval > 1 ) ) {
		stop( "Some p-values are smaller than zero or larger than one. p-value should be ranged between zero and one. Please check your p-value matrix." )
	}


	# check dimensions match between gwasPval and pgraph

	if ( !is.null(pgraph) ) {
	  if ( ncol(gwasPval) == ncol(pgraph) ) {
	    message( "Info: Prior phenotype graph is provided and will be used in the estimation." )
	    mcmcSetting$usePgraph <- TRUE
	  } else {
	    stop( "Dimensions of the GWAS p-value matrix and the prior graph matrix do not match!" )
	  }
	} else {
	  message( "Info: Uniform prior will be used for the phenotype graph in the estimation." )
	  mcmcSetting$usePgraph <- FALSE
	}


	# initialization of constants

	nBin <- nrow(gwasPval)
	nGWAS <- ncol(gwasPval)


	# report setting

	message( "Info: Number of GWAS data: ", nGWAS )


	# set zero p-values to small values to avoid log(0)

	if ( verbose >= 1 ) {
		if ( any( gwasPval < lbPval ) ) {
			message( "Info: Some SNPs have p-values close to zero." )
			message( "Info: Number of SNPs with p-values close to zero: ", length(which( gwasPval < lbPval )) )
			message( "Info: p-values for these SNPs are set to ", lbPval )

			gwasPval[ gwasPval < lbPval ] <- lbPval
		}
	}

	# probit transformation of p-values

	Y <- qnorm( 1 - gwasPval )

	Y <- t(Y)
  n_pheno <- nrow(Y)
  n_SNP <- ncol(Y)
  Varnames <- rownames(Y)

  message("")


	##############################################
	#                                            #
  #                 initialization             #
  #                                            #
	##############################################

  model = new( cGGPA2, Y )

	# initialize emission distribution

	model$mu_vec = rep( 0.4, n_pheno )
  model$sig2_vec = rep( 0.5, n_pheno )

  # initialize association matrix

  init_E = matrix( sample(x=c(0,1),size=(n_pheno*n_SNP),replace=TRUE), nrow=n_pheno )
  init_E[ Y<=0 ] <- 0 # Added for Ver_1_2_2
  model$E_mat = init_E

  # initialize the phenotype graph (25% ceiling)

  prop_nonzero <- 0.25

  init_G = diag(9, n_pheno)

  ind_UT <- as.matrix( row(init_G) < col(init_G) )
  n_UT <- length(which(ind_UT))
  n_nonzero <- ceiling( n_UT * prop_nonzero )

  init_G[ which(ind_UT)[ sample( x=seq_len(n_UT), size=n_nonzero, replace=FALSE ) ] ] <- 1

  for (i in seq_len(n_pheno-1)){
  	for (j in seq(i+1,n_pheno)){
  		init_G[j,i] = init_G[i,j]
  	}
  }

  # force in edges, if a phenotype graph is provided

  if ( !is.null(pgraph) ) {
    # symmetrize the phenotype graph matrix

    for (i in seq_len(n_pheno)){
    	for (j in seq_len(n_pheno)){
    		pgraph[i,j] = max( pgraph[i,j], pgraph[j,i] )
    	}
    }

    # force in edges

    for (i in seq_len(n_pheno)){
    	for (j in seq_len(n_pheno)){
    		init_G[i,j] = max( init_G[i,j], pgraph[i,j] )
    	}
    }
  } else {
    pgraph <- matrix( 0, n_pheno, n_pheno )
    rownames(pgraph) <- colnames(pgraph) <- colnames(gwasPval)
  }

  # initialize MRF coefficients

  init_beta = diag(-6,n_pheno)
  init_beta[ init_G == 1 ] <- 0.9

  ############
  #   added  #
  have_annot = NULL
  if (is.null(annotMat)){
    have_annot = FALSE
    n_annot = 2
    annotMat = array(0,c(n_annot,n_SNP)) # Why?
  } else {
    have_annot = TRUE
    
    annotMat = t(annotMat) # transpose annotaMat
    
    n_annot = dim(annotMat)[[1]]
  	if ( nrow(gwasPval) != ncol(annotMat) ) {
    	stop( "nrow(gwasPval) != ncol(annotMat)" )
  	}
    message( "Info: Number of Annotations: ", nrow(annotMat) )
  }	
  model$annotMat = annotMat
  model$have_annot = have_annot
  init_gamma_mat = array(0,c(n_pheno,n_annot))
  # init_gamma_mat[c(1:3),1] = init_gamma_mat[4,2] = init_gamma_mat[5,3] = init_gamma_mat[6,4] = 2 ; # Use simulation true
  model$gamma_mat = init_gamma_mat
  model$u_mat = array(0,c(n_pheno,n_annot))
  ############
  
  # Fixed Hyperparameters

  model$Beta = init_beta
  model$G_mat = init_G

  model$theta_mu = 0 ; model$tau2_mu = 10000 ; # CHECK 0.0 is ok
  model$a_sigma = 0.5 ; model$b_sigma = 0.5 ;
  model$theta_alpha = 0 ; model$tau2_alpha = 10000 ; model$stepsize_alpha = 0.1 ;
  model$a_beta = 4.0 ; model$b_beta = 2.0 ; model$stepsize_beta = 0.1 ; model$stepsize_gamma = 0.05 ; 
  model$a_betaG = 1.0 ; model$b_betaG = 1.0
  model$threshold_on = 4.3 ; # V 1.3.1  e_it = 1 if y_it > threshold_on
	
	model$a_gamma = 4.0 ; model$b_gamma = 4.0 ; 
	
  model$E_forcein_mat = pgraph

	##############################################
	#                                            #
	#                 burn-in                    #
	#                                            #
	##############################################


	  model$Initialize()

	  is_accept = array(0,c(nBurnin,9))
	  dimnames(is_accept)[[2]] = c("e_it","mu_i","sig2_ij","alpha_i","beta_ij","G_beta_ij","empty","empty","empty")
	  Count_G_ij = array(0,c(n_pheno,n_pheno))

	  draw_loglikelihood = rep(0,nBurnin)
	  draw_sum_G_ij = rep(0,nBurnin)
	  draw_beta = array(0,c(nBurnin,n_pheno,n_pheno))
	  draw_mu_vec = array(0,c(nBurnin,n_pheno))
	  draw_sig = array(0,c(nBurnin,n_pheno))
	  draw_mean_E = array(0,c(nBurnin,n_pheno))
	  draw_avg_E_ij = array(0,c(nBurnin,n_pheno))
	  draw_theta_mu = draw_tau2_mu = draw_a_beta = draw_b_beta = rep(0,nBurnin)
	  if (have_annot==TRUE){
	    draw_gamma_mat = draw_u_mat = array(0,c(nBurnin,n_pheno,n_annot))
	  }
	  # draw_normC = rep(0,nBurnin)
	  dimnames(draw_mean_E)[[2]] = Varnames

	  Starttime = Prevtime = PrevIterTime = proc.time()[3]

	  message("Burn in iterations...")
	  message("")

	  for ( temp_iter in seq_len(nBurnin) ) {

	    model$Iterate()

	  	draw_loglikelihood[temp_iter] = model$loglikelihood
	    is_accept[temp_iter,] = model$Accept
	    Count_G_ij = Count_G_ij + model$G_mat

	    draw_sum_G_ij[temp_iter] = sum(model$G_mat==1)/2
	    draw_beta[temp_iter,,] = model$Beta
	    draw_mu_vec[temp_iter,] = model$mu_vec
	    draw_sig[temp_iter,] = sqrt(model$sig2_vec)

	  	draw_mean_E[temp_iter,] = apply(model$E_mat,1,mean)

	  	draw_theta_mu[temp_iter] = model$theta_mu
	  	draw_tau2_mu[temp_iter] = model$tau2_mu
	  	draw_a_beta[temp_iter] = model$a_beta
	  	draw_b_beta[temp_iter] = model$b_beta
	  	if (have_annot==TRUE){
	  	  draw_gamma_mat[temp_iter,,] = model$gamma_mat
	  	  draw_u_mat[temp_iter,,] = model$u_mat
	  	}

	  	# draw_normC[temp_iter] = model$normC

	  	if ( verbose >= 1 ) {
	  	  if (temp_iter%%100==0) {
	  	    message( "Burn-in iteration: ",temp_iter," / ",nBurnin)

	        Currenttime = proc.time()[3]
	    		# TotalTime = Currenttime-Starttime ; AvgTime = TotalTime/temp_iter ;
	    		Last100 = Currenttime-Prevtime ; Time_to_Go = (nBurnin-temp_iter)*(Last100/100)
	    		Prevtime = Currenttime
	    		# print(paste("Elapsed=",round(TotalTime/60,1),"min. Avg for 100 iter=",round(AvgTime*100/60,1),"min, Est. Time to go=",round(Time_to_Go/60,1),"min" ))
	    		message(paste("The last 100 iter=",round(Last100/60,1),"min, Est. Time to go=",round(Time_to_Go/60,1),"min" ))
	  	  }
	  	}

	  	if ( verbose >= 3 ) {
	      if (temp_iter%%100==0) {
	    		temp_mean_E = apply(model$E_mat,1,mean)
	    		print_beta_i = diag(model$Beta)
	    		print_G = model$G_mat
	    		names(temp_mean_E) = names(print_beta_i) = Varnames
	    		dimnames(print_G)[[1]] = dimnames(print_G)[[2]] = Varnames

	    		message("acceptance rate:")
	    		print(round(apply(is_accept[seq_len(temp_iter),],2,mean),3))
	        message(sum(model$G_mat==1)/2," edges connected among possible ",n_pheno*(n_pheno-1)/2," edges")
	    		# print(paste("p(0|beta,G) = 1/C = ",round(1/model$normC,3)))
	    		message("proportion of associated SNPs:")
	    		print(round(temp_mean_E,4))
	    		message("MRF coefficient:")
	    		print(round(print_beta_i,2))
	    		message("estimated graph:")
	    		print(print_G)
	    		message("theta_mu & tau2_mu:")
	    		print(round(c(model$theta_mu,model$tau2_mu),2)) #,model$a_beta,model$b_beta),2))
	      }
	  	}

	  }

	  P_hat_ij = Count_G_ij / nBurnin
	  diag(P_hat_ij) = 0
	  dimnames(P_hat_ij)[[1]] = dimnames(P_hat_ij)[[2]] = Varnames

	  ### Final est. with s.e.

	  est_beta = sd_beta = array(0,c(n_pheno,n_pheno))
	  H95_beta = L95_beta = array(0,c(n_pheno,n_pheno))
	  est_mu_vec = sd_mu_vec = rep(0,n_pheno)
	  est_sigma1 = sd_sigma1 = rep(0,n_pheno)

	  for (i in seq_len(n_pheno)){
	    est_mu_vec[i] = mean(draw_mu_vec[,i])
	    sd_mu_vec[i] = sd(draw_mu_vec[,i])
	    est_sigma1[i] = mean(draw_sig[,i])
	    sd_sigma1[i] = sd(draw_sig[,i])
	    for (j in i:n_pheno){
	      est_beta[i,j] = mean(draw_beta[,i,j])
	      sd_beta[i,j] = sd(draw_beta[,i,j])
	      H95_beta[i,j] = quantile(draw_beta[,i,j],probs=0.975)
	      L95_beta[i,j] = quantile(draw_beta[,i,j],probs=0.025)
	    }
	  }

	  if ( verbose >= 3 ) {
	    ## is_accept ##
	    message("acceptance rate:")
	    print(round(apply(is_accept,2,mean),3))

	    ## Beta ##
	    message("MRF coefficient:")

	    dimnames(est_beta)[[1]] = dimnames(est_beta)[[2]] = Varnames
	    dimnames(sd_beta)[[1]] = dimnames(sd_beta)[[2]] = Varnames

	    est_beta_ij = est_beta
	    sd_beta_ij = sd_beta

	    print(round(est_beta_ij,2))
	    print(round(sd_beta_ij,2))

	    # Check with P_hat_ij
	    message("P_ij:")

	    print(P_hat_ij)

	    ## mu_vec ##
	    message("mu:")

	    MU = round(cbind(est_mu_vec,sd_mu_vec),2)
	    rownames(MU) = Varnames
	    print(MU)

	    ## sigma1 ##
	    message("sigma:")

	    SIGMA = round(cbind(est_sigma1,sd_sigma1),2)
	    rownames(SIGMA) = Varnames
	    print(SIGMA)

	    ## draw_mean_E ##
	    message("proportion of associated SNPs:")

	    print(round(apply(draw_mean_E,2,mean),4))
	    print(round(apply(draw_mean_E,2,sd),4))

	    ### Print from iterations after burn-in
	  }

	  model$clear.sum_E_ijt()

	  is_accept = array(0,c(nMain,9))
	  dimnames(is_accept)[[2]] = c("e_it","mu_i","sig2_ij","alpha_i","beta_ij","G_beta_ij","empty","empty","empty")
	  Count_G_ij = array(0,c(n_pheno,n_pheno))

	  draw_loglikelihood = rep(0,nMain)
	  draw_sum_G_ij = rep(0,nMain)
	  draw_beta = array(0,c(nMain,n_pheno,n_pheno))
	  draw_mu_vec = array(0,c(nMain,n_pheno))
	  draw_sig = array(0,c(nMain,n_pheno))
	  draw_mean_E = array(0,c(nMain,n_pheno))
	  draw_theta_mu = draw_tau2_mu = draw_a_beta = draw_b_beta = rep(0,nMain)
	  if (have_annot==TRUE){
	    draw_gamma_mat = draw_u_mat = array(0,c(nMain,n_pheno,n_annot))
	  }
	  # draw_normC = rep(0,nMain)
	  dimnames(draw_mean_E)[[2]] = Varnames


	##############################################
	#                                            #
	  #           main MCMC iteration              #
	  #                                            #
	##############################################

	  message("")
	  message("Main MCMC iterations...")
	  message("")

	  Starttime = Prevtime = PrevIterTime = proc.time()[3]

	  for ( temp_iter in seq_len(nMain) ) {

	  	# CurIterTime = proc.time()[3]
	  	# print(paste("Iter=",temp_iter,", For each iter, take ",round((CurIterTime-PrevIterTime)/60,1)," min"))
	  	# PrevIterTime = CurIterTime

	    model$Iterate()
	  	draw_loglikelihood[temp_iter] = model$loglikelihood
	    is_accept[temp_iter,] = model$Accept
	    Count_G_ij = Count_G_ij + model$G_mat

	    draw_sum_G_ij[temp_iter] = sum(model$G_mat==1)/2
	    draw_beta[temp_iter,,] = model$Beta
	    draw_mu_vec[temp_iter,] = model$mu_vec
	    draw_sig[temp_iter,] = sqrt(model$sig2_vec)

	  	draw_mean_E[temp_iter,] = apply(model$E_mat,1,mean)

	  	draw_theta_mu[temp_iter] = model$theta_mu
	  	draw_tau2_mu[temp_iter] = model$tau2_mu
	  	draw_a_beta[temp_iter] = model$a_beta
	  	draw_b_beta[temp_iter] = model$b_beta
	  	if (have_annot==TRUE){
	  	  draw_gamma_mat[temp_iter,,] = model$gamma_mat
	  	  draw_u_mat[temp_iter,,] = model$u_mat
	  	}
	  	
	  	# draw_normC[temp_iter] = model$normC

	  	if ( verbose >= 1 ) {
	  	  if (temp_iter%%100==0) {
	  	    message( "Main iteration: ",temp_iter," / ",nMain)

	        Currenttime = proc.time()[3]
	    		# TotalTime = Currenttime-Starttime ; AvgTime = TotalTime/temp_iter ;
	    		Last100 = Currenttime-Prevtime ; Time_to_Go = (nMain-temp_iter)*(Last100/100)
	    		Prevtime = Currenttime
	    		# print(paste("Elapsed=",round(TotalTime/60,1),"min. Avg for 100 iter=",round(AvgTime*100/60,1),"min, Est. Time to go=",round(Time_to_Go/60,1),"min" ))
	    		message(paste("The last 100 iter=",round(Last100/60,1),"min, Est. Time to go=",round(Time_to_Go/60,1),"min" ))
	  	  }
	  	}

	  	if ( verbose >= 3 ) {
	      if (temp_iter%%100==0) {
	    		temp_mean_E = apply(model$E_mat,1,mean)
	    		print_beta_i = diag(model$Beta)
	    		print_G = model$G_mat
	    		names(temp_mean_E) = names(print_beta_i) = Varnames
	    		dimnames(print_G)[[1]] = dimnames(print_G)[[2]] = Varnames

	    		message("acceptance rate:")
	    		print(round(apply(is_accept[seq_len(temp_iter),],2,mean),3))
	        message(sum(model$G_mat==1)/2," edges connected among possible ",n_pheno*(n_pheno-1)/2," edges")
	    		# print(paste("p(0|beta,G) = 1/C = ",round(1/model$normC,3)))
	    		message("proportion of associated SNPs:")
	    		print(round(temp_mean_E,4))
	    		message("MRF coefficient:")
	    		print(round(print_beta_i,2))
	    		message("estimated graph:")
	    		print(print_G)
	    		message("theta_mu & tau2_mu:")
	    		print(round(c(model$theta_mu,model$tau2_mu),2)) #,model$a_beta,model$b_beta),2))
	      }
	  	}

	  }

	  P_hat_ij = Count_G_ij / nMain
	  diag(P_hat_ij) = 0
	  dimnames(P_hat_ij)[[1]] = dimnames(P_hat_ij)[[2]] = Varnames

	  ### Final est. with s.e.
	  est_beta = sd_beta = array(0,c(n_pheno,n_pheno))
	  H95_beta = L95_beta = array(0,c(n_pheno,n_pheno))
	  est_mu_vec = sd_mu_vec = rep(0,n_pheno)
	  est_sigma1 = sd_sigma1 = rep(0,n_pheno)
	  
	  

	  for (i in seq_len(n_pheno)){
	    est_mu_vec[i] = mean(draw_mu_vec[,i])
	    sd_mu_vec[i] = sd(draw_mu_vec[,i])
	    est_sigma1[i] = mean(draw_sig[,i])
	    sd_sigma1[i] = sd(draw_sig[,i])
	    for (j in i:n_pheno){
	      est_beta[i,j] = mean(draw_beta[,i,j])
	      sd_beta[i,j] = sd(draw_beta[,i,j])
	      H95_beta[i,j] = quantile(draw_beta[,i,j],probs=0.975)
	      L95_beta[i,j] = quantile(draw_beta[,i,j],probs=0.025)
	    }
	  }
	  
	  ### Final est. with s.e. of Gammas
	  
	  if (have_annot==TRUE) {
	    nAnnot = nrow(annotMat)
	    est_gamma = sd_gamma = array(0,c(n_pheno,nAnnot))
	    
	    for (i in seq_len(n_pheno)){
	      for (j in seq(nAnnot)){
	        est_gamma[i,j] = mean(draw_gamma_mat[,i,j])
	        sd_gamma[i,j] = sd(draw_gamma_mat[,i,j])
	      }
	    }
	    
	  }

	  if ( verbose >= 3 ) {
	    ## is_accept ##
	    message("acceptance rate:")
	    print(round(apply(is_accept,2,mean),3))

	    ## Beta ##
	    message("MRF coefficient:")

	    dimnames(est_beta)[[1]] = dimnames(est_beta)[[2]] = Varnames
	    dimnames(sd_beta)[[1]] = dimnames(sd_beta)[[2]] = Varnames

	    est_beta_ij = est_beta
	    sd_beta_ij = sd_beta

	    print(round(est_beta_ij,2))
	    print(round(sd_beta_ij,2))

	    # Check with P_hat_ij
	    message("P_ij:")

	    print(P_hat_ij)

	    ## mu_vec ##
	    message("mu:")

	    MU = round(cbind(est_mu_vec,sd_mu_vec),2)
	    rownames(MU) = Varnames
	    print(MU)

	    ## sigma1 ##
	    message("sigma:")

	    SIGMA = round(cbind(est_sigma1,sd_sigma1),2)
	    rownames(SIGMA) = Varnames
	    print(SIGMA)

	    ## draw_mean_E ##
	    message("proportion of associated SNPs:")

	    print(round(apply(draw_mean_E,2,mean),4))
	    print(round(apply(draw_mean_E,2,sd),4))

	    ### Print from iterations after main MCMC
	  }

	  ## avg_prob_e_ijt
	  est_prob_e_ijt = sd_prob_e_ijt = array(0,c(n_pheno,n_pheno))
	  dimnames(est_prob_e_ijt)[[1]] = dimnames(est_prob_e_ijt)[[2]] = dimnames(sd_prob_e_ijt)[[1]] = dimnames(sd_prob_e_ijt)[[2]] = Varnames

	  for (i_pheno in seq_len(n_pheno)){
	  	for (j_pheno in i_pheno:n_pheno){
	  		est_prob_e_ijt[i_pheno,j_pheno] = mean(model$sum_E_ijt[i_pheno,j_pheno,]/nMain)
	  		est_prob_e_ijt[j_pheno,i_pheno] = est_prob_e_ijt[i_pheno,j_pheno]
	  		sd_prob_e_ijt[i_pheno,j_pheno] = sd(model$sum_E_ijt[i_pheno,j_pheno,]/nMain)
	  		sd_prob_e_ijt[j_pheno,i_pheno] = sd_prob_e_ijt[i_pheno,j_pheno]
	  	}
	  }
	  round(diag(est_prob_e_ijt),4)
	  round(apply(draw_mean_E,2,mean),4)

	  round(diag(sd_prob_e_ijt),4)

	  Sum_E_ijt = model$sum_E_ijt

	  ###### Draw graph #######

	  for (i in seq_len(n_pheno-1)){
	  	for (j in (i+1):n_pheno){
	  		P_hat_ij[j,i] = P_hat_ij[i,j]
	  	}
	  }
	  diag(P_hat_ij)=0
	  dimnames(P_hat_ij)[[1]] = dimnames(P_hat_ij)[[2]] = Varnames


	##############################################
	#                                            #
	#               result summary               #
	#                                            #
	##############################################


  if (have_annot==TRUE){
    mcmcResult <- list(
      loglik = draw_loglikelihood,
      mu = draw_mu_vec,
      sigma = draw_sig,
      mean_E = draw_mean_E,
      G = draw_sum_G_ij,
      beta = draw_beta,
      theta_mu = draw_theta_mu,
      tau2_mu = draw_tau2_mu,
      a_beta = draw_a_beta,
      b_beta = draw_b_beta, 
      draw_gamma_mat = draw_gamma_mat,
      draw_u_mat = draw_u_mat
      # draw_normC = draw_normC
    )
    
    mcmcSummary <- list(
      P_hat_ij = P_hat_ij,
      Sum_E_ijt = Sum_E_ijt,
      est_beta = est_beta,
      sd_beta = sd_beta,
      est_mu_vec = est_mu_vec,
      sd_mu_vec = sd_mu_vec,
      est_sigma1 = est_sigma1,
      sd_sigma1 = sd_sigma1,
      est_prob_e_ijt = est_prob_e_ijt,
      sd_prob_e_ijt = sd_prob_e_ijt,
      est_gamma = est_gamma,
      sd_gamma = sd_gamma
    ) 
    
    
	 } else { 
	   mcmcResult <- list(
	     loglik = draw_loglikelihood,
	     mu = draw_mu_vec,
	     sigma = draw_sig,
	     mean_E = draw_mean_E,
	     G = draw_sum_G_ij,
	     beta = draw_beta,
	     theta_mu = draw_theta_mu,
	     tau2_mu = draw_tau2_mu,
	     a_beta = draw_a_beta,
	     b_beta = draw_b_beta
	     # draw_normC = draw_normC
	   )
	   
	   mcmcSummary <- list(
	     P_hat_ij = P_hat_ij,
	     Sum_E_ijt = Sum_E_ijt,
	     est_beta = est_beta,
	     sd_beta = sd_beta,
	     est_mu_vec = est_mu_vec,
	     sd_mu_vec = sd_mu_vec,
	     est_sigma1 = est_sigma1,
	     sd_sigma1 = sd_sigma1,
	     est_prob_e_ijt = est_prob_e_ijt,
	     sd_prob_e_ijt = sd_prob_e_ijt
	   ) 
	 }
	   
   # check names of phenotypes and annotations
	 if (is.null(colnames(gwasPval))) {
	   colnames(gwasPval) <- paste0("Pheno ",seq(nGWAS))
	 } 
	 
	 annotMat =  t(annotMat)
	 
	 if (is.null(colnames(annotMat)) & have_annot ) {
	   colnames(annotMat) <- paste0("Annot ",seq(nAnnot))
	 }
	  
	  

	# return object by creating GGPAannot class object

	  new( "GGPA2", fit = mcmcResult, summary = mcmcSummary, setting = mcmcSetting, gwasPval = gwasPval, pgraph = pgraph, annotMat = annotMat)
		
}
