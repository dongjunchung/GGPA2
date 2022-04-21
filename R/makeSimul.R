
makesimul = function( true_G, Beta_ij, Alpha_i, mu1, sigma1, true_Gamma = NULL, A_mat = NULL, n_SNP=20000 ) {
	
		# n_SNP <- 20000 ; # number of SNPs, t=1,...,n_SNP
		n_pheno = dim(true_G)[[1]] ; # number of phenotypes, i=1,...,n_pheno
		iter.gibbs <- 1000 # Gibbs interation to generate initial "e"

		true_Beta <- matrix( 0, n_pheno, n_pheno ) ; 		# MRF coefficient to generate e
		for ( i in 1:(n_pheno-1) ) {
		  for ( j in (i+1):n_pheno ) {
		    if ( true_G[i,j] > 0 ) {
		      true_Beta[i,j] <- true_Beta[j,i] <- Beta_ij[i,j]
		    } #
		  } #
		} # for (i)
		 
		diag(true_Beta) = Alpha_i
		
		emat = NULL 
		if ( (is.null(true_Gamma)) || (is.null(A_mat)) ){
			emat = Gibbs_e_it_no_ann(true_Beta, true_G, n_SNP, iter.gibbs)			
		} else { 
			emat = Gibbs_e_it_ann(true_Beta, true_G, n_SNP, iter.gibbs, true_Gamma, A_mat)
		} # if ... else ... 
		
		Y_mat <- matrix( NA, n_SNP, n_pheno )
		for ( i in 1:n_pheno ) {
		  Y_mat[ emat[i,] == 0, i ] <- rnorm( length(which( emat[i,] == 0 )), 0, 1 )
		  Y_mat[ emat[i,] == 1, i ] <- exp( rnorm( length(which( emat[i,] == 1 )), mu1[i], sigma1[i] ) )
		} # for (i)

		pmat = 1 - pnorm( Y_mat )
			
		return(list(pmat=pmat,Y_mat=Y_mat,true_E_mat=emat,true_beta=true_Beta,true_G=true_G,true_mu1=mu1,true_sigma1=sigma1,true_Gamma=true_Gamma,A_mat=A_mat)) 
		
} # makesimul




