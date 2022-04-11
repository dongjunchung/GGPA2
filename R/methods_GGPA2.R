
# generic methods for "GGPA2" class

setMethod(
  f = "get_fit",
  signature = "GGPA2",
  definition = function(x) x@fit
)

setMethod(
  f = "get_summary",
  signature = "GGPA2",
  definition = function(x) x@summary
)

setMethod(
  f = "get_setting",
  signature = "GGPA2",
  definition = function(x) x@setting
)

setMethod(
  f = "get_gwasPval",
  signature = "GGPA2",
  definition = function(x) x@gwasPval
)

setMethod(
  f = "get_pgraph",
  signature = "GGPA2",
  definition = function(x) x@pgraph
)

# GGPAannot model fit summary

setMethod(
    f="show",
    signature="GGPA2",
    definition=function( object ) {

    # summary of GGPAannot fit

    # constants

		nBin <- nrow(get_gwasPval(object))
		nGWAS <- ncol(get_gwasPval(object))

		# estimates

		est_mu_vec = sd_mu_vec = rep(0,nGWAS)
    est_sigma1 = sd_sigma1 = rep(0,nGWAS)

    for (i in seq_len(nGWAS)){
      est_mu_vec[i] = mean(get_fit(object)$mu[,i])
      sd_mu_vec[i] = sd(get_fit(object)$mu[,i])
      est_sigma1[i] = mean(get_fit(object)$sigma[,i])
      sd_sigma1[i] = sd(get_fit(object)$sigma[,i])
    }

    est_mean_E = colMeans2(get_fit(object)$mean_E)
    sd_mean_E = colSds(get_fit(object)$mean_E)

    MU = round(cbind(est_mu_vec,sd_mu_vec),2)
    rownames(MU) = colnames(get_gwasPval(object))
    colnames(MU) <- c( "estimate", "SE" )

    SIGMA = round(cbind(est_sigma1,sd_sigma1),2)
    rownames(SIGMA) = colnames(get_gwasPval(object))
    colnames(SIGMA) <- c( "estimate", "SE" )

    EMAT = round(cbind(est_mean_E,sd_mean_E),2)
    rownames(EMAT) = colnames(get_gwasPval(object))
    colnames(EMAT) <- c( "estimate", "SE" )

		# output

    cat( "Summary: GGPAannot model fitting results (class: GGPAannot)\n" )
    cat( "--------------------------------------------------\n" )
    cat( "Data summary:\n" )
    cat( "\tNumber of GWAS data: ", nGWAS , "\n", sep="" )
		cat( "\tNumber of SNPs: ", nBin , "\n", sep="" )
		cat( "Use a prior phenotype graph? " )
		if ( get_setting(object)$usePgraph == TRUE ) {
		  cat( "YES\n" )
		} else {
		  cat( "NO\n" )
		}
		cat( "mu\n", sep="" )
		print(MU)
		cat( "sigma\n", sep="" )
		print(SIGMA)
		cat( "Proportion of associated SNPs\n", sep="" )
		print(EMAT)
    cat( "--------------------------------------------------\n" )
  }
)

# phenotype graph

setMethod(
  f="plot",
  signature=c("GGPA2","missing"),
  definition=function( x, y=NULL, pCutoff = 0.5, betaCI = 0.95, nodesize=20, labelsize=8, textsize=12,Names=NULL,... ) {

    P_hat_ij <- get_summary(x)$P_hat_ij
    draw_beta <- get_fit(x)$beta

    # calculate posterior probabilities
    if (is.null(Names)) {
      Names = dimnames(P_hat_ij)[[1]] 
    }
    n_pheno = dim(draw_beta)[[2]]

    P_hat = P_lb_beta = edge_weights = matrix( 0, n_pheno, n_pheno )
    
    for (i in seq_len(n_pheno)){
    	for (j in seq_len(n_pheno)){
    		P_hat[i,j] = mean( draw_beta[,i,j] > 0 )
    		P_lb_beta[i,j] = quantile( draw_beta[,i,j], probs = ( 1 - betaCI ) / 2 )
    		edge_weights[i,j] = mean(draw_beta[,i,j])
    	}
    }
    dimnames(P_hat)[[1]] = dimnames(P_lb_beta)[[1]] = Names
    dimnames(P_hat)[[2]] = dimnames(P_lb_beta)[[2]] = Names

    # 
    #edge_weights = t(get_summary(fit)$est_beta)
    
    # graph construction

    adjmat <- round(P_hat,2) > pCutoff & P_lb_beta > 0
    if ( !is.null(Names) ) {
      Names <- seq_len(n_pheno)
    } else {
      rownames(adjmat) <- colnames(adjmat) <- Names
    }

    #edge_weights <- round(edge_weights[adjmat],1)
    edge_weights <- round(edge_weights*(adjmat*1),2)
    
    
    tmp <- network::network(adjmat, directed = FALSE)
    tmp %e% "weight" <- edge_weights
    #ggnet2(tmp, label = TRUE, mode = "circle",color="lightblue", 
    #       edge.label = "weight", node.size = nodesize, edge.label.size = labelsize,
    #       label.size = textsize)
    
    ggnet2(tmp, label = TRUE, mode = "circle",color="lightblue", 
           edge.label = "weight", size = nodesize, edge.label.size = labelsize,
           label.size = textsize)
    
  }
)

# local FDR

setMethod(
    f="fdr",
    signature="GGPA2",
    definition=function( object, i=NULL, j=NULL ) {
        # return marginal FDR

    # constants

		nBin <- nrow(get_gwasPval(object))
		nGWAS <- ncol(get_gwasPval(object))

		if ( is.null(i) & is.null(j) ) {
		  #message( "Info: Marginal local FDR matrix is returned." )

  		fdrmat <- matrix( NA, nBin, nGWAS )
  		colnames(fdrmat) <- colnames(get_gwasPval(object))

  		for (i_phen in seq_len(nGWAS)){
    		Prob_e_ijt = get_summary(object)$Sum_E_ijt[i_phen,i_phen,] / length(get_fit(object)$loglik)
        fdrmat[,i_phen] <- 1 - Prob_e_ijt
  		}
		} else if ( !is.null(i) & !is.null(j) ) {
		  #message( "Info: Local FDR vector for specified i & j pair is returned." )

		  fdrmat <- rep( NA, nBin )

		  Prob_e_ijt = get_summary(object)$Sum_E_ijt[i,j,] / length(get_fit(object)$loglik)
      fdrmat <- 1 - Prob_e_ijt
		} else {
		  stop( "Both of i and j should be either NULL or numeric!" )
		}

    return(fdrmat)
  }
)

# parameter estimates

setMethod(
    f="estimates",
    signature="GGPA2",
    definition=function( object, ... ) {
        # return parameter estimates

		return(get_summary(object))
  }
)
