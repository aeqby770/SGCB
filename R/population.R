# Native population CV2 inference shared with the fixed-data analyses.
#' @import data.table
NULL

utils::globalVariables(c(".","Q","epsilon","left","log_sigma","moment_increment",
    "relative_shift","right","sharing_key","tail_mass","weight"))

.pairPopulationGG <- function(state, threads, control, previous=NULL) {
    atoms <- data.table::copy(state$atoms)
    atoms[,sharing_key:=.GRP,by=.(epsilon,relative_shift)]
    atoms[,tail_mass:=sum(weight),by=sharing_key]
    pairs <- data.table::CJ(left=seq_len(nrow(atoms)),right=seq_len(nrow(atoms)))
    initial <- switch(as.character(is.null(previous)),
        "TRUE"=atoms$weight[pairs$left]*atoms$weight[pairs$right]*
            (1-state$rho+state$rho*(atoms$sharing_key[pairs$left]==
            atoms$sharing_key[pairs$right])/atoms$tail_mass[pairs$left]),
        "FALSE"={
            stopifnot(identical(unname(as.matrix(atoms[seq_len(nrow(previous$atoms)),
                .(log_sigma,Q,epsilon,relative_shift)])),
                unname(as.matrix(previous$atoms[,.(log_sigma,Q,epsilon,relative_shift)]))))
            continued <- numeric(nrow(pairs))
            continued[(previous$pairs$left-1L)*nrow(atoms)+previous$pairs$right] <-
                previous$pairs$weight
            continued
        })
    ng <- length(state$positive)/2L
    usable <- state$usable[seq_len(ng)] & state$usable[ng+seq_len(ng)]
    elapsed <- system.time(fit <- fit_paired_population_cpp(state$log_likelihood,
        usable,initial,mixsqp::mixsqp,threads,max(control$maximum_columns,nrow(pairs)),control$pair_tolerance,
        column_batch=control$pair_column_batch,sparse_output=TRUE,block_size=32L,symmetric=FALSE))
    stopifnot(fit$converged)
    pairs <- data.table::data.table(left=fit$left,right=fit$right,weight=fit$weights)
    pairs <- data.table::rbindlist(list(
        pairs[weight>0,.(left,right,weight=weight/2)],
        pairs[weight>0,.(left=right,right=left,weight=weight/2)]))[
        ,.(weight=sum(weight)),by=.(left,right)]
    state$atoms <- atoms
    state$pairs <- pairs
    state$pair_fit <- fit
    state$pair_seconds <- elapsed[["elapsed"]]
    state
}

.refinePopulationGG <- function(state, x, group, centers, factor, threads, control) {
    old <- data.table::copy(state$atoms)
    offsets <- seq(-1,1,by=1/factor)
    shape_offsets <- switch(control$positive_family,GG=offsets,Gamma=0,lognormal=0)
    local <- data.table::CJ(base_index=seq_len(nrow(centers)),
        scale_offset=offsets,shape_offset=shape_offsets)
    fine <- unique(data.table::data.table(
        log_sigma=centers$log_sigma[local$base_index]+local$scale_offset*state$scale_step,
        Q=centers$Q[local$base_index]+local$shape_offset*state$shape_step))
    fine[,Q:=switch(control$positive_family,GG=Q,Gamma=exp(log_sigma),lognormal=0)]
    fine <- switch(control$positive_family,
        GG=unique(rbind(fine,unique(fine[,.(log_sigma,Q=0)]),
            unique(fine[,.(log_sigma,Q=exp(log_sigma))]))),unique(fine))
    fine <- fine[1+2*exp(log_sigma)*Q>0]
    fine <- fine[!old[epsilon==0,.(log_sigma,Q)],on=.(log_sigma,Q)]
    fine[,c("epsilon","relative_shift","shift","moment_increment"):=0]
    atoms <- data.table::rbindlist(list(
        old[,.(log_sigma,Q,epsilon,relative_shift,shift,moment_increment)],fine))
    temporary <- tempfile("sgcb-refinement-")
    connection <- file(temporary,open="w+b")
    on.exit({close(connection);file.remove(temporary)},add=TRUE)
    grid_time <- system.time(dictionary <- streamed_refinement_dictionary_cpp(
        x,group,fine$log_sigma,fine$Q,state$log_likelihood,state$usable,
        connection,threads,control$pooled_tolerance))
    initial <- c(old$weight,rep(0,nrow(fine)))
    prior_time <- system.time(population <- fit_sparse_population_dictionary_cpp(
        matrix(numeric(),0L,0L),mixsqp::mixsqp,threads,max(control$maximum_columns,nrow(atoms)),
        control$pooled_tolerance,initial,FALSE,state$usable,TRUE,
        prepared_dictionary=dictionary$dictionary))
    stopifnot(population$converged)
    atoms[,weight:=population$weights]
    active <- atoms$weight>0 | seq_len(nrow(atoms))<=nrow(old) |
        (isTRUE(control$preserve_nested) & atoms$epsilon==0 & (atoms$Q==0 | abs(atoms$Q-exp(atoms$log_sigma))<1e-10))
    retained <- which(active)
    fine_retained <- retained[retained>nrow(old)]-nrow(old)
    retained_time <- system.time(likelihood <- joint_population_likelihood_cpp(
        x,group,fine$log_sigma[fine_retained],fine$Q[fine_retained],
        threads,state$log_likelihood))
    list(atoms=atoms[active],log_likelihood=likelihood$log_likelihood,
        positive=state$positive,usable=state$usable,totals=state$totals,
        quadrature_error=state$quadrature_error,pooled_fit=population,
        scale_step=state$scale_step/factor,shape_step=state$shape_step/factor,
        grid_seconds=grid_time[["elapsed"]]+retained_time[["elapsed"]],prior_seconds=prior_time[["elapsed"]],
        temporary_bytes=dictionary$temporary_bytes)
}

.populationGG <- function(x, group, threads, alpha=.05, family="GG",
    truth=rep(NA_real_,nrow(x)),kernel_cache=NULL) {
    control <- list(positive_family=family,pair_prior="joint_positive_frequency_exchange_symmetrized",
        shape_grid=seq(-6,6,by=.125),
        scale_step=.125,scale_margin=log(2),
        epsilon_grid=stats::plogis(seq(-7,0,by=1)),shift_grid=seq(1,8,by=1),
        refinement=2L,maximum_columns=1000L,cdf_step=.025,probability_tolerance=1e-8,
        structural_mass=.5,
        quadrature_tolerance=1e-5,pooled_tolerance=1e-6,pair_tolerance=1e-6,
        pair_column_batch=256L,likelihood_block_rows=256L)
    start_time <- proc.time()[[3]]
    positive_log <- log(x)
    positive_log[x==0] <- NA_real_
    positive_counts <- cbind(rowSums(x[,group==0,drop=FALSE]>0),
        rowSums(x[,group==1,drop=FALSE]>0))
    control$preserve_nested <- any(x==0)
    # No two-group CV2 posterior is available without four positives per group.
    switch(as.character(any(rowSums(positive_counts>=4L)==2L)),
        "FALSE"={
            missing <- rep(NA_real_,nrow(x))
            genes <- data.frame(effect=missing,se=missing,cdf_zero_lower=missing,
                cdf_zero_upper=missing,cdf_truth_lower=missing,cdf_truth_upper=missing,
                available=0,undefined_mass=missing,quadrature_tail_mass=missing,
                p=missing,covered=missing,Q0_mean=missing,Q1_mean=missing,
                log_sigma0_mean=missing,log_sigma1_mean=missing,
                shared_shape_posterior=missing,shared_shape_bayes_factor=missing,
                frequency_variance=missing,positive_model_variance=missing,
                log_cv2_0_mean=missing,log_cv2_1_mean=missing,
                lower=missing,upper=missing,width=missing,interval_covered=missing)
            return(list(genes=genes,score_type="unavailable_DV",covariance="unavailable",atoms=data.table::data.table(),
                pairs=data.table::data.table(),control=control,
                timing=c(total=proc.time()[[3]]-start_time),
                quadrature_error=rep(NA_real_,2L*nrow(x)),
                convergence=data.frame(stage="no_usable_gene_pairs",kkt_gap=NA_real_)))
        })
    log_scale <- .5*log(cbind(matrixStats::rowVars(positive_log[,group==0,drop=FALSE],na.rm=TRUE),
        matrixStats::rowVars(positive_log[,group==1,drop=FALSE],na.rm=TRUE)))
    scale_range <- range(log_scale[positive_counts>=4],finite=TRUE)
    scale_grid <- seq(floor((scale_range[1]-control$scale_margin)/control$scale_step)*control$scale_step,
        ceiling((scale_range[2]+control$scale_margin)/control$scale_step)*control$scale_step,
        by=control$scale_step)
    atoms <- switch(family,
        GG=unique(rbind(data.table::CJ(log_sigma=scale_grid,Q=control$shape_grid),
            data.table::data.table(log_sigma=scale_grid,Q=exp(scale_grid)))),
        Gamma=data.table::data.table(log_sigma=scale_grid,Q=exp(scale_grid)),
        lognormal=data.table::data.table(log_sigma=scale_grid,Q=0))
    atoms <- atoms[1+2*exp(log_sigma)*Q>0]
    seed_index <- switch(family,GG=which(abs(atoms$Q-exp(atoms$log_sigma))<1e-10 |
        atoms$Q==0),1L)
    base_grid_time <- system.time(likelihood <- joint_population_likelihood_cpp(
        x,group,atoms$log_sigma[seed_index],atoms$Q[seed_index],threads))
    extra_seed <- logical(nrow(atoms))
    seed_time <- system.time(switch(family,GG={
        gamma_index <- which(abs(atoms$Q-exp(atoms$log_sigma))<1e-10)
        lognormal_index <- which(atoms$Q==0)
        gamma_seed <- fit_sparse_population_dictionary_cpp(
            likelihood$log_likelihood[,match(gamma_index,seed_index),drop=FALSE],mixsqp::mixsqp,threads,
            max(control$maximum_columns,length(gamma_index)),control$pooled_tolerance,
            NULL,FALSE,likelihood$usable)
        lognormal_seed <- fit_sparse_population_dictionary_cpp(
            likelihood$log_likelihood[,match(lognormal_index,seed_index),drop=FALSE],mixsqp::mixsqp,threads,
            max(control$maximum_columns,length(lognormal_index)),control$pooled_tolerance,
            NULL,FALSE,likelihood$usable)
        stopifnot(gamma_seed$converged,lognormal_seed$converged)
        extra_seed[gamma_index[gamma_seed$weights>0]] <- TRUE
        extra_seed[lognormal_index[lognormal_seed$weights>0]] <- TRUE
    }))
    control$nested_seed_support <- switch(family,GG=list(gamma=sum(gamma_seed$weights>0),
        lognormal=sum(lognormal_seed$weights>0),gamma_KKT=gamma_seed$kkt_gap,
        lognormal_KKT=lognormal_seed$kkt_gap,seconds=seed_time[["elapsed"]]))
    base_temporary <- tempfile("sgcb-base-")
    base_connection <- file(base_temporary,open="w+b")
    on.exit({close(base_connection);file.remove(base_temporary)},add=TRUE)
    base_stream_time <- system.time(base_dictionary <- streamed_refinement_dictionary_cpp(
        x,group,atoms$log_sigma,atoms$Q,matrix(numeric(),2L*nrow(x),0L),
        likelihood$usable,base_connection,threads,control$pooled_tolerance))
    base_grid_time[["elapsed"]] <- base_grid_time[["elapsed"]]+base_stream_time[["elapsed"]]
    base_prior_time <- system.time(population <- fit_sparse_population_dictionary_cpp(
        matrix(numeric(),0L,0L),mixsqp::mixsqp,threads,
        max(control$maximum_columns,nrow(atoms)),control$pooled_tolerance,
        NULL,FALSE,likelihood$usable,TRUE,prepared_dictionary=base_dictionary$dictionary))
    rm(base_dictionary)
    close(base_connection)
    invisible(file.remove(base_temporary))
    on.exit(NULL)
    stopifnot(population$converged)
    atoms[,weight:=population$weights]
    base_atoms <- atoms[weight>0 | extra_seed]
    base_kkt <- population$kkt_gap
    rm(likelihood)
    expansion <- data.table::CJ(base_index=seq_len(nrow(base_atoms)),
        epsilon=control$epsilon_grid,relative_shift=control$shift_grid)
    atoms <- data.table::rbindlist(list(base_atoms[,.(log_sigma,Q,epsilon=0,relative_shift=0)],
        cbind(base_atoms[expansion$base_index,.(log_sigma,Q)],expansion[,.(epsilon,relative_shift)])))
    atoms[,shift:=relative_shift*exp(log_sigma)]
    atoms[,moment_increment:=log1p(epsilon*expm1(2*shift))-2*log1p(epsilon*expm1(shift))]
    prior_cache <- switch(as.character(is.null(kernel_cache)),
        "TRUE"=list(loglike=NULL,parameters=NULL,reversed=FALSE,
            accuracy=c(control$quadrature_tolerance,control$probability_tolerance)),
        "FALSE"={
            list(loglike=kernel_cache$loglike,parameters=kernel_cache$parameters,
                reversed=kernel_cache$first_group!=group[1L],accuracy=kernel_cache$tolerance)
        })
    mixture_temporary <- tempfile("sgcb-mixture-")
    mixture_connection <- file(mixture_temporary,open="w+b")
    on.exit({close(mixture_connection);file.remove(mixture_temporary)},add=TRUE)
    message("DV mixture: ",nrow(atoms)," states, ",2L*nrow(x)," rows, streamed")
    mixture_grid_time <- system.time(likelihood <- population_mixture_likelihood_cpp(
        x,group,atoms$log_sigma,atoms$Q,atoms$epsilon,atoms$relative_shift,threads,
        control$quadrature_tolerance,control$probability_tolerance,
        previous_loglike=prior_cache$loglike,previous_parameters=prior_cache$parameters,
        previous_reversed=prior_cache$reversed,previous_tolerance=prior_cache$accuracy[1],
        previous_tail_probability=prior_cache$accuracy[2],sparse_output=TRUE,
        connection=mixture_connection,sparse_tolerance=control$pooled_tolerance,
        sparse_block_rows=control$likelihood_block_rows))
    control$kernel_reuse <- list(columns=likelihood$matched_columns,cells=likelihood$reused_cells)
    rm(prior_cache)
    switch(as.character(is.null(kernel_cache)),"FALSE"={
        kernel_cache$loglike <- NULL
        kernel_cache$parameters <- NULL
    })
    invisible(gc(FALSE))
    initial <- c(base_atoms$weight,rep(0,nrow(expansion)))
    mixture_prior_time <- system.time(population <- fit_sparse_population_dictionary_cpp(
        matrix(numeric(),0L,0L),mixsqp::mixsqp,threads,
        max(control$maximum_columns,nrow(atoms)),control$pooled_tolerance,
        initial,FALSE,likelihood$usable,TRUE,prepared_dictionary=likelihood$dictionary))
    stopifnot(population$converged)
    atoms[,weight:=population$weights]
    active <- atoms$weight>0 | (isTRUE(control$preserve_nested) & atoms$epsilon==0 &
        (atoms$Q==0 | abs(atoms$Q-exp(atoms$log_sigma))<1e-10))
    retained_time <- system.time(retained_likelihood <- read_streamed_population_columns_cpp(
        mixture_connection,2L*nrow(x),nrow(atoms),which(active),control$likelihood_block_rows))
    mixture_grid_time[["elapsed"]] <- mixture_grid_time[["elapsed"]]+retained_time[["elapsed"]]
    close(mixture_connection)
    invisible(file.remove(mixture_temporary))
    on.exit(NULL)
    state <- list(atoms=atoms[active],log_likelihood=retained_likelihood,
        positive=likelihood$positive,usable=likelihood$usable,totals=tabulate(group+1L,nbins=2L),
        quadrature_error=likelihood$quadrature_error,pooled_fit=population,
        scale_step=control$scale_step,shape_step=control$scale_step)
    rm(likelihood,retained_likelihood)
    invisible(gc(FALSE))
    switch(as.character(!control$preserve_nested),"TRUE"={
        control$structural_mass <- 1
    state$atoms[,sharing_key:=.GRP,by=.(epsilon,relative_shift)]
    initialization_time <- system.time(initialization <- joint_population_cv_cpp(
        state$log_likelihood,state$positive,state$usable,state$atoms$log_sigma,
        state$atoms$Q,state$atoms$weight,state$totals,rep(NA_real_,nrow(x)),
        control$cdf_step,control$probability_tolerance,alpha,threads,TRUE,
        state$atoms$moment_increment,state$atoms$sharing_key,
        structural_mass=control$structural_mass,shared_shape_only=TRUE))
    state$rho <- initialization$shared_shape_probability
    rm(initialization)
    # Refinement depends only on pooled weights; carry the original starting
    # density directly into the final dictionary and optimize the pair prior once.
    coarse <- state
    coarse$atoms[,tail_mass:=sum(weight),by=sharing_key]
    coarse$pairs <- data.table::CJ(left=seq_len(nrow(coarse$atoms)),right=seq_len(nrow(coarse$atoms)))
    coarse$pairs[,weight:=coarse$atoms$weight[left]*coarse$atoms$weight[right]*
        (1-coarse$rho+coarse$rho*(coarse$atoms$sharing_key[left]==coarse$atoms$sharing_key[right])/
            coarse$atoms$tail_mass[left])]
    refined <- .refinePopulationGG(state,x,group,base_atoms,control$refinement[1],threads,control)
    continued <- .pairPopulationGG(refined,threads,control,coarse)
    control$kernel_reuse$refinement_temporary_bytes <- refined$temporary_bytes
    active_pairs <- continued$pairs[weight>0]
    posterior_time <- system.time(posterior <- joint_population_cv_cpp(
        continued$log_likelihood,continued$positive,continued$usable,
        continued$atoms$log_sigma,continued$atoms$Q,continued$atoms$weight,
        continued$totals,truth,control$cdf_step,
        control$probability_tolerance,alpha,threads,FALSE,
        continued$atoms$moment_increment,NULL,active_pairs$left-1L,
        active_pairs$right-1L,active_pairs$weight,structural_mass=control$structural_mass))

        pp <- continued$pairs[weight>0]
        K <- posterior$log_positive_moment_ratio
        phi <- K+log(-expm1(-K))
        base_K <- K-continued$atoms$moment_increment
        base_phi <- base_K+log(-expm1(-base_K))
        comparison <- cbind(phi[pp$right]-phi[pp$left],
            base_phi[pp$right]-phi[pp$left],
            phi[pp$right]-base_phi[pp$left],
            base_phi[pp$right]-base_phi[pp$left])
        intervals_time <- system.time(intervals <- joint_logcv_reference_summary_cpp(
            continued$log_likelihood,continued$positive,continued$usable,phi,
            continued$atoms$epsilon>0 & continued$atoms$relative_shift>0,
            pp$left-1L,pp$right-1L,pp$weight,comparison,
            matrix(0,nrow(pp),4),matrix(Inf,nrow(pp),4),matrix(TRUE,nrow(pp),4),
            truth,alpha,control$probability_tolerance,control$cdf_step,threads,TRUE))
        genes <- as.data.frame(posterior$genes)
        genes$lower <- intervals[,"lower"]
        genes$upper <- intervals[,"upper"]
        genes$width <- intervals[,"width"]
        genes$interval_covered <- intervals[,"interval_covered"]
        genes$covered <- genes$interval_covered
        atoms <- data.table::copy(continued$atoms)
        atoms[,c("positive_index","pi"):=list(seq_len(.N),1)]
        control$reference <- "paired_positive_population_posterior"
        timing <- c(base_grid=base_grid_time[["elapsed"]],base_prior=base_prior_time[["elapsed"]],
            mixture_grid=mixture_grid_time[["elapsed"]],mixture_prior=mixture_prior_time[["elapsed"]],
            pair_initialization=initialization_time[["elapsed"]],
            refined_grid=refined$grid_seconds,refined_prior=refined$prior_seconds,
            joint_pair=continued$pair_seconds,posterior=posterior_time[["elapsed"]],
            intervals=intervals_time[["elapsed"]],total=proc.time()[[3]]-start_time)
        return(list(genes=genes,atoms=atoms,positive_atoms=continued$atoms,pairs=pp,
            score_type="two_sided_positive_population_tail_bound",
            covariance="paired_population_posterior",control=control,timing=timing,
            quadrature_error=continued$quadrature_error,
            convergence=data.frame(stage=c("base","mixture","refined","joint_pair"),
                kkt_gap=c(base_kkt,state$pooled_fit$kkt_gap,refined$pooled_fit$kkt_gap,
                    continued$pair_fit$kkt_gap))))
    })
    refined <- .refinePopulationGG(state,x,group,base_atoms,control$refinement[1],threads,control)
    control$kernel_reuse$refinement_temporary_bytes <- refined$temporary_bytes
    ng <- nrow(x)
    pa <- data.table::copy(refined$atoms)
    ll <- refined$log_likelihood
    ll[!refined$usable,] <- 0
    positive <- refined$positive
    totals <- as.integer(refined$totals)
    pi_grid <- seq.int(0L,max(totals))/max(totals)
    np <- length(pi_grid)
    nt <- nrow(pa)
    log_frequency <- matrix(stats::dbinom(rep(positive,np),
        size=rep(rep(totals,each=ng),np),prob=rep(pi_grid,each=2L*ng),log=TRUE),2L*ng)
    joint_grid_time <- system.time(dictionary <- product_population_dictionary_cpp(
        ll,log_frequency,threads,control$pair_tolerance))
    joint_prior_time <- system.time(joint_prior <- fit_sparse_population_dictionary_cpp(
        matrix(numeric(),0L,0L),mixsqp::mixsqp,threads,max(control$maximum_columns,nt*np),
        control$pair_tolerance,NULL,FALSE,column_batch=control$pair_column_batch,
        prepared_dictionary=dictionary$dictionary))
    stopifnot(joint_prior$converged)
    grid <- data.table::as.data.table(expand.grid(theta=seq_len(nt),pi=seq_len(np)))
    js <- grid[joint_prior$weights>0]
    q_safe <- data.table::fifelse(pa$Q==0,1,pa$Q)
    log_variance <- exp(2*pa$log_sigma)*
        data.table::fifelse(pa$Q==0,1,trigamma(1/q_safe^2)/q_safe^2)+
        pa$epsilon*(1-pa$epsilon)*pa$shift^2
    gamma_index <- which(pa$epsilon==0 & abs(pa$Q-exp(pa$log_sigma))<1e-10)
    normal_index <- which(pa$epsilon==0 & pa$Q==0)
    closure <- switch(family,GG={
        gamma_map <- gamma_index[max.col(-abs(outer(log(log_variance),
            log(log_variance[gamma_index]),"-")),ties.method="first")]
        normal_map <- normal_index[max.col(-abs(outer(log(log_variance),
            log(log_variance[normal_index]),"-")),ties.method="first")]
        data.table::rbindlist(list(data.table::data.table(theta=gamma_map[js$theta],pi=js$pi),
            data.table::data.table(theta=normal_map[js$theta],pi=js$pi)))
    },Gamma={
        gamma_map <- gamma_index[max.col(-abs(outer(log(log_variance),
            log(log_variance[gamma_index]),"-")),ties.method="first")]
        data.table::data.table(theta=gamma_map[js$theta],pi=js$pi)
    },lognormal={
        normal_map <- normal_index[max.col(-abs(outer(log(log_variance),
            log(log_variance[normal_index]),"-")),ties.method="first")]
        data.table::data.table(theta=normal_map[js$theta],pi=js$pi)
    })
    js <- unique(data.table::rbindlist(list(js,closure)))
    marginal <- joint_prior$weights[(js$pi-1L)*nt+js$theta]
    marginal <- marginal/sum(marginal)
    joint_ll <- ll[,js$theta,drop=FALSE]+log_frequency[,js$pi,drop=FALSE]
    joint_pair_time <- system.time(joint_fit <- fit_paired_population_cpp(
        joint_ll,rep(TRUE,ng),numeric(),mixsqp::mixsqp,threads,
        max(control$maximum_columns,nrow(js)^2),control$pair_tolerance,
        column_batch=control$pair_column_batch,marginal_start=marginal,
        sparse_output=TRUE,symmetric=FALSE))
    stopifnot(joint_fit$converged)
    pairs <- data.table::rbindlist(list(
        data.table::data.table(left=joint_fit$left,right=joint_fit$right,weight=joint_fit$weights/2),
        data.table::data.table(left=joint_fit$right,right=joint_fit$left,weight=joint_fit$weights/2)))[
        ,.(weight=sum(weight)),by=.(left,right)][weight>0]
    state_used <- sort(unique(c(pairs$left,pairs$right)))
    js <- data.table::copy(js[state_used])
    pairs[,c("left","right"):=list(match(left,state_used),match(right,state_used))]
    theta_used <- sort(unique(js$theta))
    pa <- pa[theta_used]
    js[,theta:=match(theta,theta_used)]
    joint_ll <- joint_ll[,state_used,drop=FALSE]
    rm(ll,log_frequency,dictionary)
    invisible(gc(FALSE))
    membership_time <- system.time(membership <- joint_state_responsibilities_cpp(
        joint_ll,js$theta-1L,js$pi-1L,nrow(pa),np,
        pairs$left-1L,pairs$right-1L,pairs$weight,threads))
    location_time <- system.time(location <- population_mixture_mean_likelihood_cpp(
        x,group,pa$log_sigma,pa$Q,pa$epsilon,pa$relative_shift,threads,
        control$quadrature_tolerance,control$probability_tolerance))
    information_time <- system.time(information <- population_joint_profile_information_cpp(
        log(x),group,location$log_mean,membership$positive,pa$log_sigma,pa$Q,
        pa$epsilon,pa$relative_shift,threads,switch(family,GG=0L,Gamma=1L,lognormal=2L)))
    rm(location)
    positive_influence <- attr(information,"sample_influence_log_positive_cv2")
    identified_theta <- is.finite(information[,"log_positive_cv2_sd"])
    base_influence <- attr(information,"sample_influence_log_base_cv2")
    identified_base <- is.finite(information[,"log_base_cv2_sd"])
    S <- information[,"log_positive_cv2"]
    K <- pmax(S,0)+log1p(exp(-abs(S)))
    rp0 <- membership$frequency[seq_len(ng),,drop=FALSE]
    rp1 <- membership$frequency[ng+seq_len(ng),,drop=FALSE]
    mass0 <- colSums(rp0)
    mass1 <- colSums(rp1)
    exposure <- totals[1]*mass0+totals[2]*mass1
    effective_samples <- exposure^2/pmax(totals[1]*mass0^2+totals[2]*mass1^2,.Machine$double.xmin)
    correction <- sqrt(effective_samples/pmax(effective_samples-1,1))
    score_pi <- matrix(0,np,ncol(x))
    score_pi[,group==0L] <- sweep(crossprod(rp0,1*(x[,group==0L,drop=FALSE]>0)),1,
        drop(crossprod(rp0,positive[seq_len(ng)]/totals[1])),"-")
    score_pi[,group==1L] <- sweep(crossprod(rp1,1*(x[,group==1L,drop=FALSE]>0)),1,
        drop(crossprod(rp1,positive[ng+seq_len(ng)]/totals[2])),"-")
    frequency_influence <- sweep(score_pi,1,
        correction/data.table::fifelse(pi_grid*exposure>0,pi_grid*exposure,1),"*")
    positive_influence[!is.finite(positive_influence)] <- 0
    base_influence[!is.finite(base_influence)] <- 0
    K_influence <- sweep(positive_influence,1,plogis(S),"*")
    base_S <- information[,"log_base_cv2"]
    base_K <- pmax(base_S,0)+log1p(exp(-abs(base_S)))
    base_K_influence <- sweep(base_influence,1,plogis(base_S),"*")
    T <- K[js$theta]-log(pi_grid[js$pi])
    base_T <- base_K[js$theta]-log(pi_grid[js$pi])
    phi <- T+log(-expm1(-T))
    base_phi <- base_T+log(-expm1(-base_T))
    state_influence <- sweep(K_influence[js$theta,,drop=FALSE]-
        frequency_influence[js$pi,,drop=FALSE],1,1/(-expm1(-T)),"*")
    base_state_influence <- sweep(base_K_influence[js$theta,,drop=FALSE]-
        frequency_influence[js$pi,,drop=FALSE],1,1/(-expm1(-base_T)),"*")
    pair_influence <- rbind(
        state_influence[pairs$right,,drop=FALSE]-state_influence[pairs$left,,drop=FALSE],
        base_state_influence[pairs$right,,drop=FALSE]-state_influence[pairs$left,,drop=FALSE],
        state_influence[pairs$right,,drop=FALSE]-base_state_influence[pairs$left,,drop=FALSE],
        base_state_influence[pairs$right,,drop=FALSE]-base_state_influence[pairs$left,,drop=FALSE])
    influence_scale <- pmax(matrixStats::rowMaxs(abs(pair_influence)),.Machine$double.xmin)
    scaled_influence <- pair_influence/influence_scale
    variance <- rowSums(scaled_influence^2)
    variance_of_variance <- totals[1]*matrixStats::rowVars(scaled_influence[,group==0L,drop=FALSE]^2)+
        totals[2]*matrixStats::rowVars(scaled_influence[,group==1L,drop=FALSE]^2)
    contrast_sd <- matrix(influence_scale*sqrt(variance),nrow(pairs),4L)
    contrast_df <- matrix(data.table::fifelse(variance_of_variance>0,
        2*variance^2/variance_of_variance,Inf),nrow(pairs),4L)
    contrast_mean <- cbind(phi[pairs$right]-phi[pairs$left],
        base_phi[pairs$right]-phi[pairs$left],
        phi[pairs$right]-base_phi[pairs$left],
        base_phi[pairs$right]-base_phi[pairs$left])
    pi_valid <- pi_grid[js$pi[pairs$left]]>0 & pi_grid[js$pi[pairs$right]]>0
    contrast_identified <- cbind(
        (pairs$left==pairs$right | (identified_theta[js$theta[pairs$left]] &
            identified_theta[js$theta[pairs$right]])) & pi_valid,
        identified_theta[js$theta[pairs$left]] & identified_base[js$theta[pairs$right]] & pi_valid,
        identified_base[js$theta[pairs$left]] & identified_theta[js$theta[pairs$right]] & pi_valid,
        (pairs$left==pairs$right | (identified_base[js$theta[pairs$left]] &
            identified_base[js$theta[pairs$right]])) & pi_valid)
    mixed <- pa$epsilon[js$theta]>0 & pa$relative_shift[js$theta]>0
    reference_time <- system.time(posterior <- joint_logcv_reference_summary_cpp(
        joint_ll,positive,refined$usable,phi,mixed,pairs$left-1L,pairs$right-1L,pairs$weight,
        contrast_mean,contrast_sd,contrast_df,contrast_identified,truth,alpha,
        control$probability_tolerance,control$cdf_step,threads,TRUE))
    pairs[,c("reference_mean","reference_sd","reference_df","reference_identified"):=
        list(contrast_mean[,1],contrast_sd[,1],contrast_df[,1],contrast_identified[,1])]
    state_atoms <- data.table::copy(pa[js$theta])
    state_atoms[,c("positive_index","pi","weight"):=list(js$theta,pi_grid[js$pi],0)]
    pair_mass <- data.table::data.table(index=c(pairs$left,pairs$right),
        value=c(pairs$weight,pairs$weight)/2)[,.(weight=sum(value)),by=index]
    state_atoms[pair_mass$index,weight:=pair_mass$weight]
    genes <- as.data.frame(posterior)
    genes$covered <- genes$interval_covered
    control$frequency_grid <- pi_grid
    control$structural_mass <- sum(joint_prior$weights[grid$pi==np])
    control$joint_candidate_states <- nrow(grid)
    control$occupied_joint_states <- nrow(js)
    control$reference <- "sample_cluster_studentized_logcv2_mixture_bounds"
    control$reference_df <- "variance_contribution_moment_match"
    timing <- c(base_grid=base_grid_time[["elapsed"]],base_prior=base_prior_time[["elapsed"]],
        mixture_grid=mixture_grid_time[["elapsed"]],mixture_prior=mixture_prior_time[["elapsed"]],
        refined_grid=refined$grid_seconds,refined_prior=refined$prior_seconds,
        joint_grid=joint_grid_time[["elapsed"]],joint_prior=joint_prior_time[["elapsed"]],
        joint_pair=joint_pair_time[["elapsed"]],membership=membership_time[["elapsed"]],
        location=location_time[["elapsed"]],information=information_time[["elapsed"]],
        reference=reference_time[["elapsed"]],total=proc.time()[[3]]-start_time)
    list(genes=genes,atoms=state_atoms,positive_atoms=pa,pairs=pairs,control=control,
        score_type="two_sided_studentized_reference_tail_bound",covariance="sample_cluster_joint_parameter_influence",
        reference=list(mean=contrast_mean,sd=contrast_sd,df=contrast_df,
            identified=contrast_identified,mixed=mixed),
        timing=timing,quadrature_error=refined$quadrature_error,
        convergence=data.frame(stage=c("base","mixture","refined","joint_marginal","joint_pair"),
            kkt_gap=c(base_kkt,state$pooled_fit$kkt_gap,refined$pooled_fit$kkt_gap,
                joint_prior$kkt_gap,joint_fit$kkt_gap)))
}

# Arithmetic means use a finite-first-moment GG nuisance fit and a contrast
# population. The CV population is fitted separately on finite second moments;
# exact mixed-GG likelihood values are reused within the same call.
.populationMean <- function(x, group, threads,kernel_cache=NULL) {
    direction <- ifelse(group[1L]==0L,1,-1)
    group <- as.integer(group!=group[1L])
    control <- list(shape_grid=seq(-6,6,by=.125),scale_step=.125,
        scale_margin=log(2),epsilon_grid=numeric(),
        shift_grid=seq(1,8,by=1),maximum_columns=2000L,
        pooled_tolerance=1e-6,pair_tolerance=1e-6,weight_tolerance=1e-6,
        quadrature_tolerance=1e-9,probability_tolerance=1e-10,
        effect_step=.01,pair_block_size=32L,pair_column_batch=64L,
        mean_step=.125,mean_bandwidth=.125,slab_mean=0,slab_sd=1,
        slab_scales=2^seq(-3,2),
        tail_probability=1e-12,mean_likelihood="four_cumulant_normal_mixture",
        residual_likelihood="conditional_first_two_moments")
    started <- proc.time()[["elapsed"]]
    positive_log <- log(x)
    positive_log[x==0] <- NA_real_
    log_scale <- .5*log(cbind(
        matrixStats::rowVars(positive_log[,group==0,drop=FALSE],na.rm=TRUE),
        matrixStats::rowVars(positive_log[,group==1,drop=FALSE],na.rm=TRUE)))
    scale_range <- range(log_scale,finite=TRUE)
    scale_grid <- seq(
        floor((scale_range[1]-control$scale_margin)/control$scale_step)*control$scale_step,
        ceiling((scale_range[2]+control$scale_margin)/control$scale_step)*control$scale_step,
        by=control$scale_step)
    atoms <- data.table::CJ(log_sigma=scale_grid,Q=control$shape_grid)
    atoms <- atoms[1+exp(log_sigma)*Q>0]
    atoms[,c("epsilon","relative_shift"):=0]
    base_grid_time <- system.time(likelihood <- population_mixture_likelihood_cpp(
        x,group,atoms$log_sigma,atoms$Q,atoms$epsilon,atoms$relative_shift,threads,
        control$quadrature_tolerance,control$probability_tolerance,
        minimum_positive=1L))
    base_prior_time <- system.time(population <- fit_sparse_population_dictionary_cpp(
        likelihood$log_likelihood,mixsqp::mixsqp,threads,
        max(control$maximum_columns,nrow(atoms)),control$pooled_tolerance,
        NULL,FALSE,likelihood$usable))
    stopifnot(population$converged)
    atoms[,weight:=population$weights]
    base_atoms <- atoms[weight>0]
    base_kkt <- population$kkt_gap
    rm(likelihood)
    expansion <- data.table::CJ(base_index=seq_len(nrow(base_atoms)),
        epsilon=control$epsilon_grid,relative_shift=control$shift_grid)
    atoms <- data.table::rbindlist(list(
        base_atoms[,.(log_sigma,Q,epsilon,relative_shift)],
        cbind(base_atoms[expansion$base_index,.(log_sigma,Q)],
            expansion[,.(epsilon,relative_shift)])))
    mixture_grid_time <- system.time(likelihood <- population_mixture_likelihood_cpp(
        x,group,atoms$log_sigma,atoms$Q,atoms$epsilon,atoms$relative_shift,threads,
        control$quadrature_tolerance,control$probability_tolerance,
        minimum_positive=1L))
    switch(as.character(is.null(kernel_cache)),"FALSE"={
        kernel_cache$loglike <- likelihood$log_likelihood
        kernel_cache$parameters <- as.matrix(atoms[,.(log_sigma,Q,epsilon,relative_shift)])
        kernel_cache$first_group <- group[1L]
        kernel_cache$tolerance <- c(control$quadrature_tolerance,control$probability_tolerance)
    })
    initial <- c(base_atoms$weight,rep(0,nrow(expansion)))
    mixture_prior_time <- system.time(population <- fit_sparse_population_dictionary_cpp(
        likelihood$log_likelihood,mixsqp::mixsqp,threads,
        max(control$maximum_columns,nrow(atoms)),control$pooled_tolerance,
        initial,FALSE,likelihood$usable))
    stopifnot(population$converged)
    atoms[,weight:=population$weights]
    # Estimate empirical nuisance support from all group-wise likelihoods
    # before learning its cross-group dependence.
    active <- atoms$weight>0
    ll <- likelihood$log_likelihood[,active,drop=FALSE]
    atoms <- atoms[active]
    rm(likelihood)
    # No positive observation contributes the constant likelihood one.
    # Such a gene remains in the fitting family and has no mean test.
    ll[is.na(ll)] <- 0
    ng <- nrow(x)
    pair_time <- system.time(paired <- fit_paired_population_cpp(ll,rep(TRUE,ng),
        numeric(),mixsqp::mixsqp,threads,
        as.integer(min(.Machine$integer.max,max(control$maximum_columns,nrow(atoms)^2))),
        control$pair_tolerance,marginal_start=atoms$weight,
        sparse_output=TRUE,block_size=control$pair_block_size,
        column_batch=control$pair_column_batch))
    stopifnot(paired$converged)
    pairs <- data.table::data.table(left=paired$left,right=paired$right,
        weight=paired$weights)[weight>0]
    marginal <- data.table::data.table(index=c(pairs$left,pairs$right),
        weight=rep(pairs$weight,2))[,.(weight=sum(weight)/2),by="index"]
    atoms[,weight:=0]
    atoms[marginal$index,weight:=marginal$weight]
    retained <- which(atoms$weight>0)
    pairs[,c("left","right"):=list(match(left,retained),match(right,retained))]
    ll <- ll[,retained,drop=FALSE]
    atoms <- atoms[retained]
    # Retain fitted states and add the nested pure-lognormal support.
    profile_start_pairs <- data.table::copy(pairs[,.(left,right,weight)])
    pure <- unique(atoms[,.(log_sigma,Q=0,epsilon=0,relative_shift=0,weight=0)])
    atoms <- unique(rbind(atoms,pure),by=c("log_sigma","Q","epsilon","relative_shift"))

    original_pairs <- profile_start_pairs[,.(left,right)]
    gaussian <- which(atoms$Q==0 & atoms$epsilon==0)
    lookup <- gaussian[match(atoms$log_sigma,atoms$log_sigma[gaussian])]
    gaussian_pairs <- data.table::data.table(left=lookup[original_pairs$left],right=lookup[original_pairs$right])
    mean_pairs <- unique(data.table::rbindlist(list(original_pairs,
        original_pairs[,.(left=right,right=left)],
        data.table::data.table(left=seq_len(nrow(atoms)),right=seq_len(nrow(atoms))),
        gaussian_pairs,gaussian_pairs[,.(left=right,right=left)])))

    moment_time <- system.time(location <- population_mixture_mean_likelihood_cpp(
        x,group,atoms$log_sigma,atoms$Q,atoms$epsilon,atoms$relative_shift,threads,
        control$quadrature_tolerance,control$probability_tolerance))
    ll <- location$log_likelihood
    ll[is.na(ll)] <- 0
    totals <- rep(c(sum(group==0L),sum(group==1L)),each=ng)
    positive <- location$positive
    control$frequency_fit <- population_frequency_structure_cpp(positive,
        as.integer(tabulate(group+1L,nbins=2L)))
    control$structural_mass <- control$frequency_fit$mass
    log_continuous <- log1p(-control$structural_mass)+
        lbeta(positive+.5,totals-positive+.5)-lbeta(.5,.5)
    log_spike <- ifelse(positive==totals,log(control$structural_mass),-Inf)
    continuous <- exp(log_continuous-matrixStats::rowLogSumExps(cbind(log_continuous,log_spike)))
    log_frequency <- continuous*(digamma(positive+.5)-digamma(totals+1))
    log_total_mean <- location$log_mean+log_frequency
    frequency_variance <- continuous*(trigamma(positive+.5)-trigamma(totals+1))+
        continuous*(1-continuous)*(digamma(positive+.5)-digamma(totals+1))^2
    conditional_variance <- sweep(location$log_mean_variance,1,frequency_variance,"+")
    common_pairs <- unique(rbind(original_pairs,original_pairs[,.(left=right,right=left)]))
    common <- (log_total_mean[seq_len(ng),common_pairs$left,drop=FALSE]+
        log_total_mean[ng+seq_len(ng),common_pairs$right,drop=FALSE])/2
    limits <- range(common,finite=TRUE)
    A <- seq(floor(limits[1]/control$mean_step)*control$mean_step-control$mean_step,
        ceiling(limits[2]/control$mean_step)*control$mean_step+control$mean_step,
        by=control$mean_step)
    rm(common)
    location <- location["quadrature_error"]
    parent_output <- Sys.getenv("SGCB_MEAN_PARENT_OUTPUT")
    if(nzchar(parent_output)) {
        parent_state <- list(ll=ll,mean=log_total_mean,variance=conditional_variance,
            positive=positive,pairs=mean_pairs,A=A,control=control,
            atoms=atoms,direction=direction,
            profile_start_pairs=profile_start_pairs,gene_id=rownames(x),
            samples=colnames(x),group=group,phase="paired_nuisance_ready")
        saveRDS(parent_state,parent_output,compress=FALSE)
    }
    # Learn one continuous nonzero-effect distribution from the paired
    # location-invariant likelihood, without labels or external cached fits.
    profile_time <- system.time({
        profile_start <- numeric(nrow(atoms)^2)
        profile_start[(profile_start_pairs$left-1L)*nrow(atoms)+profile_start_pairs$right] <-
            profile_start_pairs$weight
        profile_pair <- fit_paired_population_cpp(ll,rep(TRUE,ng),profile_start,
            mixsqp::mixsqp,threads,as.integer(max(control$maximum_columns,nrow(atoms)^2)),
            control$pair_tolerance,sparse_output=TRUE,
            block_size=control$pair_block_size,initialize_support=TRUE,
            column_batch=control$pair_column_batch)
        stopifnot(profile_pair$converged)
        profile_pairs <- data.table::data.table(left=profile_pair$left,
            right=profile_pair$right,weight=profile_pair$weights)[weight>0]
        profile_marginal <- data.table::data.table(
            index=c(profile_pairs$left,profile_pairs$right),
            weight=rep(profile_pairs$weight,2))[,.(weight=sum(weight)/2),by="index"]
        profile_weights <- numeric(nrow(atoms))
        profile_weights[profile_marginal$index] <- profile_marginal$weight
        flat_effect <- paired_mean_location_expectation_cpp(ll,log_total_mean,positive,
            profile_pairs$left,profile_pairs$right,profile_pairs$weight,threads)
        limits <- range(flat_effect,finite=TRUE)
        grid <- seq(floor(limits[1]/control$mean_step),
            ceiling(limits[2]/control$mean_step))*control$mean_step
        slab_means <- c(0,grid[grid!=0])
        # Integration support includes the Gaussian tails; center candidates
        # retain the original mean grid. Use the contrast FFT's own spacing.
        padding <- stats::qnorm(control$tail_probability/2,lower.tail=FALSE)*
            max(control$slab_scales)
        limits <- range(slab_means)+c(-padding,padding)
        grid <- seq(floor(limits[1]/control$effect_step),
            ceiling(limits[2]/control$effect_step))*control$effect_step
        effect_grid <- c(0,grid[grid!=0])
        point <- population_mixture_mean_grid_cpp(x,group,atoms$log_sigma,
            atoms$Q,atoms$epsilon,atoms$relative_shift,profile_weights,ll,
            rep(0,length(effect_grid)),control$effect_step,control$structural_mass,
            control$probability_tolerance,threads,profile_pairs$left,
            profile_pairs$right,profile_pairs$weight,slab_mean=effect_grid)
        point <- point[c("log_likelihood","usable","checks")]
        slab <- mean_effect_profile_cpp(point$log_likelihood,effect_grid,
            control$slab_scales,point$usable,threads,slab_means)
        stopifnot(max(abs(slab$profile$prior_grid_mass-1))<control$probability_tolerance)
        control$slab_mean <- slab$parameters$mean
        control$slab_sd <- slab$parameters$sd
        stopifnot(max(point$checks[,"contrast_mass_error"])<1e-7)
        point_checks <- point$checks
    })
    if(nzchar(parent_output)) {
        parent_state$control <- control
        parent_state$profile_pair_fit <- profile_pair
        parent_state$effect_profile <- slab
        parent_state$phase <- "effect_profile_ready"
        saveRDS(parent_state,parent_output,compress=FALSE)
    }
    rm(point)
    invisible(gc(FALSE))
    joint_kernel_time <- system.time(joint_kernel <- anchored_location_dictionary_cpp(
        ll,log_total_mean,conditional_variance,positive,mean_pairs$left,mean_pairs$right,
        A,control$mean_bandwidth,control$slab_mean,control$slab_sd,threads,
        atoms$log_sigma,atoms$Q,control$tail_probability,
        anchor=ifelse(direction==1,0,-1),symmetric=FALSE))
    joint_fit_time <- system.time(joint_fit <- fit_shared_location_cpp(joint_kernel$dictionary,
        threads,max(control$maximum_columns,ng),control$weight_tolerance,maximum_seconds=0))
    stopifnot(joint_fit$converged)
    posterior_time <- system.time(summary <- anchored_location_summary_cpp(
        joint_kernel$dictionary,joint_fit$weights,.05,threads,control$probability_tolerance))
    density_error <- max(abs(summary[,"density"]/joint_fit$density-1))
    null_error <- max(abs(summary[,"posterior_null"]-joint_fit$posterior_null))
    stopifnot(density_error<1e-8,null_error<1e-8)
    posterior_null <- joint_fit$posterior_null
    available <- positive[seq_len(ng)]>0L & positive[ng+seq_len(ng)]>0L
    score <- ifelse(available,posterior_null,Inf)
    total_null <- sum(posterior_null)
    p <- cumsum(posterior_null[order(score)])[data.table::frank(score,ties.method="max")]/
        ifelse(total_null>0,total_null,1)
    p[!available] <- 1
    moments <- summary
    genes <- data.frame(effect=direction*moments[,"effect"],
        se=sqrt(pmax(0,moments[,"variance"])),p=p,available=as.integer(available),
        posterior_null=posterior_null,
        lower=ifelse(rep(direction==1,nrow(moments)),moments[,"lower"],-moments[,"upper"]),
        upper=ifelse(rep(direction==1,nrow(moments)),moments[,"upper"],-moments[,"lower"]))
    shared_weights <- rowSums(matrix(joint_fit$weights,ncol=2L))
    mean_pairs[,weight:=colSums(matrix(shared_weights,nrow=length(A)))]
    positive_pairs <- mean_pairs[weight>0]
    final_marginal <- data.table::data.table(index=c(positive_pairs$left,positive_pairs$right),
        weight=rep(positive_pairs$weight,2))[,.(weight=sum(weight)/2),by="index"]
    atoms[,weight:=0]
    atoms[final_marginal$index,weight:=final_marginal$weight]
    location_prior <- data.table::data.table(state=which(shared_weights>0),weight=shared_weights[shared_weights>0])
    location_prior[,c("A_index","pair_index"):=list(
        (state-1L)%%length(A)+1L,(state-1L)%/%length(A)+1L)]
    control$mean_likelihood <- "reference_anchored_conditional_mean"
    control$nuisance_structure <- "shared_reference_level_and_GG_shape"
    control$reference_level <- "API control group"
    control$prior_null <- joint_fit$prior_null
    convergence <- data.frame(stage=c("mean_base","mean_mixture","mean_pair","mean_profile_pair","mean_shared_reference"),
        kkt_gap=c(base_kkt,population$kkt_gap,paired$kkt_gap,profile_pair$kkt_gap,joint_fit$kkt_gap))
    mean_pairs[,c("left","right"):=list(ifelse(rep(direction==1,.N),left,right),
        ifelse(rep(direction==1,.N),right,left))]
    slab$parameters$mean <- direction*slab$parameters$mean
    control$slab_mean <- direction*control$slab_mean
    list(genes=genes,atoms=atoms,pairs=mean_pairs,control=control,A_grid=A,
        location_prior=location_prior,
        effect_profile=slab$parameters,quadrature_error=location$quadrature_error,
        integration_checks=c(maximum_contrast_mass_error=max(point_checks[,"contrast_mass_error"]),
            profile_null_weight_KKT=slab$null_weight_kkt_gap,
            posterior_density_identity_error=density_error,posterior_null_identity_error=null_error),
        timing=c(base_grid=base_grid_time[["elapsed"]],base_prior=base_prior_time[["elapsed"]],
            mixture_grid=mixture_grid_time[["elapsed"]],mixture_prior=mixture_prior_time[["elapsed"]],
            paired_prior=pair_time[["elapsed"]],location=moment_time[["elapsed"]],
            effect_profile=profile_time[["elapsed"]],joint_kernel=joint_kernel_time[["elapsed"]],
            joint_prior=joint_fit_time[["elapsed"]],posterior=posterior_time[["elapsed"]],
            total=proc.time()[["elapsed"]]-started),
        convergence=convergence,
        posterior_convergence=data.frame(quantity=c("shared_nuisance","global_null_weight"),
            kkt_gap=c(joint_fit$nuisance_gap,joint_fit$class_gap)),
        rng_state=NULL,endpoint="log2 ratio of unconditional arithmetic means",
        score_type="posterior_null_mass_under_reference_anchored_population",
        support="finite first moment; GG conditional location approximation and positive-frequency uncertainty; shared prior over the API control mean and paired GG shapes; global point-null fraction and Gaussian nonzero effects independent of that baseline; finite learned nuisance catalogue")
}
