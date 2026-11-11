known_depth_inputs <- function(rep_id,design,template) {
    id <- 384+rep_id; job <- design$jobs[id,]; sim <- revision_matrix(job,design,template)
    g <- sim$group; exposure <- ifelse(g==0,1,.25)
    raw <- round(sweep(sim$x,2,exposure,"*")); G <- nrow(raw)
    # The sibling stream supplies new biological alternatives; no streams are selected by outcomes.
    alt_stream <- parallel::nextRNGSubStream(design$streams[[id]])
    assign(".Random.seed",alt_stream,.GlobalEnv)
    mu <- sample(2^seq(-1,4,length.out=G)); phi <- rep(c(.25,.5,1),length.out=G)
    label <- sample(rep(c("null","null","mean_up","mean_down","null","spread"),length.out=G))
    m1 <- mu*ifelse(label=="mean_up",2,ifelse(label=="mean_down",.5,1))
    phi1 <- phi*ifelse(label=="spread",4,1)
    mus <- cbind(mu,m1)[,g+1,drop=FALSE]; phis <- cbind(phi,phi1)[,g+1,drop=FALSE]
    latent <- matrix(rgamma(length(mus),shape=1/phis,scale=mus*phis),nrow=G)
    raw_alt <- matrix(rpois(length(latent),sweep(latent,2,exposure,"*")),nrow=G)
    dimnames(raw_alt) <- dimnames(raw)
    list(id=id,raw=raw,raw_alt=raw_alt,g=g,exposure=exposure,G=G,mu=mu,m1=m1,alt_stream=alt_stream)
}
