# Independent matrices for calibration and paired same-target comparisons.
revision_matrix <- function(job, design, template) {
    assign(".Random.seed",design$streams[[job$id]],envir=.GlobalEnv)
    sample_sizes <- as.integer(strsplit(job$size,"_",fixed=TRUE)[[1]])
    group <- rep(0:1,sample_sizes)
    G <- 480L
    row <- template[sample.int(nrow(template),G,replace=TRUE),]
    label <- rep("observed_null",G)
    p0 <- p1 <- row$probability
    c0 <- c1 <- row$cv2
    mu0 <- mu1 <- row$mean_positive*p0
    if (job$panel=="mixed") {
        label <- rep(c("observed_null","observed_null","observed_null","mean_only","spread_only","reinforcing","compensated","detection_only"),each=G/8)
        i <- label=="mean_only"; mu1[i] <- 2*mu0[i]
        i <- label=="spread_only"; p0[i] <- p1[i] <- 0.85; c0[i] <- 0.5; c1[i] <- 2
        i <- label=="reinforcing"; p0[i] <- 0.8; p1[i] <- 0.4; c0[i] <- 0.5; c1[i] <- 2
        i <- label=="compensated"; p0[i] <- 0.65; p1[i] <- 0.85; c0[i] <- 0.5; c1[i] <- (1+c0[i])*p1[i]/p0[i]-1
        i <- label=="detection_only"; p0[i] <- 0.85; p1[i] <- 0.50; c0[i] <- c1[i] <- 0.5
    }
    tech <- job$panel %in% c("depth","equal_depth")
    if (tech) {
        label <- rep("same_latent_dispersion",G)
        mu0 <- mu1 <- 2^seq(-1,4,length.out=G)
        order <- sample.int(G); mu0 <- mu1 <- mu0[order]
        phi <- rep(c(0.25,0.5,1),length.out=G)
        s0 <- 1; s1 <- if(job$panel=="depth") 0.25 else 1
        exposure <- ifelse(group==0,s0,s1)
        latent <- matrix(rgamma(G*length(group),shape=1/phi,scale=mu0*phi),nrow=G)
        counts <- matrix(rpois(length(latent),sweep(latent,2,exposure,"*")),nrow=G)
        x <- sweep(counts,2,exposure,"/")
        cv0 <- phi+1/(s0*mu0); cv1 <- phi+1/(s1*mu1)
        p0 <- 1-(1+s0*mu0*phi)^(-1/phi); p1 <- 1-(1+s1*mu1*phi)^(-1/phi)
    } else {
        positive_mean <- cbind(mu0/p0,mu1/p1)[,group+1L,drop=FALSE]
        positive_cv <- cbind(c0,c1)[,group+1L,drop=FALSE]
        pr <- cbind(p0,p1)[,group+1L,drop=FALSE]
        # Both positive families have the specified arithmetic means and CV squared.
        logvar <- log1p(positive_cv)
        ln <- matrix(rlnorm(length(pr),log(positive_mean)-logvar/2,sqrt(logvar)),nrow=G)
        ga <- matrix(rgamma(length(pr),shape=1/positive_cv,scale=positive_mean*positive_cv),nrow=G)
        family_gamma <- rep(c(TRUE,FALSE),length.out=G)
        x <- ln; x[family_gamma,] <- ga[family_gamma,]
        x <- x*(matrix(runif(length(pr)),nrow=G)<pr)
        cv0 <- (1+c0)/p0-1; cv1 <- (1+c1)/p1-1
    }
    rownames(x) <- sprintf("g%03d",seq_len(G)); colnames(x) <- sprintf("s%03d",seq_along(group))
    truth <- data.frame(id=rownames(x),state=label,mean_effect=log2(mu1/mu0),dv_effect=log2(cv1/cv0),
        positive0=p0,positive1=p1,count0=rowSums(x[,group==0,drop=FALSE]>0),count1=rowSums(x[,group==1,drop=FALSE]>0))
    list(x=x,group=group,truth=truth)
}

# A two-channel baseline for the same arithmetic mean and total log-CV-squared targets.
moment_workflow <- function(x,groups) {
    samples <- lapply(0:1,function(g) {
        y <- x[,groups==g,drop=FALSE]; n <- ncol(y); m <- rowMeans(y); v <- matrixStats::rowVars(y)
        delta <- y-m
        influence <- (delta^2-v)/v-2*delta/m
        list(m=m,v=v,cv=v/m^2,n=n,npos=rowSums(y>0),
             scale=matrixStats::rowVars(influence)/n,meanvar=v/n)
    })
    a <- samples[[1]]; b <- samples[[2]]
    welch <- function(effect,s0,s1) {
        se <- sqrt(s0+s1)
        df <- (s0+s1)^2/(s0^2/(a$n-1)+s1^2/(b$n-1))
        2*stats::pt(-abs(effect)/se,df)
    }
    pmean <- welch(b$m-a$m,a$meanvar,b$meanvar)
    e <- log(b$cv/a$cv); pdv <- welch(e,a$scale,b$scale)
    # Match the SGCB minimum support and retain the complete adjustment family.
    pmean[a$npos<1 | b$npos<1 | !is.finite(pmean)] <- NA_real_
    ok <- a$npos>=4 & b$npos>=4 & is.finite(pdv) & is.finite(e)
    pdv[!ok] <- NA_real_
    bounds <- stats::qt(.975,(a$scale+b$scale)^2/(a$scale^2/(a$n-1)+b$scale^2/(b$n-1)))*sqrt(a$scale+b$scale)/log(2)
    lo <- e/log(2)-bounds; hi <- e/log(2)+bounds; lo[!ok] <- hi[!ok] <- NA_real_
    data.frame(mean_p=pmean,mean_q=stats::p.adjust(pmean,"BH",n=nrow(x)),dv_p=pdv,
        dv_q=stats::p.adjust(pdv,"BH",n=nrow(x)),estimate=e/log(2),lower=lo,upper=hi,
        mean_available=is.finite(pmean),dv_available=ok)
}
