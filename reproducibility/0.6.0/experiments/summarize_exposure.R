root <- normalizePath(Sys.getenv("SGCB_R2_ROOT"),mustWork=TRUE); library(data.table)
out <- file.path(root,"results"); paths <- file.path(out,"exposure_final",sprintf("exposure_%02d.tsv",1:24))
stopifnot(all(file.exists(paths)))
x <- rbindlist(lapply(paths,fread)); stopifnot(nrow(x)==96)
x <- rbind(copy(x)[,scope:="all matrices"],copy(x[replicate>8])[,scope:="held-out replicates 9-24"])
s <- x[,.(replicates=.N,any_datasets=sum(any),any_probability=mean(any),any_lower=binom.test(sum(any),.N)$conf.int[1],
    any_upper=binom.test(sum(any),.N)$conf.int[2],mean_FDP=mean(FDP),mean_FDP_se=sd(FDP)/sqrt(.N),
    mean_power=if(all(is.na(power))) NA_real_ else mean(power),power_se=sd(power)/sqrt(.N),false_calls=sum(fp),true_calls=sum(tp)),
    by=.(scope,condition,arm)]
fwrite(s,file.path(out,"exposure_summary.tsv"),sep="\t",na="")
print(s[scope=="all matrices",.(condition,arm,replicates,any_datasets,false_calls,true_calls,mean_FDP,mean_power)])