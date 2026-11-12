root <- normalizePath(Sys.getenv("SGCB_R2_ROOT"),mustWork=TRUE)
library(data.table);library(ggplot2);library(patchwork)
p <- file.path(root,"results"); dest <- file.path(root,"manuscript","figures")
a <- fread(file.path(p,"independent_operating_characteristics.tsv"))
outer <- fread(file.path(p,"outer_refit_summary.tsv"))
a[,size:=factor(size,c("10_10","10_50","50_50"),c("10:10","10:50","50:50"))]
paint <- c(SGCB="#167F83",Moments="#D06C4A")
legend_labels <- c(SGCB="SGCB",Moments="Two-moment workflow")
base <- theme_classic(base_size=11,base_family="Arial")+theme(
    axis.line=element_line(linewidth=.35,colour="#445151"),axis.ticks=element_line(linewidth=.3),
    text=element_text(colour="#263D42"),axis.text=element_text(colour="#364A4D"),
    legend.position="bottom",legend.title=element_blank(),plot.tag=element_text(face="bold",size=15),
    panel.grid.major.y=element_line(colour="#E6EAE7",linewidth=.30),plot.margin=margin(6,9,5,5))
methods <- scale_colour_manual(values=paint,labels=legend_labels)
pd <- position_dodge(width=.35)
operating <- function(ep,which_metric,ylabel) {
   d <- a[panel=="mixed" & state=="all" & method%in%names(paint) & endpoint==ep & metric==which_metric]
   ggplot(d,aes(size,value,colour=method,group=method))+
       geom_linerange(aes(ymin=lower,ymax=upper),position=pd,linewidth=1.05)+
       geom_point(aes(shape=method),position=pd,size=3,stroke=1.0,fill="white")+
       scale_shape_manual(values=c(SGCB=21,Moments=17),labels=legend_labels)+methods+
       scale_y_continuous(labels=scales::label_percent(accuracy=1),limits=c(0,NA),expand=expansion(mult=c(.04,.11)))+
       labs(x="Samples per group",y=ylabel)+base
}
g1 <- operating("mean","FDP","Mean false-discovery proportion")+
    geom_hline(yintercept=.05,linetype=3,colour="#6D746C",linewidth=.45)
g2 <- operating("dv","power","Power for observed CV² changes")
sel <- outer[state=="reinforcing"]
sel[,size:=factor(size,c("10_10","50_50"),c("10:10","50:50"))]
compare <- melt(sel,id.vars="size",measure.vars=c("conditional_coverage","outer_coverage"),
    variable.name="interval",value.name="coverage")
compare[,interval:=factor(interval,levels=c("conditional_coverage","outer_coverage"),labels=c("Conditional","Background refit"))]
g3 <- ggplot(compare,aes(size,coverage,shape=interval,colour=interval))+
    geom_hline(yintercept=.95,colour="#68756B",linetype=3,linewidth=.5)+
    geom_point(position=pd,size=3.4,stroke=1.0,fill="white")+
    scale_colour_manual(values=c("Conditional"="#315A78","Background refit"="#BA7241"))+
    scale_shape_manual(values=c(21,17))+scale_y_continuous(labels=scales::label_percent(accuracy=1),limits=c(0,1.02))+
    labs(x="Samples per group",y="Coverage for reinforcing changes")+base
mu <- exp(seq(log(.5),log(16),length.out=120))
curves <- rbindlist(list(data.table(mu=mu,source="Latent abundance",cv=.5),
    data.table(mu=mu,source="Observed, full depth",cv=.5+1/mu),
    data.table(mu=mu,source="Observed, quarter depth",cv=.5+4/mu)))
g4 <- ggplot(curves,aes(mu,cv,colour=source,linetype=source))+
    geom_line(linewidth=.85)+scale_x_log10(breaks=c(.5,1,2,4,8,16))+
    scale_y_log10(breaks=c(.5,1,2,4,8))+
    scale_colour_manual(values=c("Latent abundance"="#53685E","Observed, full depth"="#187F87","Observed, quarter depth"="#C56849"),labels=c("Latent","Full depth","Quarter depth"),guide=guide_legend(nrow=2,byrow=TRUE))+
    scale_linetype_manual(values=c("Latent abundance"="dotted","Observed, full depth"="solid","Observed, quarter depth"="longdash"),labels=c("Latent","Full depth","Quarter depth"),guide=guide_legend(nrow=2,byrow=TRUE))+
    labs(x="Latent mean abundance",y="CV² at fixed latent dispersion")+base+
    theme(legend.text=element_text(size=8))
fig <- (g1+g2+g3+g4)+plot_layout(ncol=2,guides="keep")+plot_annotation(tag_levels="a")
ggsave(file.path(dest,"Figure_S9.pdf"),fig,width=7.15,height=6.4,device=grDevices::cairo_pdf)
ggsave(file.path(dest,"Figure_S9.png"),fig,width=7.15,height=6.4,dpi=300,device=ragg::agg_png)
