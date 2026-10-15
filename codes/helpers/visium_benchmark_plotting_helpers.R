# Runtime/count layouts extracted from the 12 October 2026 Visium report.
# These functions draw only; their data are supplied by the plotting entry point.
runtime_plot <- function() {
  shown <- performance
  memory_max <- max(performance$peak_rss_gib)
  runtime_max <- max(performance$elapsed_minutes,1)
  memory_breaks <- pretty(c(0,memory_max*1.25),n=4)
  runtime_breaks <- pretty(c(0,runtime_max*1.5),n=4)
  p <- ggplot(shown,aes(peak_rss_gib,elapsed_minutes,colour=condition,
    shape=condition,fill=condition))+
    geom_point(size=3.8,stroke=0.8)+
    geom_text(aes(label=sprintf("%.2f min\n%.2f GiB",elapsed_minutes,peak_rss_gib),
      x=peak_rss_gib+ifelse(condition=="spatialdm_permutation",-1,1)*max(memory_breaks)*0.055,
      y=pmax(elapsed_minutes,max(runtime_breaks)*0.04),
      hjust=ifelse(condition=="spatialdm_permutation",1,0)),
      size=2.6,lineheight=0.95,vjust=0.5,show.legend=FALSE)+
    scale_colour_manual(values=colours,breaks=methods,labels=labels,drop=FALSE)+
    scale_shape_manual(values=shapes,breaks=methods,labels=labels,drop=FALSE)+
    scale_fill_manual(values=colours,breaks=methods,labels=labels,drop=FALSE)+
    scale_x_continuous(limits=range(memory_breaks),
      breaks=memory_breaks,
      labels=scales::label_comma(),expand=expansion(mult=c(0.035,0.02)))+
    scale_y_continuous(limits=range(runtime_breaks),
      breaks=runtime_breaks,
      labels=scales::label_comma(),expand=expansion(mult=c(0.035,0.02)))+
    labs(x="Peak RAM (GiB)",y="Runtime (min)",colour=NULL,shape=NULL,fill=NULL)+
    guides(colour=guide_legend(nrow=2,byrow=TRUE),
      shape=guide_legend(nrow=2,byrow=TRUE),fill=guide_legend(nrow=2,byrow=TRUE))+
    theme_bw(base_size=11,base_family="sans")+
    theme(aspect.ratio=1,panel.border=element_rect(colour="black",fill=NA,linewidth=0.7),
      panel.grid=element_blank(),axis.ticks=element_line(linewidth=0.5),
      axis.text=element_text(face="bold",colour="black",size=11),
      axis.title=element_text(face="bold",size=12),plot.title=element_text(face="bold",size=12),
      legend.position="bottom",legend.text=element_text(size=10),
      legend.key.size=grid::unit(0.18,"in"),legend.spacing.x=grid::unit(0.05,"in"),
      legend.margin=margin(0,0,0,0),legend.box.spacing=grid::unit(0.04,"in"),
      plot.margin=margin(4,4,4,4))
  p
}
linked_plot <- function(ids,title,show_names=TRUE) {
  shown <- counts[counts$lr_id %in% ids,]
  shown$method <- factor(shown$condition,levels=methods)
  shown$x <- match(shown$lr_id,ids)
  dot_max <- max(c(1,counts$n_sig[counts$lr_id %in% ids]),na.rm=TRUE)
  available <- shown[is.finite(shown$n_sig),]
  missing <- shown[!is.finite(shown$n_sig),]
  has_missing <- nrow(missing)>0L
  segments <- do.call(rbind,lapply(methods,function(method) {
    x <- shown[shown$condition==method,]
    x <- x[order(x$x),]
    x$x_next <- c(x$x[-1],NA_real_)
    x$y_next <- c(x$n_sig[-1],NA_real_)
    x[is.finite(x$n_sig)&is.finite(x$y_next),]
  }))
  dot_breaks <- pretty(c(0,dot_max),n=5)
  dot_breaks <- dot_breaks[dot_breaks>=0 & dot_breaks<=dot_max*1.07]
  offsets <- setNames(c(-0.18,0,0.18),methods)
  p <- ggplot(available,aes(x,n_sig,colour=method))+
    (if(has_missing) annotate("rect",xmin=0.45,xmax=length(ids)+0.55,ymin=-0.17*dot_max,
      ymax=-0.035*dot_max,fill="#F2F2F2",colour=NA))+
    geom_hline(yintercept=0,colour="#777777",linewidth=0.35)+
    geom_segment(data=segments,aes(xend=x_next,yend=y_next,linetype=method),linewidth=0.65)+
    geom_point(aes(shape=method,fill=method),size=if(show_names)2.6 else 0.65,stroke=if(show_names)0.75 else 0.25)+
    geom_point(data=missing,aes(x=x+offsets[as.character(method)],y=-0.11*dot_max,
      colour=method,shape=method,fill=method),inherit.aes=FALSE,
      size=if(show_names)2.1 else 0.65,stroke=0.65,show.legend=FALSE)+
    scale_colour_manual(values=colours,breaks=methods,labels=labels,drop=FALSE)+
    scale_shape_manual(values=shapes,breaks=methods,labels=labels,drop=FALSE)+
    scale_fill_manual(values=colours,breaks=methods,labels=labels,drop=FALSE)+
    scale_linetype_manual(values=setNames(rep("solid",3),methods),breaks=methods,labels=labels,drop=FALSE)+
    scale_y_continuous(breaks=c(if(has_missing) -0.11*dot_max,dot_breaks),
      labels=c(if(has_missing) "NA",scales::label_comma()(dot_breaks)),
      limits=c(if(has_missing) -0.18*dot_max else -0.035*dot_max,1.07*dot_max),
      expand=expansion(mult=0))+
    labs(title=title,x=if(show_names)NULL else "LR rank (10k blisa hotspot count)",
      y="Significant spots",colour=NULL,shape=NULL,fill=NULL,linetype=NULL)+
    guides(colour=guide_legend(nrow=1),shape=guide_legend(nrow=1),
      fill=guide_legend(nrow=1),linetype=guide_legend(nrow=1))+
    theme_bw(base_size=11)+
    theme(text=element_text(face="bold",colour="#222222"),
      panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),
      panel.grid.major.y=element_line(colour="#E5E5E5",linewidth=0.3),
      panel.border=element_rect(colour="#555555",linewidth=0.5),
      axis.text.x=element_text(angle=if(show_names)55 else 0,hjust=if(show_names)1 else 0.5,
        vjust=1,size=11,face="bold",colour="#222222"),
      axis.text.y=element_text(size=11,face="bold",colour="#222222"),
      axis.title.y=element_text(size=12,face="bold",margin=margin(r=5)),
      legend.position="bottom",legend.justification="left",legend.text=element_text(size=11,face="bold"),
      legend.key.width=grid::unit(0.8,"cm"),legend.key.height=grid::unit(0.4,"cm"),
      legend.margin=margin(0,0,0,0),legend.box.spacing=grid::unit(1,"mm"),
      plot.title=element_text(face="bold",size=13,margin=margin(b=6)),
      plot.margin=margin(6,8,4,5))
  if(show_names) p <- p+scale_x_continuous(breaks=seq_along(ids),
    labels=gsub("__"," - ",ids,fixed=TRUE),limits=c(0.45,length(ids)+0.55),expand=expansion(mult=0))
  else p <- p+scale_x_continuous(limits=c(0.45,length(ids)+0.55),expand=expansion(mult=0))
  p
}
