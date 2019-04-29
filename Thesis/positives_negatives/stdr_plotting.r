library(ggplot2)

for( i in 1:length(colnames(chi)) ){
  gene <- colnames(chi[i])
  y_axis_name <- paste0(gene," mutations chi-like frequencies")

  tmp_fig <- ggplot(data=chi, aes(x=rownames(chi), y=log2(chi[[gene]]), color=rownames(chi)))+
    geom_bar(stat="identity", fill="white")+
    geom_text(aes(label=signif(log2(chi[[gene]]), digits=5)), vjust=1.6, color="black", size=3.5)+
    xlab("Cancer types")+
    ylab(y_axis_name)+
    theme(legend.position="none")

  tmp_path <- paste0("C:/Users/INTEL/Documents/Atom/MC_Thesis/positives_negatives/figures/",gene,"_chi_dist.png")
  print(tmp_fig)
  invisible(dev.print(png, tmp_path, width=1366, height=650))
  invisible(dev.off())
}
