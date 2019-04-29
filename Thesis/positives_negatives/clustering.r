for( i in 1:length(colnames(chi)) ){
  values <- chi[[i]]
  model <- kmeans(values, 2)
  df <- data.frame(names=rownames(chi), values=values, clusters=model$cluster, stringsAsFactors=FALSE)

  gene <- colnames(chi[i])
  y_axis_name <- paste0(gene," chi values")

  tmp_fig <- ggplot(data=df, aes(x=names, y=log2(df[[2]]), color=clusters))+
      geom_point()+
      xlab("Cancer types")+
      ylab(y_axis_name)+
      labs(colour="Clusters")

  tmp_path <- paste0("C:/Users/INTEL/Documents/Atom/MC_Thesis/positives_negatives/figures/",gene,"_clustering.png")
  print(tmp_fig)
  invisible(dev.print(png, tmp_path, width=1366, height=650))
  invisible(dev.off())
}
