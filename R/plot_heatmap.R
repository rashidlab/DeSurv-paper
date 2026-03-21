plot_heatmap <- function(
    x, k,k_row,
    results = NULL, data = NULL,   # only needed if you want to rebuild Xtemp
    Xtemp = NULL,
    palette = c("darkblue","blue","white","red","darkred"),
    scale_rows = TRUE
) {
  if (missing(k) || length(k) != 1L || !is.numeric(k) || k < 2) stop("`k` must be a single integer >= 2.")
  if (length(x$clusCol) < k || is.null(x$clusCol[[k]])) {
    avail <- paste(seq_along(x$clusCol), collapse = ", ")
    stop(sprintf("clusCol[[%d]] not available. Available K: %s", k, avail))
  }
  
  if (is.null(Xtemp)) {
    if (is.null(results) || is.null(data)) stop("Provide `Xtemp`, or `results`+`data` to rebuild it.")
    rebuilt <- build_expression_matrix(results, data, facs = x$facs, weight = FALSE)
    Xtemp <- rebuilt$Xtemp
  }
  if (!is.matrix(Xtemp)) Xtemp <- as.matrix(Xtemp)
  
  if (isTRUE(scale_rows)) {
    cn <- colnames(Xtemp); Xtemp <- t(apply(Xtemp, 1, scale)); colnames(Xtemp) <- cn
  }
  
  dend_col <- as.dendrogram(x$clusCol[[k]]$consensusTree)
  samp_order <- colnames(Xtemp)[labels(dend_col)]
  
  clus_order = x$clusCol[[k]]$consensusClass[samp_order]
  clus_order = c(clus_order[clus_order==1],
                 clus_order[clus_order==2],
                 clus_order[clus_order==3],
                 clus_order[clus_order==4],
                 clus_order[clus_order==5])
  samp_order = names(clus_order)
  
  if (!length(samp_order)) stop("No overlap between dendrogram labels and Xtemp columns.")
  Xtemp <- Xtemp[, samp_order, drop = FALSE]
  
  gene_anno <- NULL
  if (!is.null(x$clusRow)) {
    if (!is.null(x$clusRow[[k_row]])) {
      dend_row <- as.dendrogram(x$clusRow[[k_row]]$consensusTree)
      gene_order <- names(x$clusRow[[k_row]]$consensusClass)[labels(dend_row)]
      
      clus_order = x$clusRow[[k_row]]$consensusClass[gene_order]
      clus_order = c(clus_order[clus_order==1],
                     clus_order[clus_order==2],
                     clus_order[clus_order==3],
                     clus_order[clus_order==4],
                     clus_order[clus_order==5])
      gene_order = names(clus_order)
      genes = intersect(gene_order,rownames(Xtemp))
      Xtemp <- Xtemp[genes, , drop = FALSE]
      gcl <- x$clusRow[[k_row]]$consensusClass
      gene_anno <- data.frame(Gene = names(gcl), GeneCluster = factor(as.integer(gcl)))
      rownames(gene_anno) = gene_anno$Gene
      gene_anno = gene_anno[genes,]
    }
  }
  
  scl <- x$clusCol[[k]]$consensusClass
  sample_anno <- data.frame(Sample = names(scl), Cluster = factor(as.integer(scl)))
  rownames(sample_anno) = sample_anno$Sample
  sample_anno = sample_anno[samp_order,]
  
  df <- tibble::as_tibble(Xtemp, rownames = "Gene") |>
    tidyr::pivot_longer(-Gene, names_to = "Sample", values_to = "Expression")
  
  df$Sample <- factor(df$Sample, levels = colnames(Xtemp))
  df$Gene   <- factor(df$Gene,   levels = rownames(Xtemp))
  sample_anno$Sample <- factor(sample_anno$Sample, levels = colnames(Xtemp))
  if (!is.null(gene_anno)) gene_anno$Gene <- factor(gene_anno$Gene, levels = rownames(Xtemp))
  
  rng <- max(abs(df$Expression), na.rm = TRUE)
  min_val <- -rng; max_val <- rng
  
  p <- ggplot2::ggplot(df, ggplot2::aes(Sample, Gene, fill = Expression)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(
      colors = palette,
      values = scales::rescale(c(min_val, -3, 0, 3, max_val)),
      limits = c(min_val, max_val)
    ) +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_tile(
      data = sample_anno,
      ggplot2::aes(x = Sample, y = 0, fill = Cluster),
      inherit.aes = FALSE, height = 1
    ) +
    ggnewscale::new_scale_fill()
  
  if (!is.null(gene_anno) && nrow(gene_anno) > 0L) {
    p <- p + ggplot2::geom_tile(
      data = gene_anno,
      ggplot2::aes(x = 0, y = Gene, fill = GeneCluster),
      inherit.aes = FALSE, width = 1
    )
  }
  
  p +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(size = 6)) +
    ggplot2::labs(title = x$title, x = NULL, y = NULL)
}
