construct_global_enrichment_map <- function(ora_res,
                                            which = c("GO", "KEGG"),
                                            top_n = 15,
                                            padj_max = 0.05,
                                            q_max = 1,
                                            sim = "jaccard",
                                            edge_min = 0.25,
                                            label_top = 30,
                                            layout = "fr",
                                            seed = 1,
                                            out_file = NULL) {
  which <- match.arg(which)
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("igraph must be installed to build enrichment maps.")
  }
  if (!requireNamespace("ggraph", quietly = TRUE)) {
    stop("ggraph must be installed to plot enrichment maps.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 must be installed to plot enrichment maps.")
  }
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("ggrepel must be installed to label enrichment maps.")
  }
  if (sim != "jaccard") {
    stop("Only 'jaccard' similarity is currently supported.")
  }

  split_gene_ids <- function(gene_string) {
    if (is.null(gene_string) || is.na(gene_string) || gene_string == "") {
      return(character(0))
    }
    parts <- unlist(strsplit(gene_string, "[/;]"))
    parts <- trimws(parts)
    parts <- parts[parts != ""]
    unique(parts)
  }

  `%||%` <- function(x, y) {
    if (is.null(x)) y else x
  }

  enrich_list <- if (which == "GO") ora_res$enrich_GO else ora_res$enrich_KEGG
  if (is.null(enrich_list) || length(enrich_list) == 0) {
    stop("No enrichment results found for ", which, ".")
  }
  if (!is.list(enrich_list)) {
    stop("Expected a list of enrichResult objects for ", which, ".")
  }

  term_rows <- list()
  for (i in seq_along(enrich_list)) {
    x <- enrich_list[[i]]
    if (is.null(x)) {
      message("Factor ", i, " is NULL; skipping.")
      next
    }
    df <- as.data.frame(x)
    if (nrow(df) == 0) {
      message("Factor ", i, " has 0 rows; skipping.")
      next
    }
    if (!"p.adjust" %in% names(df)) {
      message("Factor ", i, " missing p.adjust; skipping.")
      next
    }
    if (!"geneID" %in% names(df)) {
      message("Factor ", i, " missing geneID; skipping.")
      next
    }
    if (!"ID" %in% names(df)) {
      message("Factor ", i, " missing ID; skipping.")
      next
    }
    if (!"Description" %in% names(df)) {
      df$Description <- df$ID
    }

    df <- df[df$p.adjust <= padj_max, , drop = FALSE]
    if ("qvalue" %in% names(df)) {
      df <- df[df$qvalue <= q_max, , drop = FALSE]
    }
    if (nrow(df) == 0) {
      message("Factor ", i, " has no terms after filtering; skipping.")
      next
    }

    if ("Count" %in% names(df)) {
      ord <- order(df$p.adjust, -df$Count, na.last = TRUE)
    } else {
      ord <- order(df$p.adjust, na.last = TRUE)
    }
    df <- df[ord, , drop = FALSE]
    if (nrow(df) > top_n) {
      df <- df[seq_len(top_n), , drop = FALSE]
    }

    genes_split <- lapply(df$geneID, split_gene_ids)
    counts <- if ("Count" %in% names(df)) df$Count else rep(NA_integer_, nrow(df))
    idx <- length(term_rows) + 1
    term_rows[[idx]] <- data.frame(
      term_id = as.character(df$ID),
      description = as.character(df$Description),
      p.adjust = df$p.adjust,
      Count = counts,
      factor_index = i,
      stringsAsFactors = FALSE
    )
    term_rows[[idx]]$genes <- I(genes_split)
  }

  if (length(term_rows) == 0) {
    message("No terms selected after filtering for ", which, ".")
    empty_graph <- igraph::make_empty_graph()
    empty_plot <- ggplot2::ggplot() + ggplot2::theme_void()
    return(list(
      graph = empty_graph,
      nodes = data.frame(),
      edges = data.frame(),
      plot = empty_plot
    ))
  }

  term_df <- do.call(rbind, term_rows)
  term_split <- split(term_df, term_df$term_id)
  nodes_list <- lapply(term_split, function(df) {
    df <- df[order(df$p.adjust, df$factor_index), , drop = FALSE]
    genes_union <- unique(unlist(df$genes))
    factor_membership <- unique(df$factor_index)
    best_row <- df[1, , drop = FALSE]
    out <- data.frame(
      term_id = as.character(best_row$term_id),
      description = as.character(best_row$description),
      best_padj = best_row$p.adjust,
      best_factor = best_row$factor_index,
      Count = best_row$Count,
      n_factors = length(factor_membership),
      stringsAsFactors = FALSE
    )
    out$factor_membership <- I(list(paste0("F", factor_membership)))
    out$genes <- I(list(genes_union))
    out$gene_count <- length(genes_union)
    out
  })
  nodes <- do.call(rbind, nodes_list)
  nodes$shared <- nodes$n_factors > 1
  nodes$size_metric <- -log10(pmax(nodes$best_padj, .Machine$double.eps))
  nodes$size_metric <- pmin(nodes$size_metric, 10)
  factor_levels <- sort(unique(nodes$best_factor))
  nodes$best_factor <- factor(
    nodes$best_factor,
    levels = factor_levels,
    labels = paste0("F", factor_levels)
  )
  nodes$label <- NA_character_
  if (label_top > 0 && nrow(nodes) > 0) {
    top_idx <- order(nodes$best_padj, nodes$term_id)
    top_idx <- top_idx[seq_len(min(label_top, length(top_idx)))]
    nodes$label[top_idx] <- nodes$description[top_idx]
  }
  nodes$name <- nodes$term_id

  edges <- tryCatch({
    gene_sets <- nodes$genes
    gene_counts <- vapply(gene_sets, length, integer(1))
    eligible <- which(gene_counts >= 2)
    if (length(eligible) < 2) {
      return(data.frame(
        from = character(),
        to = character(),
        weight = numeric(),
        overlap = integer(),
        stringsAsFactors = FALSE
      ))
    }

    # Inverted index limits pairwise work compared to full O(n^2) scans.
    gene_map <- list()
    for (i in eligible) {
      genes <- unique(gene_sets[[i]])
      genes <- genes[!is.na(genes) & genes != ""]
      if (length(genes) == 0) next
      for (g in genes) {
        gene_map[[g]] <- c(gene_map[[g]], i)
      }
    }

    pair_counts <- new.env(hash = TRUE, parent = emptyenv())
    for (g in names(gene_map)) {
      terms <- unique(gene_map[[g]])
      if (length(terms) < 2) next
      combs <- utils::combn(terms, 2)
      for (k in seq_len(ncol(combs))) {
        i <- combs[1, k]
        j <- combs[2, k]
        key <- paste(i, j, sep = "|")
        pair_counts[[key]] <- (pair_counts[[key]] %||% 0L) + 1L
      }
    }

    keys <- ls(pair_counts)
    if (length(keys) == 0) {
      return(data.frame(
        from = character(),
        to = character(),
        weight = numeric(),
        overlap = integer(),
        stringsAsFactors = FALSE
      ))
    }

    edge_list <- vector("list", length(keys))
    edge_idx <- 0
    for (key in keys) {
      parts <- strsplit(key, "\\|")[[1]]
      i <- as.integer(parts[1])
      j <- as.integer(parts[2])
      overlap <- pair_counts[[key]]
      if (min(gene_counts[i], gene_counts[j]) < 2) next
      union <- gene_counts[i] + gene_counts[j] - overlap
      weight <- overlap / union
      if (weight >= edge_min) {
        edge_idx <- edge_idx + 1
        edge_list[[edge_idx]] <- data.frame(
          from = nodes$term_id[i],
          to = nodes$term_id[j],
          weight = weight,
          overlap = overlap,
          stringsAsFactors = FALSE
        )
      }
    }
    if (edge_idx == 0) {
      data.frame(
        from = character(),
        to = character(),
        weight = numeric(),
        overlap = integer(),
        stringsAsFactors = FALSE
      )
    } else {
      do.call(rbind, edge_list[seq_len(edge_idx)])
    }
  }, error = function(e) {
    message("Similarity computation failed: ", e$message)
    data.frame(
      from = character(),
      to = character(),
      weight = numeric(),
      overlap = integer(),
      stringsAsFactors = FALSE
    )
  })

  if (nrow(edges) == 0) {
    message("No edges met the similarity threshold; plotting nodes only.")
  }

  graph <- igraph::graph_from_data_frame(
    d = edges,
    directed = FALSE,
    vertices = nodes
  )

  plot <- tryCatch({
    set.seed(seed)
    p <- ggraph::ggraph(graph, layout = layout)
    if (nrow(edges) > 0) {
      p <- p + ggraph::geom_edge_link(
        ggplot2::aes(width = weight, alpha = weight),
        color = "grey50",
        show.legend = TRUE
      )
    }
    p <- p +
      ggraph::geom_node_point(
        ggplot2::aes(size = size_metric, fill = best_factor, color = shared),
        shape = 21,
        stroke = 0.8
      ) +
      ggrepel::geom_text_repel(
        ggplot2::aes(label = label),
        size = 3,
        max.overlaps = Inf
      ) +
      ggplot2::scale_fill_brewer(palette = "Set2") +
      ggplot2::scale_color_manual(
        values = c(`FALSE` = "grey60", `TRUE` = "black")
      ) +
      ggraph::scale_edge_width(range = c(0.2, 1.5)) +
      ggraph::scale_edge_alpha(range = c(0.2, 0.8)) +
      ggplot2::labs(
        fill = "Best factor",
        color = "Shared term",
        size = "-log10(adj p)",
        edge_width = "Jaccard",
        edge_alpha = "Jaccard",
        caption = "Shared terms have black outlines."
      ) +
      ggplot2::guides(
        color = ggplot2::guide_legend(override.aes = list(fill = "white"))
      ) +
      ggplot2::theme_void()
    p
  }, error = function(e) {
    message("Plotting failed: ", e$message)
    ggplot2::ggplot() + ggplot2::theme_void()
  })

  if (!is.null(out_file)) {
    padj_tag <- gsub("\\.", "p", format(padj_max, trim = TRUE))
    edge_tag <- gsub("\\.", "p", format(edge_min, trim = TRUE))
    suffix <- paste0(which, "_top", top_n, "_padj", padj_tag, "_edge", edge_tag)
    base <- tools::file_path_sans_ext(out_file)
    ggplot2::ggsave(
      filename = paste0(base, "_", suffix, ".pdf"),
      plot = plot,
      width = 10,
      height = 8
    )
    ggplot2::ggsave(
      filename = paste0(base, "_", suffix, ".png"),
      plot = plot,
      width = 10,
      height = 8,
      dpi = 300
    )
  }

  list(
    graph = graph,
    nodes = nodes,
    edges = edges,
    plot = plot
  )
}
