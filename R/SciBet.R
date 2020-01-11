.onLoad <- function(libname, pkgname) {
  op <- options()
  op.scibet <- list(
    scibet.n_threads = 0
  )
  toset <- !(names(op.scibet) %in% names(op))
  if(any(toset)) options(op.scibet[toset])
  invisible()
}

#' Select informative genes of the training set in a supervised manner.
#' @name SelectGene
#' @usage SelectGene(expr, k=1000, r = FALSE)
#' @param expr The expression dataframe, with rows being cells, and columns being genes. The last column should be "label".
#' @param k Number of genes to select according to the total entropy differneces among cell types.
#' @param r Whether to reture a datafram of the gene list. The default is FALSE.
#' @return A list of genes informative genes.
#' @export
SelectGene <- function(expr, k = 1000, r = FALSE) {
  as_df <- r
  n_threads <- getOption("scibet.n_threads")
  labels <- factor(expr$label)
  labels_in <- as.integer(labels) - 1
  out <- DsGene(expr, labels_in, as_df, k, 0, n_threads)
  if (r) {
    out <- as.data.frame(out)
    colnames(out) <- c(levels(labels), "Total")
    rownames(out) <- colnames(expr)[-ncol(expr)]
    out <- out %>% tibble::rownames_to_column(var = "gene") %>% as.tibble() %>%
      dplyr::arrange(desc(Total))
    type_ds_sub <- out %>% dplyr::select(-gene, -Total)
    type_label <- colnames(type_ds_sub)
    get_label <- function(.x) {
      type_label[which(type_ds_sub[.x, ] == max(type_ds_sub[.x,
                                                            ]))[1]]
    }
    gene_label <- tibble(gene_num = 1:nrow(out), gene = out$gene) %>%
      dplyr::mutate(label = purrr::map_chr(gene_num, get_label))

    gene_label <- gene_label %>%
      dplyr::select(-gene_num) %>%
      tidyr::nest(-label) %>%
      dplyr::mutate(
        data = purrr::map(
          .x = data,
          .f = function(.x) {
            .x %>% dplyr::mutate(flag = 1:nrow(.))}))

    gene_label <- Reduce(rbind, gene_label$data)
    out <- out %>% dplyr::inner_join(gene_label, by = "gene") %>%
      dplyr::arrange(flag, desc(Total))
    return(out$gene[1:k])
  }
  return(colnames(expr)[out + 1])
}

#' Train SciBet model and generate a "Bet" function for cell type prediction.
#' @name Learn
#' @usage Learn(expr, geneset, k=1000, a=5)
#' @param expr The reference expression dataframeThe expression dataframe, with rows being cells, and columns being genes. The last column should be "label".
#' @param geneset A user-defined set of genes for prediction.
#' @param k See [SelectGene()] for details.
#' @param a For each cell type, select a genes by perform E-test in the one-vs-rest manner. If users set a>0, the param k will be neglected.
#' @return A function that takes two inputs: `expr` and `result`. See [Bet()] for details.
#' @export
Learn <- function(expr, geneset=NULL, k=1000, a=0){
  n_threads <- getOption("scibet.n_threads")
  labels <- factor(expr$label)
  labels_in <- as.integer(labels) - 1
  if(is.null(geneset)){
    geneset <- DsGene(expr, labels_in, FALSE, k, a, n_threads)
  }
  prob <- GenProb(expr, labels_in, geneset, n_threads)
  genes <- colnames(expr)[geneset + 1]
  labell <- levels(labels)
  rm(expr, labels, geneset, labels_in)
  function(expr, result="list"){
    have_genes <- which(genes %in% colnames(expr))
    expra <- log1p(as.matrix(expr[, genes[have_genes]])) / log(2)
    switch(result,
      list = labell[Gambler(expra, prob[have_genes, ], FALSE, n_threads) + 1],
      table = {
        out <- Gambler(expra, prob[have_genes, ], TRUE, n_threads)
        colnames(out) <- labell
        return(out)
      }
    )
  }
}

#' Classify cells of a given query dataset using a reference dataset.
#' @description SciBet main function. Train SciBet with the reference dataset to assign cell types for the query dataset.
#' @name SciBet
#' @usage SciBet(train, test, k=1000, result=c("list", "table"))
#' @param train The reference dataset, with rows being cells, and columns being genes. The last column should be "label".
#' @param test The query dataset. Rows should be cells and columns should be genes.
#' @param k Number of genes to select according to the total entropy differneces among cell types.
#' @examples SciBet(train.matr, query.matr)
#' @export
SciBet <- function(train, test, k=1000, result="list"){
  Learn(train, NULL, k, 0)(test, result)
}

#' Compute expression entropy.
#' @name Entropy
#' @usage Entropy(expr, window=120, low = 2000)
#' @param expr The expression dataframe. Rows should be cells and columns should be genes.
#' @param window The window size for expression value discretization.
#' @param low The lower limit for normalizing expression entropy
#' @return A dataframe..
#' @export
Entropy <- function(expr, window=120, low = 2000){
  n_threads <- getOption("scibet.n_threads")
  expr <- as.data.frame(expr)

  ent_res <- tibble(
    gene = colnames(expr),
    mean.expr = colMeans(expr)
  ) %>%
    dplyr::filter(mean.expr < 6000)

  expr <- expr[,ent_res$gene]
  out <- GenEntr(expr, window, n_threads)

  ent_res %>%
    dplyr::mutate(entropy = out) %>%
    dplyr::mutate(fit = 0.18*log(0.03*mean.expr + 1)) -> ent_res

  ent_res %>%
    dplyr::filter(mean.expr > low) %>%
    dplyr::mutate(k = entropy/fit) %>%    #linear normalization of expression entropy
    dplyr::pull(k) %>%
    quantile(0.75) %>%
    as.numeric() -> k

  ent_res <- ent_res %>% dplyr::mutate(norm_ent = entropy/k)

  return(ent_res)
}

#' S-E curve.
#' @name SEM
#' @usage SEM(ent_res)
#' @param ent_res The computing result of expression entropy implemented with Entropy function.
#' @return A figure.
#' @export
SEM <- function(ent_res){
  ent_res %>%
    ggplot(aes(mean.expr, norm_ent)) +
    geom_point(colour = "#1E90FF") +
    geom_line(aes(mean.expr, fit), colour = 'black', lwd = 0.9) +
    theme_classic() +
    theme(
      legend.position = 'none',
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 15),
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 15),
      axis.text.y = element_text(color="black"),
      axis.text.x = element_text(color="black")
    ) +
    labs(
      x = "mean expression (E)",
      y = "expression entropy (S)"
    ) +
    ylim(0,1) -> p

  return(p)
}


#' Point plot for entropy reduction.
#' @name DsPlot
#' @usage DsPlot(ent_res, cutoff = 0, gene_num = 10)
#' @param ent_res The computing result of expression entropy implemented with Entropy function.
#' @param cutoff The cutoff of entropy reduction.
#' @param gene_num Top genes with maximal entropy reduction.
#' @return A figure.
#' @export
DsPlot <- function(ent_res, cutoff = 0, gene_num = 10){
  ent_res %>%
    dplyr::mutate(ds = fit - norm_ent) %>%
    dplyr::arrange(desc(ds)) %>%
    dplyr::filter(ds > cutoff) %>%
    dplyr::mutate(sig = as.character(1:nrow(.))) %>%
    dplyr::mutate(sig = ifelse(sig %in% as.character(1:gene_num), sig, as.character(0))) %>%
    dplyr::mutate(Gene = ifelse(sig %in% as.character(1:gene_num), gene, NA)) %>%
    dplyr::arrange(gene) %>%
    ggplot(aes(1:nrow(.), ds)) +
    geom_point(aes(colour = sig)) +
    scale_colour_manual(values = c("grey50", rep('red', gene_num))) +
    geom_text(aes(label = Gene), vjust = 0.5, hjust = -0.1) +
    theme_classic() +
    theme(
      legend.position = 'none',
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 15),
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 15),
      axis.text.y = element_text(color="black"),
      axis.text.x = element_text(color="black")
    ) +
    labs(
      x = "Genes",
      y = "S-reduction"
    )  -> p

  return(p)
}

#' Expression patterns of informative genes across cell types.
#' @name Marker_heatmap
#' @usage Marker_heatmap(expr, gene)
#' @param expr The expression dataframe. Rows should be cells, columns should be genes and last column should be "label".
#' @param gene A vector of informative genes.
#' @return A figure.
#' @export
Marker_heatmap <- function(expr, gene){
  expr <- expr[,c(gene,'label')]
  type_expr <- expr %>%
    tidyr::nest(-label) %>%
    dplyr::rename(expr = data) %>%
    dplyr::mutate(colmeans = purrr::map(
      .x = expr,
      .f = function(.x){colMeans(.x)}))

  type_expr$colmeans %>%
    as.data.frame() %>%
    tibble::remove_rownames() %>%
    t() %>%
    as.data.frame() %>%
    tibble::remove_rownames() -> type_mean_expr

  rownames(type_mean_expr) <- type_expr$label
  colnames(type_mean_expr) <- colnames(expr)[-ncol(expr)]

  sub_expr <- type_mean_expr
  sub_expr <- sub_expr %>%
    as.tibble() %>%
    dplyr::mutate_all(funs((. - mean(.))/sd(.))) %>%
    t()
  colnames(sub_expr) <- type_expr$label
  get_label <- function(num){
    v <- sub_expr[num,]
    colnames(sub_expr)[which(v == max(v))]
  }
  sub_expr <- sub_expr %>%
    tibble::as.tibble() %>%
    dplyr::mutate(group = purrr::map_chr(1:length(gene), get_label))
  sub_expr <- as.data.frame(sub_expr)
  rownames(sub_expr) <- gene
  sub_expr <- sub_expr %>%
    dplyr::mutate(gene = gene) %>%
    tidyr::gather(key = 'cell_type', value = 'zscore', -group, -gene) %>%
    dplyr::arrange(group, desc(zscore))
  sub_expr %>%
    ggplot(aes(factor(gene, levels = unique(sub_expr$gene)),
               factor(cell_type, levels = sort(unique(sub_expr$cell_type), decreasing = T)))) +
    geom_point(aes(size = zscore, colour = zscore)) +
    theme(
      strip.text.x = element_blank(),
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 13),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 13),
      axis.text.y = element_text(color="black"),
      axis.text.x = element_text(color="black", angle = -90, hjust = 0),
      panel.background = element_rect(colour = "black", fill = "white"),
      panel.grid = element_line(colour = "grey", linetype = "dashed"),
      panel.grid.major = element_line(
        colour = "grey",
        linetype = "dashed",
        size = 0.2
      )
    ) +
    facet_grid(. ~ group, scales = "free", space = "free") +
    scale_colour_distiller(palette = "RdYlBu") +
    labs(
      x = '',
      y = ''
    ) -> p

  return(p)
}

#' Heatmap of classification result.
#' @name Confusion_heatmap
#' @usage Confusion_heatmap(ori, prd)
#' @param ori A vector of the original labels for each cell in the test set.
#' @param prd A vector of the predicted labels for each cell in the test set..
#' @return A heatmap for the confusion matrix of the classification result.
#' @export
Confusion_heatmap <- function(ori, prd){
  tibble(
    ori = ori,
    prd = prd
  ) %>%
    dplyr::count(ori, prd) %>%
    tidyr::spread(key = prd, value = n) -> cross.validation.filt

  cross.validation.filt[is.na(cross.validation.filt)] = 0
  cross.validation.filt[,-1] <- round(cross.validation.filt[,-1]/rowSums(cross.validation.filt[,-1]),2)
  cross.validation.filt <- cross.validation.filt %>%
    tidyr::gather(key = 'prd', value = 'Prob', -ori)

  cross.validation.filt %>%
    ggplot(aes(ori,prd,fill = Prob)) +
    geom_tile() +
    theme(axis.title = element_text(size = 0)) +
    theme(axis.text = element_text(size = 10)) +
    theme(legend.title = element_text(size = 0)) +
    theme(legend.text = element_text(size = 10)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank()) +
    theme(axis.text.y = element_text(color="black"),
          axis.text.x = element_text(color="black",angle = 45, hjust = 1)) +
    scale_fill_viridis() -> p

  return(p)
}

#' Calculating the confidence score C for false positive control.
#' @name conf_score
#' @usage conf_score(ref, query, null_expr, gene_num = 500)
#' @param ref The reference dataset. Rows should be cells, columns should be genes and last column should be "label".
#' @param query The query dataset. Rows should be cells and columns should be genes.
#' @param gene_num The number of common markers of reference set used for false positive control.
#' @return A vector of confidence scores.
#' @export
conf_score <- function(ref, query, null_expr, gene_num)
{
  n_threads <- getOption("scibet.n_threads")
  labels <- factor(ref$label)
  labels_in <- as.integer(labels) - 1
  genes <- Reduce(intersect, list(colnames(ref), colnames(query), names(null_expr)))
  ref <- log1p(as.matrix(ref[, genes])) / log(2)
  query <- log1p(as.matrix(query[, genes])) / log(2)
  a <- NullTest(ref, query, null_expr[genes], labels_in, gene_num, n_threads)
  b <- NullTest(ref, ref, null_expr[genes], labels_in, gene_num, n_threads)
  prob <- a/max(b)
  prob[prob > 1] <- 1
  return(prob)
}

#' Heatmap for the confusion matrix of the classification with the false postive control.
#' @name Confusion_heatmap_negctrl
#' @usage Confusion_heatmap_negctrl(res, cutoff = 0.4)
#' @param res Classification result.
#' @param cutoff The cutoff of confifence score C.
#' @return A heatmap for the confusion matrix of the classification result with the false postive control.
#' @export
Confusion_heatmap_negctrl <- function(res, cutoff = 0.4){
  res %>%
    dplyr::mutate(prd = ifelse(c_score < cutoff, 'unassigned', prd)) %>%
    dplyr::count(ori, prd) %>%
    tidyr::spread(key = prd, value = n) -> cla.res

  cla.res[is.na(cla.res)] = 0
  cla.res[,-1] <- round(cla.res[,-1]/rowSums(cla.res[,-1]),2)
  cla.res <- cla.res %>% tidyr::gather(key = 'prd', value = 'Prob', -ori)
  label <- cla.res$ori %>% unique()
  cla.res %>%
    ggplot(aes(prd, factor(ori, levels = c(label[-3],'Neg.cell')), fill = Prob)) +
    geom_tile(colour = 'white', lwd = 0.5) +
    theme(axis.title = element_text(size = 12)) +
    theme(axis.text = element_text(size = 12)) +
    theme(legend.title = element_text(size = 0)) +
    theme(legend.text = element_text(size = 12)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank()) +
    theme(axis.text.y = element_text(color="black"),
          axis.text.x = element_text(color="black", angle = 50, hjust = 1)) +
    scale_fill_material('blue')
}

#' Make predictions with the trained SciBet model.
#' @name Bet
#' @usage Bet(expr, result=c("list", "table"))
#' @param expr The expression matrix or dataframe for prediction. Rows should be cells and columns should be genes.
#' @param result Return a "list" of predicted labels, or a "table" of probabilities of each tested cell belonging to each label.
#' @return A vector or a dataframe for the classification result.
#' @export
Bet <- function(expr, result="list"){
}

#' Export model as matrix.
#' @details If you don't plan to export the models for usage on other platforms, simply saving the Bet function in a RData file would also work in R.
#' @name ExportModel
#' @usage ExportModel(Bet)
#' @param Bet A prediction function generated by SciBet.
#' @return A matrix.
#' @export
ExportModel <- function(Bet){
  prob <- environment(Bet)$prob
  rownames(prob) <- environment(Bet)$genes
  colnames(prob) <- environment(Bet)$labell
  return(prob)
}

#' Generate Bet function from a model matrix.
#' @name LoadModel
#' @usage LoadModel(x, genes, labels)
#' @param x A SciBet model in the format of a matrix.
#' @param genes (Optional).
#' @param labels (Optional).
#' @return A Bet function.
#' @export
LoadModel <- function(x, genes=NULL, labels=NULL){
  n_threads <- getOption("scibet.n_threads")
  prob <- x
  if (is.null(genes))
    genes <- rownames(x)
  if (is.null(labels))
    labels <- colnames(x)
  function(expr, result="list"){
    have_genes <- which(genes %in% colnames(expr))
    expra <- log1p(as.matrix(expr[, genes[have_genes]])) / log(2)
    switch(result,
           list = labels[Gambler(expra, prob[have_genes, ], FALSE, n_threads) + 1],
           table = {
             out <- Gambler(expra, prob[have_genes, ], TRUE, n_threads)
             colnames(out) <- labels
             return(out)
           }
    )
  }
}

#' Process scibet.core
#' @name pro.core
#' @usage pro.core(scibet.core)
#' @param scibet.core A SciBet core
#' @return A processed SciBet core
#' @export
pro.core <- function(scibet.core){
  cell.type <- unname(unlist(scibet.core[,1]))
  scibet.core <- as.data.frame(t(scibet.core[,-1]))
  colnames(scibet.core) <- cell.type
  return(as.matrix(scibet.core))
}
