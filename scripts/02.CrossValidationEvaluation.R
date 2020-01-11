library(SingleCellExperiment)
library(reticulate)
library(tidyverse)
library(Seurat)
library(scmap)
library(fmsb)
library(scibet)

file.name <- "GSE75140.rds.gz"
file_path <- "/home/pauling/projects/01_classifier/01_data/03_FeatureSelectEval/new/matrix"
out.path <- "/home/pauling/projects/01_classifier/01_data/15_revise/01.crossvalidation"

expr <- readr::read_rds(file.path(file_path, file.name))
expr[,-ncol(expr)] <- 1000000*expr[,-ncol(expr)]/rowSums(expr[,-ncol(expr)])
expr <- as.data.frame(expr)


scmap_ck <- function(expr_train, expr_test, gene_num = 500){

  train_label <- expr_train$label
  test_label <- expr_test$label

  expr_train <- expr_train %>% dplyr::select(-label) %>% t()
  expr_test <- expr_test %>% dplyr::select(-label) %>% t()

  ann <- data.frame(cell_type1 = train_label)
  sce <- SingleCellExperiment(assays = list(normcounts = as.matrix(expr_train)), colData = ann)
  logcounts(sce) <- log2(normcounts(sce) + 1)
  rowData(sce)$feature_symbol <- rownames(sce)
  sce <- sce[!duplicated(rownames(sce)), ]

  tx_sce <- SingleCellExperiment(assays = list(normcounts = as.matrix(expr_test)))
  logcounts(tx_sce) <- log2(normcounts(tx_sce) + 1)
  rowData(tx_sce)$feature_symbol <- rownames(tx_sce)

  sce <- selectFeatures(sce, n_features = gene_num, suppress_plot = T)
  sce <- indexCluster(sce)
  scmapCluster_results <- scmapCluster(
    projection = tx_sce,
    threshold = 0,
    index_list = list(
      yan = metadata(sce)$scmap_cluster_index
    )
  )

  tibble(
    ori = as.character(test_label),
    prd = unlist(scmapCluster_results$combined_labs),
    prob = scmapCluster_results$scmap_cluster_siml[,1],
    method = 'scmap')

}
SciBet_ck <- function(expr_train, expr_test, gene_num = 500){

  prd <- SciBet(expr_train, expr_test[,-ncol(expr_test)], k = gene_num, a = 0)

  tibble(
    ori = as.character(expr_test$label),
    prd = prd,
    method = 'SciBet')
}
Seurat3 <- function(expr_train, expr_test, gene_num = 500){

  data.frame(
    celltype = expr_train$label,
    tech = 'xx'
  ) -> metadata

  data.frame(
    celltype = expr_test$label,
    tech = 'yy'
  ) -> metadata1

  ori <- expr_test$label

  X_train <- as.matrix(t(expr_train[,-ncol(expr_train)]))
  X_test <- as.matrix(t(expr_test[,-ncol(expr_test)]))

  X_train <- X_train

  matr <- cbind(X_train, X_test)
  metadata <- rbind(metadata, metadata1)
  colnames(matr) <- as.character(1:ncol(matr))
  rownames(metadata) <- as.character(1:nrow(metadata))

  ttest <- CreateSeuratObject(counts = matr, meta.data = metadata)
  ttest.list <- SplitObject(object = ttest, split.by = "tech")

  for (i in 1:length(x = ttest.list)) {
    ttest.list[[i]] <- NormalizeData(object = ttest.list[[i]], verbose = FALSE)
    ttest.list[[i]] <- FindVariableFeatures(object = ttest.list[[i]],
                                            selection.method = "vst", nfeatures = gene_num, verbose = FALSE)
  }

  anchors <- FindTransferAnchors(reference = ttest.list[[1]],
                                 query = ttest.list[[2]],
                                 dims = 1:30,
                                 k.filter = 100,
                                 features = VariableFeatures(ttest.list[[1]]))

  predictions <- TransferData(anchorset = anchors,
                              refdata = ttest.list[[1]]$celltype,
                              dims = 1:30)

  tibble(
    ori = ori,
    prd = predictions$predicted.id,
    #prob = predictions$prediction.score.max,
    method = 'Seurat'
  )
}

pipe_fun <- function(ID, gene_num){

  train <- expr[ID,]
  test <- expr[-ID,]

  ck1 <- scmap_ck(train, test, gene_num)
  ck2 <- SciBet_ck(train, test, gene_num)
  ck3 <- Seurat3(train, test, gene_num)

  ck1 %>%
    dplyr::bind_rows(ck2) %>%
    dplyr::bind_rows(ck3) -> res

  return(res)
}

res <- list()
ran <- list()
fi.res <- list()

gene.numbers <- c(500)
for (j in 1:length(gene.numbers)) {
  for(i in 1:50){
    tibble(num = 1:nrow(expr),
           label = expr$label) %>%
      dplyr::group_by(label) %>%
      dplyr::sample_frac(0.7) %>%
      dplyr::ungroup() %>%
      dplyr::pull(num) -> ran[[i]]
  }

  for (i in 1:50) {
    res[[i]] <- pipe_fun(ran[[i]], gene.numbers[j])
    print(i)
  }

  tmp <- Reduce(rbind, res)
  tmp <- tmp %>% dplyr::mutate(gene_num = gene.numbers[j])
  fi.res[[j]] <- tmp
}

fi.res %>% readr::write_rds(path = file.path(out.path, file.name), compress = "gz")
