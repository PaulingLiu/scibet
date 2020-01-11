#To speed up, we run these datasets seperately.

library(SingleCellExperiment)
library(reticulate)
library(tidyverse)
library(Seurat)
library(scmap)
library(fmsb)
library(scibet)

path_1 <- '/home/pauling/projects/01_classifier/01_data/10_CrossValidation_NegControl/expr/'     #file path
path_2 <- '/home/pauling/projects/01_classifier/01_data/15_revise/05.negative.contral/'      #result path

GSE <- 'GSE75748.rds.gz'   # file name of each dataset

expr <- readr::read_rds(paste(path_1, GSE, sep = ''))
null_expr <- readr::read_rds('/home/pauling/projects/01_classifier/01_data/10_CrossValidation_NegControl/null_expr.rds.gz')

scmap_ck <- function(expr_train, expr_test, gene_num = 500){ 
  
  train_label <- expr_train$label
  test_label <- expr_test$label
  
  expr_train <- expr_train %>% dplyr::select(-label) %>% t()
  expr_test <- expr_test %>% dplyr::select(-label) %>% t()
  
  expr_train <- 2^expr_train - 1
  expr_test <- 2^expr_test - 1
  
  ann <- data.frame(cell_type1 = train_label)
  sce <- SingleCellExperiment(assays = list(normcounts = as.matrix(expr_train)), colData = ann)
  logcounts(sce) <- log2(normcounts(sce) + 1)
  rowData(sce)$feature_symbol <- rownames(sce)
  sce <- sce[!duplicated(rownames(sce)), ]
  
  tx_sce <- SingleCellExperiment(assays = list(normcounts = as.matrix(expr_test)))
  logcounts(tx_sce) <- log2(normcounts(tx_sce) + 1)
  rowData(tx_sce)$feature_symbol <- rownames(tx_sce)
  
  sce <- selectFeatures(sce, n_features = gene_num,  suppress_plot = T)
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
  
  expr_train[,-ncol(expr_train)] <- 2^expr_train[,-ncol(expr_train)] - 1
  expr_test[,-ncol(expr_test)] <- 2^expr_test[,-ncol(expr_test)] - 1
  
  prd <- SciBet(expr_train, expr_test)
  c_score <- conf_score(ref = expr_train, query = expr_test[,-ncol(expr_test)], null_expr = null_expr, gene_num = gene_num)
  
  tibble(
    ori = as.character(expr_test$label),
    prd = prd,
    prob = c_score,
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
  
  X_train <- 2^X_train-1
  X_test <- 2^X_test-1
  
  matr <- cbind(X_train, X_test)
  metadata <- rbind(metadata, metadata1)
  colnames(matr) <- as.character(1:ncol(matr))
  rownames(metadata) <- as.character(1:nrow(metadata))
  
  ttest <- CreateSeuratObject(counts = matr, meta.data = metadata)
  ttest.list <- SplitObject(object = ttest, split.by = "tech")
  
  for (i in 1:length(x = ttest.list)) {
    ttest.list[[i]] <- NormalizeData(object = ttest.list[[i]], verbose = FALSE)
    ttest.list[[i]] <- FindVariableFeatures(object = ttest.list[[i]], 
                                            k.filter = 100,
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
    prob = predictions$prediction.score.max,
    method = 'Seurat'
  )
}
pipe_fun <- function(.x, .y, ID, gene_num){
  
  over_gene <- intersect(colnames(.x), colnames(.y))
  matr1 <- .x[,over_gene]
  matr2 <- .y[,over_gene]
  matr2$label <- 'Neg'
  
  train <- matr1[ID,]
  test <- matr1[-ID,] %>% dplyr::bind_rows(matr2)
  
  
  
  train[,-ncol(train)] <- log2(train[,-ncol(train)] + 1)
  test[,-ncol(test)] <- log2(test[,-ncol(test)] + 1)
  
  ck1 <- scmap_ck(train, test, gene_num)
  ck2 <- SciBet_ck(train, test, gene_num)
  ck3 <- Seurat3(train, test, gene_num)
  
  ck1 %>%
    dplyr::bind_rows(ck2) %>%
    dplyr::bind_rows(ck3) -> res
  
  return(res)
}

expr_1 <- expr %>%
  dplyr::filter(type == 'expr_1') %>%
  dplyr::select(-type)

expr_2 <- expr %>%
  dplyr::filter(type == 'expr_2') %>%
  dplyr::select(-type)


res <- list()
ran <- list()
fi.res <- list()

gene.numbers <- c(500)
for (j in 1:length(gene.numbers)) {
  for(i in 1:50){
    tibble(num = 1:nrow(expr_1),
           label = expr_1$label) %>%
      dplyr::group_by(label) %>%
      dplyr::sample_frac(0.7) %>%
      dplyr::ungroup() %>%
      dplyr::pull(num) -> ran[[i]]
  }
  
  for (i in 1:50) {
    res[[i]] <- pipe_fun(expr_1, expr_2, ran[[i]], gene.numbers[j])
    print(i)
  }
  
  tmp <- Reduce(rbind, res)
  tmp <- tmp %>% dplyr::mutate(gene_num = gene.numbers[j])
  fi.res[[j]] <- tmp
}

fi.res %>% readr::write_rds(paste(path_2, GSE, sep = ''), compress = 'gz')
