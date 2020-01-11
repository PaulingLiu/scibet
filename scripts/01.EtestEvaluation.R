#To speed up, we run these 14 datasets seperately.

library(SingleCellExperiment)
library(tidyverse)
library(reticulate)
library(tidyverse)
library(ggplot2)
library(scmap)
library(fmsb)
library(scibet)

path_1 <- '/home/pauling/projects/01_classifier/01_data/15_revise/02.CrossValidationAddData/'       #file path
path_2 <- '/home/pauling/projects/01_classifier/01_data/15_revise/04.EtestReRun/'      #result path

GSE <- 'GSE64016.rds.gz'               # file name of each dataset

expr <- readr::read_rds(paste(path_1, GSE, sep = ''))
expr[,-ncol(expr)] <- 1000000*expr[,-ncol(expr)]/rowSums(expr[,-ncol(expr)])
expr <- as.data.frame(expr)


DsGene_fun <- function(expr){
  etestgene <- SelectGene(expr, k=6000, r=T)
  tibble(gene = etestgene)
}

M3Gene_fun <- function(expr){
  ann <- data.frame(cell_type1 = expr$label)
  expr <- expr %>% dplyr::select(-label)
  expr <- t(expr)

  sce <- SingleCellExperiment(assays = list(normcounts = as.matrix(expr)), colData = ann)
  logcounts(sce) <- log2(normcounts(sce) + 1)
  rowData(sce)$feature_symbol <- rownames(sce)
  sce <- sce[!duplicated(rownames(sce)), ]
  sce <- selectFeatures(sce, suppress_plot = T, n_features = 5000)
  tibble(
    gene = rownames(sce),
    score = rowData(sce)[['scmap_scores']]) %>%
    dplyr::arrange(desc(score)) -> M3Gene

  M3Gene
}

Ftest_fun <- function(expr){
  expr[,-ncol(expr)] <- log2(expr[,-ncol(expr)] + 1)
  tibble(gene = colnames(expr)[-ncol(expr)]) %>%
    dplyr::mutate(pval = purrr::map_dbl(
      .x = gene,
      .f = function(.x){
        tibble(
          gene = unlist(expr[, .x]),
          label = expr$label
        ) %>%
          aov(gene~label, data = .) %>%
          summary() -> res

        res[[1]][[5]][1]
      }
    )
    ) -> FtestGene

  FtestGene %>%
    dplyr::arrange(pval) %>%
    tibble::as.tibble()
}

scmap_ck <- function(expr_train, expr_test, gene){

  train_label <- expr_train$label
  test_label <- expr_test$label

  expr_test <- expr_test %>% dplyr::select(-label) %>% t()

  expr_test <- 2^expr_test - 1

  tx_sce <- SingleCellExperiment(assays = list(normcounts = as.matrix(expr_test)))
  logcounts(tx_sce) <- log2(normcounts(tx_sce) + 1)
  rowData(tx_sce)$feature_symbol <- rownames(tx_sce)

  type_expr <- expr_train %>%
    tibble::as_tibble() %>%
    tidyr::nest(-label) %>%
    dplyr::rename(expr = data) %>%
    dplyr::mutate(colmeans = purrr::map(
      .x = expr,
      .f = function(.x){colMedians(as.matrix(.x))}))

  type_expr$colmeans %>%
    as.data.frame() %>%
    tibble::remove_rownames() %>%
    t() %>%
    as.data.frame() %>%
    tibble::remove_rownames() -> type_mean_expr

  rownames(type_mean_expr) <- type_expr$label
  colnames(type_mean_expr) <- colnames(expr_train)[-ncol(expr_train)]
  type_mean_expr <- t(type_mean_expr)

  tibble(num = 1:nrow(type_mean_expr),
         Gene = colnames(expr_train)[-ncol(expr_train)]) -> gene_index

  cal_ck <- function(.x){
    gene_index %>%
      dplyr::filter(Gene %in% gene[1:.x]) %>%
      dplyr::arrange(Gene) %>%
      dplyr::pull(num) -> ID

    type_mean_expr <- type_mean_expr[ID,]

    scmapCluster_results <- scmapCluster(
      projection = tx_sce,
      threshold = 0,
      index_list = list(
        yan = as.data.frame(type_mean_expr)
      )
    )

    CohenKappa <- tibble(
      ori = as.character(test_label),
      prd = unlist(scmapCluster_results$combined_labs))

    CohenKappa %>%
      dplyr::filter(ori == prd) %>%
      nrow(.) -> accuracy

    accuracy/nrow(CohenKappa)
  }

  tibble(
    gene_num = c(30,50,70,100,200,500,700,1000,1500,2000,5000)
  ) %>%
    dplyr::mutate(ck = purrr::map_dbl(gene_num, cal_ck)) %>%
    dplyr::mutate(method = 'scmap')
}
svm_ck <- function(expr_train, expr_test, gene){
  use_python('/home/pauling/anaconda3/bin/python')
  MNB <- import('sklearn.naive_bayes')
  MNB <- MNB$MultinomialNB
  rf <- import('sklearn.ensemble')
  svm <- import('sklearn.svm')
  RF <- rf$RandomForestClassifier
  SVM <- svm$LinearSVC

  X_label <- expr_train$label
  Y_label <- expr_test$label
  clf <- SVM(max_iter = as.integer(9000))
  cal_ck <- function(.x){
    gene_num <- .x
    X_train <- as.matrix(expr_train[,gene[1:gene_num]])
    X_test <- as.matrix(expr_test[,gene[1:gene_num]])
    clf$fit(X_train, X_label)
    prd <- clf$predict(X_test)
    label <- X_label %>% unique() %>% as.character()

    CohenKappa <- tibble(
      ori = as.character(Y_label),
      prd = prd)

    CohenKappa %>%
      dplyr::filter(ori == prd) %>%
      nrow(.) -> accuracy

    accuracy/nrow(CohenKappa)
  }

  tibble(
    gene_num = c(30,50,70,100,200,500,700,1000,1500,2000,5000)
  ) %>%
    dplyr::mutate(ck = purrr::map_dbl(gene_num, cal_ck)) %>%
    dplyr::mutate(method = 'SVM')
}
RF_ck <- function(expr_train, expr_test, gene){
  use_python('/home/pauling/anaconda3/bin/python')
  MNB <- import('sklearn.naive_bayes')
  MNB <- MNB$MultinomialNB
  rf <- import('sklearn.ensemble')
  svm <- import('sklearn.svm')
  RF <- rf$RandomForestClassifier
  SVM <- svm$LinearSVC
  X_label <- expr_train$label
  Y_label <- expr_test$label
  clf <- RF(n_estimators = as.integer(100))
  cal_ck <- function(.x){
    X_train <- as.matrix(expr_train[,gene[1:.x]])
    X_test <- as.matrix(expr_test[,gene[1:.x]])
    clf$fit(X_train, X_label)
    prd <- clf$predict(X_test)
    label <- X_label %>% unique() %>% as.character()
    CohenKappa <- tibble(
      ori = as.character(Y_label),
      prd = prd)

    CohenKappa %>%
      dplyr::filter(ori == prd) %>%
      nrow(.) -> accuracy

    accuracy/nrow(CohenKappa)
  }

  tibble(
    gene_num = c(30,50,70,100,200,500,700,1000,1500,2000,5000)
  ) %>%
    dplyr::mutate(ck = purrr::map_dbl(gene_num, cal_ck)) %>%
    dplyr::mutate(method = 'RF')
}

pipe_fun <- function(.x){

  tibble(num = 1:nrow(.x),
         label = .x$label) %>%
    dplyr::group_by(label) %>%
    dplyr::sample_frac(0.7) %>%
    dplyr::ungroup() %>%
    dplyr::pull(num) -> ID

  train <- .x[ID,]
  test <- .x[-ID,]

  DsGene <- DsGene_fun(train)
  M3Gene <- M3Gene_fun(train)
  FtestGene <- Ftest_fun(train)
  RandomGene <- tibble(gene = sample(colnames(.x[,-ncol(.x)]),length(DsGene$gene)))

  tibble(
    Method = c('DsGene', 'M3Gene', 'FtestGene', 'RandomGene'),
    data = list(DsGene, M3Gene, FtestGene, RandomGene)
  ) -> GeneSet

  train[,-ncol(train)] <- log2(train[,-ncol(train)] + 1)
  test[,-ncol(test)] <- log2(test[,-ncol(test)] + 1)

  get_ck <- function(data){
    gene <- data$gene
    train <- train[,c(gene,'label')]
    test <- test[,c(gene,'label')]

    ck1 <- scmap_ck(train, test, gene)
    ck2 <- svm_ck(train, test, gene)
    ck3 <- RF_ck(train, test, gene)

    ck1 %>%
      dplyr::bind_rows(ck2) %>%
      dplyr::bind_rows(ck3)
  }

  GeneSet %>%
    dplyr::mutate(ck = purrr::map(.x = data, .f = get_ck)) %>%
    dplyr::select(-data)
}

tmp <- tibble(data = list(expr))
res <- list()
for (i in 1:50) {
  res[[i]] <- tmp %>%
    dplyr::mutate(res = purrr::map(data, pipe_fun)) %>%
    dplyr::select(-data)
}

res %>% readr::write_rds(paste(path_2, GSE, sep = ''), compress = 'gz')
