ModuleScore <- function (project, 
                         features, 
                         assay = 'GeneScoreMatrix', 
                         reducedDims, 
                         nbin = 24, 
                         ctrl = 100, 
                         name = "Cluster", 
                         seed = 1){
  if (!is.null(x = seed)) {
    set.seed(seed = seed)
  }
  if (length(getImputeWeights(project)) == 0){ # if there are no imputed weights, add them 
    project <- addImputeWeights(project, reducedDims = reducedDims) 
  }
  assay.atac <- getMatrixFromProject(project, useMatrix = 'GeneScoreMatrix')
  matGS <- imputeMatrix(assay(assay.atac), getImputeWeights(project)) %>% as(., "dgCMatrix") #get imputed gene score
  # assay.data <- assay(assay.atac) # get score data for any matrix that is not ImputedGeneMatrix
  assay.data <- matGS # so I don't need to go back and change everything in the code for dev phase
    
  all.features.input <- rowData(assay.atac)$name
  assay.data@Dimnames[[1]] <- all.features.input
  
  # features.old <- features # not sure why we have this
  if (is.null(x = features)) {
    stop("Missing input feature list")
  }
  features <- lapply(X = features, FUN = function(x) {
    missing.features <- setdiff(x = x, y = all.features.input)
    if (length(x = missing.features) > 0) {
      warning("The following features are not present in the object: ", 
              paste(missing.features, collapse = ", "), ", not searching for symbol synonyms", 
      call. = FALSE, immediate. = TRUE)
    }
    return(intersect(x = x, y = all.features.input))
  }) # checking for missing genes
  cluster.length <- length(x = features)

  if (!all(LengthCheck(values = features))) {
    warning(paste("Could not find enough features in the object from the following feature lists:", 
                  paste(names(x = which(x = !LengthCheck(values = features)))), 
                  "Attempting to match case..."))
    features <- lapply(X = features.old, FUN = CaseMatch, 
                       match = all.features.input)
  } # making sure there are enough genes per list 
  if (!all(LengthCheck(values = features))) {
    stop(paste("The following feature lists do not have enough features present in the object:", 
               paste(names(x = which(x = !LengthCheck(values = features)))), 
               "exiting..."))
  } # making sure there are enough genes per list 
  
  # pool <- pool %||% all.features.input
  # data.avg <- Matrix::rowMeans(x = assay.data[pool, ])
  data.avg <- Matrix::rowMeans(x = assay.data) # we are going to hardcode this to always use all genes as pool
  data.avg <- data.avg[order(data.avg)] #ordering the genes based on average expression for binning
  
  data.cut <- cut_number(x = data.avg + rnorm(n = length(data.avg))/1e+30, 
                         n = nbin, labels = FALSE, right = FALSE) # get binned data, adding normally distributed error incorporating number of genes
  names(x = data.cut) <- names(x = data.avg)
  
  # getting (random) control features from bins
  ctrl.use <- vector(mode = "list", length = cluster.length)
  # for each gene (in each list), get n (ctrl) control genes from the same bin. Same bin ~ similar expression levels
  for (i in 1:cluster.length) { # for each list entry ~ group of genes
    features.use <- features[[i]]
    for (j in 1:length(x = features.use)) { # for each gene
      ctrl.use[[i]] <- c(ctrl.use[[i]], names(x = sample(x = data.cut[which(x = data.cut == 
                                                                              data.cut[features.use[j]])], size = ctrl, replace = FALSE)))
      
      }
  }
  
  ctrl.use <- lapply(X = ctrl.use, FUN = unique) #make list unique 
  ctrl.scores <- matrix(data = numeric(length = 1L), nrow = length(x = ctrl.use), 
                        ncol = ncol(x = assay.data))
  for (i in 1:length(ctrl.use)) { # score for all control features (average expression for all genes per cell)
    features.use <- ctrl.use[[i]]
    ctrl.scores[i, ] <- Matrix::colMeans(x = assay.data[features.use,])
  }
  features.scores <- matrix(data = numeric(length = 1L), nrow = cluster.length, 
                            ncol = ncol(x = assay.data))
  for (i in 1:cluster.length) {  # score for features in list (average expression for all genes per cell)
    features.use <- features[[i]]
    data.use <- assay.data[features.use, , drop = FALSE]
    features.scores[i, ] <- Matrix::colMeans(x = data.use)
  }
  features.scores.use <- features.scores - ctrl.scores # subtract control scores from feature scores (~background correction)
  # rownames(x = features.scores.use) <- paste0(name, 1:cluster.length)
  # features.scores.use <- as.data.frame(x = t(x = features.scores.use))
  # rownames(x = features.scores.use) <- colnames(x = assay.data)
  project <- addCellColData(ArchRProj = project, 
                              data = t(features.scores.use), 
                              name = name, 
                              cells = colnames(x = assay.data),
                              force = T)
  # object[[colnames(x = features.scores.use)]] <- features.scores.use
  return(project)
}

################################################################################
# Check the length of components of a list
#
# @param values A list whose components should be checked
# @param cutoff A minimum value to check for
#
# @return a vector of logicals
#
LengthCheck <- function(values, cutoff = 0) {
  return(vapply(
    X = values,
    FUN = function(x) {
      return(length(x = x) > cutoff)
    },
    FUN.VALUE = logical(1)
  ))
}

################################################################################
#' Match the case of character vectors
#'
#' @param search A vector of search terms
#' @param match A vector of characters whose case should be matched
#'
#' @return Values from search present in match with the case of match
#'
#' @export
#' @concept utilities
#'
#' @examples
#' data("pbmc_small")
#' cd_genes <- c('Cd79b', 'Cd19', 'Cd200')
#' CaseMatch(search = cd_genes, match = rownames(x = pbmc_small))
#'
CaseMatch <- function(search, match) {
  search.match <- sapply(
    X = search,
    FUN = function(s) {
      return(grep(
        pattern = paste0('^', s, '$'),
        x = match,
        ignore.case = TRUE,
        perl = TRUE,
        value = TRUE
      ))
    }
  )
  return(unlist(x = search.match))
}

