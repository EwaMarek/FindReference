# na wejscie: dane ekspresji, info o unikalnych probkach
# na wyjsciu: ranking genow, info która kontrola przypisana któremu IRowi, macierz fold-changow, lista probek IR, Entrez ID wszystkich genow,
#             tabela z wartosciami wskaznika dla wszystkich eksperymentow

#' @title Create a stability ranking for genes or miRNAs
#'
#' @param all_exp_data List of expression matrices. The first element of the output of \code{rep_elim} function.
#'
#' @param all_uniq_samples List of character matrices with information about unique samples in each experiment. The second element
#' of the output of \code{rep_elim} function.
#'
#' @param miRNA Logical value indicating if stability ranking should be performed for miRNA or genes. Default value is FALSE what means
#' that ranking will be created for genes.
#'
#' @return Function returns a list. It's content depends on miRNA argument.
#'
#' If miRNA == FALSE, the list contains as folows:
#' \itemize{
#' \item $GeneRanking: stability ranking for genes
#' \item $controls: pairs of treated and control samples
#' \item $FC_data: matrices with fold changes
#' \item $samples: names of treated samples
#' \item $EntrezId: all Entrez Ids
#' }
#' else if miRNA == TRUE the list contains as folows:
#' \itemize{
#' \item $miRNA_ranking: stability ranking for miRNAs
#' \item $controls: pairs of treated and control samples
#' \item $FC_data: matrices with fold changes
#' \item $samples: names of treated samples
#' \item $MIMATids: all MIMAT ids
#' }
#'
#' @description
#' \code{create_ranking} function creates stability ranking of genes or miRNAs expression levels.
#'
#' @seealso
#' \code{\link{rep_elim}}, \code{\link{lFC_in_time}}, \code{\link{get_info_about_used_exp}}
#'
#' @details
#' Note that \code{create_ranking} works properly only if arguments are the output of \code{rep_elim} function so you should use it
#' even if you know that there are no replications in your data. If there are no replications in the data \code{rep_elim} function
#' will just prepared objects to use with \code{create_ranking} function.
#'
#' The ranking is created with the following steps:
#' \itemize{
#' \item Z-scores are calculated for each expression value for each gene within an experiment.
#' \item For each treated sample a control is founded and logarithmic fold change calculated.
#' \item For each gene within an experiment gene score is calculated as ratio between logarithmic control expression
#' and logarithmic fold change.
#' \item Stability index for each gene is calculated as a weighed mean of gene scores from all experiments.
#' }
#'
#' Fold change in time for chosen genes or miRNAs could be plotted with \code{lFC_in_time} function.
#' You can also generate summary for used experimental conditions with \code{get_info_about_used_exp} function.
#'
#' It should be made into consideration that reliable ranking can be produced only with reasonable number of microarray data.
#'
#' @examples
#' \dontrun{
#' ##### Create stability ranking for genes
#'
#' # download data from ArrayExpress database
#' to_download = c("E-GEOD-67309", "E-MTAB-966")
#' my_data = downloadAE(to_download, getwd())
#'
#' # load data
#' platforms = c("Affymetrix", "Agilent")
#' loaded_data = load_multi_data(my_data, platforms)
#'
#' # normalize and annotate
#' norm_data = multi_norm_and_annot(loaded_data$raw_expression_data, platforms)
#'
#' # prepare tables for rep_elim function as shown in details
#' path_to_tables = system.file("extdata", "tables_ex3.rds", package = "FindReference")
#' my_tables = readRDS(path_to_tables)
#'
#' # eliminate replications and prepare object for create_ranking function
#' no_rep_data = rep_elim(norm_data, my_tables)
#'
#' # create ranking
#' gene_ranking = create_ranking(no_rep$noRepData, no_rep$uniqSamples, miRNA = FALSE)
#'
#' ##### Create stability ranking for miRNAs
#'
#' # download data from ArrayExpress database
#' datamiRNA = downloadAE("E-MTAB-5197", "/home/emarek/")
#'
#' # prepare table as shown in details load_miRNA help page
#' path_to_table = system.file("inst/extdata", "miRNA_ex1.rds", package = "FindReference")
#' my_table = readRDS(path_to_table)
#'
#' # load data
#' loaded_data = load_miRNA(my_table, datamiRNA[[1]]$path)
#'
#' # normalize and annotate data
#' norm_data = norm_and_annot_miRNA(loaded_data)
#'
#' # eliminate replications and prepare object for create_ranking function
#' no_rep_data = rep_elim(norm_data, my_table)
#'
#' # create ranking
#' miRNA_ranking = create_ranking(no_rep$noRepData, no_rep$uniqSamples, miRNA = TRUE)
#' }
#'
#' @rdname create_ranking
#'
#' @importFrom org.Hs.eg.db org.Hs.egSYMBOL
#' @importFrom stats sd
#'
#' @export


create_ranking = function(all_exp_data, all_uniq_samples, miRNA = FALSE){

  #############################################################################
  ############################ Calculate z-scores #############################
  #############################################################################
  z_scores_data = rep(list(list()), length(all_exp_data))

  for (i in 1:length(all_exp_data)) {
    if(class(all_exp_data[[i]]) != 'list' & class(all_exp_data[[i]]) != 'character'){
      z_scores_data[[i]] = z_scores(all_exp_data[[i]])
    }else if(class(all_exp_data[[i]]) == 'list'){
      z_scores_data[[i]] = lapply(all_exp_data[[i]], z_scores)
    }
  }


  #############################################################################
  ### Calculate fold change and find control sample for each treated sample ###
  #############################################################################
  FC_data = rep(list(list()), length(z_scores_data))
  C_for_IR = rep(list(list()), length(z_scores_data)) #FC_data
  FCnames = vector("character", length(FC_data))

  for (i in 1:length(z_scores_data)) {
    if(class(all_uniq_samples[[i]]) != "character"){
      FC_dataset = log_fold_change(z_scores_data[[i]], all_uniq_samples[[i]], TRUE)
      FC_data[[i]] = FC_dataset[[1]]
      C_for_IR[[i]] = FC_dataset[[2]]

      if(class(all_uniq_samples[[i]]) == "data.frame"){
        FCnames[i] = as.character(all_uniq_samples[[i]]$Experiment[1])
      }else{
        FCnames[i] = as.character(all_uniq_samples[[i]][[1]]$Experiment[1])
      }


    }else{
      FC_data[[i]] = "No data"
      C_for_IR[[i]] = "No data"
      FCnames[i] = "NA"
    }
  }


  names(FC_data) = FCnames

  #############################################################################
  ######################### Take controls expression ##########################
  #############################################################################
  C_expression = rep(list(list()), length(FC_data))
  for(i in 1:length(FC_data)){
    if(class(C_for_IR[[i]]) != "character"){
      C_expression[[i]] = C_log_exp(z_scores_data[[i]], C_for_IR[[i]])
    }else{
      C_expression[[i]] = NULL
    }
  }

  #############################################################################
  ########## Calculate gene score for each experiment independently ###########
  #############################################################################

  exp_score = rep(list(list()), length(FC_data))

  for(i in 1:length(FC_data)){
    exp_score[[i]] = gene_exp_score(FC_data[[i]], C_expression[[i]], C_for_IR[[i]])
  }

  #############################################################################
  ####################### Find all Entrez or MIMAT Ids ########################
  #############################################################################
  Gene_EntrezId = c() # called Gene_EntrezId but could MIMAT ids if performed for miRNA

  for (i in 1:length(FC_data)) {

    if(class(FC_data[[i]]) == 'matrix'){
      Gene_EntrezId = c(Gene_EntrezId,rownames(FC_data[[i]]))
    }else{
      for (j in 1:length(FC_data[[i]])) {
        Gene_EntrezId = c(Gene_EntrezId, rownames(FC_data[[i]][[j]]))
      }
    }
  }
  Gene_EntrezId = unique(Gene_EntrezId)

  #############################################################################
  ############# Find gene symbol for previous founded Entrez Ids ##############
  #############################################################################

  if(miRNA == FALSE){

    EG2SYM = AnnotationDbi::as.list(org.Hs.egSYMBOL)
    Symbols = unlist(EG2SYM[Gene_EntrezId])
    Symbols = Symbols[Gene_EntrezId] # 19 genów nie ma symbolu

  }

  #############################################################################
  ########## Create one table with gene scores from all experiments ###########
  #############################################################################

  all_IRsamples = c()

  for (i in 1:length(C_for_IR)) {

    if(class(all_uniq_samples[[i]]) != 'character'){

      if(class(C_for_IR[[i]]) != 'list'){
        all_IRsamples = c(all_IRsamples, paste(C_for_IR[[i]][1,], all_uniq_samples[[i]]$Experiment[1], sep = " "))
      }else{
        for (z in 1:length(C_for_IR[[i]])) {
          all_IRsamples = c(all_IRsamples, paste(C_for_IR[[i]][[z]][1,], all_uniq_samples[[i]][[z]]$Experiment[1], sep = " "))
        }
      }
    }
  }

  Gene_Score = array(NA, dim = c(length(Gene_EntrezId), length(all_IRsamples)))
  rownames(Gene_Score) = Gene_EntrezId
  colnames(Gene_Score) = all_IRsamples


  for (i in 1:length(exp_score)) {

    if(class(exp_score[[i]]) == 'matrix'){

      expID = all_uniq_samples[[i]]$Experiment[1]
      cols = which(colnames(Gene_Score) %in% paste(colnames(exp_score[[i]]), expID, sep=" "))
      Gene_Score[rownames(exp_score[[i]]), cols] = exp_score[[i]]

    }else if(class(exp_score[[i]]) == 'list'){

      for (j in 1:length(exp_score[[i]])) {

        expID = all_uniq_samples[[i]][[j]]$Experiment[1]
        cols = which(colnames(Gene_Score) %in% paste(colnames(exp_score[[i]][[j]]), expID, sep=" "))
        Gene_Score[rownames(exp_score[[i]][[j]]), cols] = exp_score[[i]][[j]]
      }
    }
  }

  ### średnia wartość współczynnika dla danego genu, czyli ranking bez wag
  # Mean_Gene_Score = rowMeans(Gene_Score, na.rm = TRUE)                               # srednia/mediana
  #
  # Gene_Ranking = data.frame(ID=Gene_EntrezId, symbol=Symbols, score=Mean_Gene_Score)
  # Gene_Ranking = Gene_Ranking[order(Gene_Ranking$score, decreasing = TRUE),]
  # Gene_Ranking$rank = 1:length(Gene_Ranking$score)


  ### ranking z wagami
  # w = apply(Gene_Score, 1, function(x) dim(Gene_Score)[2]-sum(is.na(x)))
  # W_Gene_Score = w*Mean_Gene_Score
  #
  # W_Gene_Ranking = data.frame(ID=Gene_EntrezId, symbol=Symbols, score=W_Gene_Score)
  # W_Gene_Ranking = W_Gene_Ranking[order(W_Gene_Ranking$score, decreasing = TRUE),]
  # W_Gene_Ranking$rank = 1:length(W_Gene_Ranking$score)

  ### ranking z uwzględnieniem tylko statystyki
  # source("remove_quantiles.R")
  # Gene_Score_stat = t(apply(Gene_Score, 1, remove_quantiles, q_values = c(0.05, 0.95)))
  # Gene_Ranking_stat = data.frame(ID=Gene_EntrezId, symbol=Symbols, score=Mean_Gene_Score_stat)
  # Gene_Ranking_stat = Gene_Ranking_stat[order(Gene_Ranking_stat$score, decreasing = TRUE),]
  # Gene_Ranking_stat$rank = 1:length(Gene_Ranking_stat$score)
  #


  #############################################################################
  ############################ Create Gene Ranking ############################
  #############################################################################

  # remove outliers
  Gene_Score_stat = t(apply(Gene_Score, 1, remove_quantiles, q_values = c(0.05, 0.95)))

  # calculate mean score value for each gene
  Mean_Gene_Score_stat = rowMeans(Gene_Score_stat, na.rm = TRUE)                               # srednia/mediana

  # calculate and apply weights
  w_stat = apply(Gene_Score_stat, 1, function(x){dim(Gene_Score_stat)[2]-sum(is.na(x))} )
  #w_stat = apply(Gene_Score_stat, 1, function(x){y = sqrt(dim(Gene_Score_stat)[2]-sum(is.na(x)))/abs(sd(x, na.rm=TRUE)) ; return(y)}) # zle rzeczy
  W_Gene_Score_stat = w_stat*Mean_Gene_Score_stat

  # create ranking
  if(miRNA == FALSE){
    W_Gene_Ranking_stat = data.frame(ID=Gene_EntrezId, symbol=Symbols, score=W_Gene_Score_stat)
    W_Gene_Ranking_stat = W_Gene_Ranking_stat[-which(is.na(W_Gene_Ranking_stat$symbol) == TRUE | W_Gene_Ranking_stat$score == 'NaN'),]
    W_Gene_Ranking_stat = W_Gene_Ranking_stat[order(W_Gene_Ranking_stat$score, decreasing = TRUE),]
    W_Gene_Ranking_stat$rank = 1:length(W_Gene_Ranking_stat$score)

    return(list(GeneRanking = W_Gene_Ranking_stat, controls = C_for_IR, FC_data = FC_data,
                samples = all_IRsamples, EntrezId = Gene_EntrezId))

  }else{
    W_Gene_Ranking_stat = data.frame(ID=Gene_EntrezId, score=W_Gene_Score_stat)
    W_Gene_Ranking_stat = W_Gene_Ranking_stat[-which(W_Gene_Ranking_stat$score == 'NaN'),]
    W_Gene_Ranking_stat = W_Gene_Ranking_stat[order(W_Gene_Ranking_stat$score, decreasing = TRUE),]
    W_Gene_Ranking_stat$rank = 1:length(W_Gene_Ranking_stat$score)

    return(list(miRNA_ranking = W_Gene_Ranking_stat, controls = C_for_IR, FC_data = FC_data,
                samples = all_IRsamples, MIMATids = Gene_EntrezId))
  }




}
