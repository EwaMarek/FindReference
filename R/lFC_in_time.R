# na wejscie: symbole genow do wyrysowania, EntrezID wszystkich genow wystepujacych na mikromacierzach, info o wszytkich probkach IR, FC_data,
#             ranking, skala: linear, log, czy dotyczy miRNA
# na wyjscie: wykres

#' @title Plot fold change in time
#'
#' @param genes_to_valid Character vector with gene symbols or MIMAT ids (for miRNA) which fold change in time should be plotted.
#'
#' @param dane_ranking_z_scores A list - output of \code{create_ranking} function.
#'
#' @param scale A character indicating the scale of plot. Could be 'linear' for fold change or 'log' for logarithmic fold change.
#' Default value is 'linear'.
#'
#' @param miRNA Logical indicating if plot is created for gene or miRNA ranking.
#'
#' @return Function returns a ggplot object which could be plotted with \code{plot} function as in the examples.
#'
#' @description
#' \code{lFC_in_time} function could be used to plot fold change in time for chosen genes or miRNAs by gene symbol or MIMAT ids respectively.
#'
#' @seealso
#' \code{\link{create_ranking}}
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
#' path_to_tables = system.file("inst/extdata", "tables_ex3.rds", package = "FindReference")
#' my_tables = readRDS(path_to_tables)
#'
#' # eliminate replications and prepare object for create_ranking function
#' no_rep_data = rep_elim(norm_data, my_tables)
#'
#' # create ranking
#' gene_ranking = create_ranking(no_rep$noRepData, no_rep$uniqSamples, miRNA = FALSE)
#'
#' # plot fold change in time for some genes
#' genes_to_plot = c('GAPDH', 'ACTB', 'LYPLA2', 'B2M', 'TP53')
#' genes_FC_in_time = lFC_in_time(genes_to_plot, gene_ranking)
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
#'
#' # plot fold change in time for the most stable miRNAs
#' miRNA_to_plot = miRNA_ranking$miRNA_ranking[1:9, 'ID']
#' genes_FC_in_time = lFC_in_time(miRNA_to_plot, miRNA_ranking, scale = 'linear', miRNA = TRUE)
#' }
#'
#' @rdname lFC_in_time
#'
#' @importFrom reshape2 melt
#' @importFrom stats quantile
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 stat_summary
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 facet_wrap
#'
#'
#' @export


lFC_in_time = function(genes_to_valid, dane_ranking_z_scores, scale = 'linear', all_uniq_samples, miRNA=FALSE){


  #############################################################################
  ######################## Take elements from a list ##########################
  #############################################################################

  if(miRNA == FALSE){
    Gene_EntrezId = dane_ranking_z_scores$EntrezId
    W_Gene_Ranking_stat = dane_ranking_z_scores$GeneRanking
  }else{
    W_Gene_Ranking_stat = dane_ranking_z_scores$miRNA_ranking
    Gene_EntrezId = dane_ranking_z_scores$MIMATids
  }

  all_IRsamples = dane_ranking_z_scores$samples
  FC_data = dane_ranking_z_scores$FC_data


  #############################################################################
  ################## Create matrix with log fold change data ##################
  #############################################################################

  lFC_matrix = array(NA, dim = c(length(Gene_EntrezId), length(all_IRsamples)))
  rownames(lFC_matrix) = Gene_EntrezId
  colnames(lFC_matrix) = all_IRsamples

  for (i in 1:length(FC_data)) {

    if(class(FC_data[[i]]) == 'matrix'){

      #CN = strsplit(colnames(lFC_matrix), " ")
      cols = which(colnames(lFC_matrix) %in% paste(colnames(FC_data[[i]]), names(FC_data[i]), sep = " "))
      lFC_matrix[rownames(FC_data[[i]]), cols] = FC_data[[i]][, unlist(lapply(strsplit(colnames(lFC_matrix)[cols], " "), "[[", 1))]

    }else if(class(FC_data[[i]]) == 'list'){

      for (j in 1:length(FC_data[[i]])) {

        #CN = unlist(strsplit(colnames(lFC_matrix), " "))
        cols = which(colnames(lFC_matrix) %in% paste(colnames(FC_data[[i]][[j]]), names(FC_data[i]), sep = " "))
        #cols = which(colnames(lFC_matrix) %in% colnames(FC_data[[i]][[j]]))
        lFC_matrix[rownames(FC_data[[i]][[j]]), cols] = FC_data[[i]][[j]][, unlist(lapply(strsplit(colnames(lFC_matrix)[cols], " "), "[[", 1))]
        #lFC_matrix[rownames(FC_data[[i]][[j]]), cols] = FC_data[[i]][[j]][, colnames(lFC_matrix)[cols]]
      }
    }
  }


  #############################################################################
  ############## For genes change entrez ids into gene symbols ################
  #############################################################################

  if(miRNA==FALSE){
    EG2SYM = AnnotationDbi::as.list(org.Hs.egSYMBOL)
    Symbols = unlist(EG2SYM[Gene_EntrezId])
    Symbols = Symbols[Gene_EntrezId]
    rownames(lFC_matrix) = Symbols
  }


  #############################################################################
  ######## Extract rows chosen by user and add info about experiment ##########
  #############################################################################

  ### change all unique samples structure - add elements names
  tablesNames = vector(mode = "character", length(all_uniq_samples))

  for(i in 1:length(all_uniq_samples)){

    if(class(all_uniq_samples[[i]]) == 'data.frame'){

      tablesNames[i] = as.character(all_uniq_samples[[i]]$Experiment[1])

    }else if(class(all_uniq_samples[[i]]) == 'list'){

      tablesNames[i] = as.character(all_uniq_samples[[i]][[1]]$Experiment[1])

    }else{

      tablesNames[i] = "NA"

    }
  }

  names(all_uniq_samples) = tablesNames

  ### extract information
  lFC_vs_time = melt(lFC_matrix[as.character(genes_to_valid), ])
  cells = vector(mode="character", length = dim(lFC_vs_time)[1])
  dose = vector(mode="character", length = dim(lFC_vs_time)[1])
  Time = vector(mode="character", length = dim(lFC_vs_time)[1])
  Exp = vector(mode="character", length = dim(lFC_vs_time)[1])

  for(i in 1:dim(lFC_vs_time)[1]){

    samp = as.character(lFC_vs_time[i,2])
    uniTab = all_uniq_samples[[strsplit(samp, " ")[[1]][2]]]

    if(class(uniTab) == 'data.frame'){

      currRow = uniTab[which(uniTab$internalId == strsplit(samp, " ")[[1]][1]),]
      cells[i] = as.character(currRow$CellLine)
      dose[i] = as.character(currRow$Dose)
      Time[i] = as.character(currRow$Time)
      Exp[i] = as.character(currRow$Experiment)

    }else if(class(uniTab) == 'list'){

      if(grepl("A", strsplit(samp, " ")[[1]][1])){

        currRow = uniTab[[1]][which(uniTab[[1]]$internalId == strsplit(samp, " ")[[1]][1]),]

      }else{

        currRow = uniTab[[2]][which(uniTab[[2]]$internalId == strsplit(samp, " ")[[1]][1]),]

      }

      cells[i] = as.character(currRow$CellLine)
      dose[i] = as.character(currRow$Dose)
      Time[i] = as.character(currRow$Time)
      Exp[i] = as.character(currRow$Experiment)

    }

  }

  # add info to main table
  lFC_vs_time$cells = cells
  lFC_vs_time$dose = dose
  lFC_vs_time$Time = Time
  lFC_vs_time$exp = Exp



  #lFC_vs_time$gene = Genes[as.character(lFC_vs_time$Var1),]$symbol
  ranks = genes_to_valid

  if(miRNA==FALSE){
    col_with_names = which(colnames(W_Gene_Ranking_stat) == 'symbol')
  }else{
    col_with_names = 1
    }


  for(i in 1:length(genes_to_valid)){
    ranks[i] = W_Gene_Ranking_stat[which(W_Gene_Ranking_stat[,col_with_names] == genes_to_valid[i]), 'rank']
  }
  lFC_vs_time$label = paste(genes_to_valid, 'rank:', ranks ,sep = ' ')
  lFC_vs_time$linear = 2^(lFC_vs_time$value)



  #############################################################################
  ################ Calculate additional parametres for ggplot #################
  #############################################################################

  # ile wierszy na wykresie
  N_rows = round(sqrt(length(genes_to_valid)))

  # funkcje do wyznaczenia gornego i dolnego kwartyla
  median.quartile1 <- function(x){
    out = quantile(x, probs = 0.25)
    return(out)
  }

  median.quartile3 <- function(x){
    out = quantile(x, probs = 0.75)
    return(out)
  }

  # ustalenie kolejnosci wykresow
  lev_tab = data.frame(pre_label=unique(lFC_vs_time$label))
  rozdzielone = strsplit(as.character(lev_tab$pre_label), ': ', fixed=TRUE)
  lev_tab$num = lapply(rozdzielone, function(x){as.numeric(x[[2]])})
  lev_tab = lev_tab[order(unlist(lev_tab$num)),]

  lFC_vs_time$label2 = factor(lFC_vs_time$label, levels = as.character(lev_tab$pre_label))



  #############################################################################
  #################### Create plot in scale chosen by user ####################
  #############################################################################

  if(scale=='linear'){

    lFC_vs_time = lFC_vs_time[which(is.na(lFC_vs_time$linear) == FALSE),]

    wykres = ggplot(lFC_vs_time, aes(x=Time, y=linear, group = 1))+
      stat_summary(geom = 'ribbon', fun.ymin = median.quartile1, fun.ymax = 'median', colour = 'blue', alpha=0.5)+
      stat_summary(geom = 'ribbon', fun.ymin = 'median', fun.ymax = median.quartile3, colour = 'blue', alpha=0.5)+
      stat_summary(fun.y = "median", colour = "red", size = 0.5, geom = "line")+
      theme_bw() +
      ggtitle("Median, 1. i 3. quartile of fold-change value in time") +
      xlab('Time [h]') +
      ylab('Fold change')+
      facet_wrap(  ~ label2, nrow=N_rows)

  }else if(scale=='log'){

    lFC_vs_time = lFC_vs_time[which(is.na(lFC_vs_time$value) == FALSE),]

    wykres = ggplot(lFC_vs_time, aes(x=Time, y=value, group = 1))+
      stat_summary(geom = 'ribbon', fun.ymin = median.quartile1, fun.ymax = 'median', colour = 'blue', alpha=0.5)+
      stat_summary(geom = 'ribbon', fun.ymin = 'median', fun.ymax = median.quartile3, colour = 'blue', alpha=0.5)+
      stat_summary(fun.y = "median", colour = "red", size = 0.5, geom = "line")+
      theme_bw() +
      ggtitle("Median, 1. i 3. quartile of fold-change value in time") +
      xlab('Time [h]') +
      ylab('Fold change')+
      facet_wrap(  ~ label2, nrow=N_rows)
  }

  return(wykres)
}
