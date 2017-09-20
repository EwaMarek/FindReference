# na wejście obiekt z danymi (AffyBatch, ExonFeatureSet, marrayRaw, EListRaw) oraz info o platformie
# na wyjściu Expression Set -> dane znormalizowane i z odpowiednią adnotacją

#' @title Normalization and annotation of microarray data
#'
#' @param obiekt An object of class marrayRaw, EListRaw, AffyBatch, ExonFeatureSet or list of these objects with raw expression data.
#'
#' @param platform A character indicating experiment's platform, could be 'Affymetrix' or 'Agilent'.
#'
#' @return An expression matrix with normalized and annotated expression values.
#'
#' @description
#' \code{norm_and_annot_by_platform} function normalizes and annotates data from a single microarray experiment.
#'
#' @details
#' Normalization
#'
#' For Affymetrix platform rma background correction and normalization is performed. For Agilent platform background correction
#' is performed by subtracting background values. Used normalization method for this platform is vsn.
#'
#' Annotation
#'
#' Annotation for Affybatch object is made using custom cdf file automatically downloaded from brainarray and for ExonFeatureSet
#' using huex10sttranscriptcluster.db package. Annotation for Agilent platform is made with org.Hs.eg.db package.
#'
#' The first argument could be a list due to make it possible to process data from an experiment where there are few array types within
#' one platform (Affymetrix or Agilent). To process data from many experiments at once you should consider using
#' \code{multi_norm_and_annot} function.
#'
#' @seealso
#' \code{\link{multi_norm_and_annot}}
#'
#' @examples
#' \dontrun{
#' ### download, load, normalize and annotate data from ArrayExpress
#' data = downloadAE("E-GEOD-21066", getwd())
#' loaded_experiment = load_data(data, "Affymetrix")
#' processed_data = norm_and_annot_by_platform(loaded_experiment[[1]], "Affymetrix")
#' }
#'
#' @rdname norm_and_annot_by_platform
#'
#' @importFrom AnnotationDbi revmap
#' @importFrom AnnotationDbi mapIds
#' @importFrom org.Hs.eg.db org.Hs.egACCNUM
#' @importFrom hgug4110b.db hgug4110b.db
#' @importFrom affy cleancdfname
#' @importFrom affy cdfName
#' @importFrom affy rma
#' @importFrom affy expresso
#' @importFrom utils installed.packages
#' @importFrom utils download.file
#' @importFrom utils install.packages
#' @importFrom Biobase exprs
#' @importFrom vsn vsn2
#' @importFrom Matrix.utils aggregate.Matrix
#'
#'
#' @export

norm_and_annot_by_platform = function(obiekt, platform){

  if(grepl('Agilent', platform)){

    platform = "Agilent"

  }else if(grepl('Affymetrix', platform)){

    platform = "Affymetrix"

  }else{

    platform = platform
    wyjscie = paste("Only 'Agilent' and 'Affymetrix' platforms are supported. You have chosen:", as.character(platform), sep=' ')
  }


  ### Agilent

  if(platform == 'Agilent'){

    if(class(obiekt) == 'list'){

      wyjscie = lapply(obiekt, processAgilent)

    }else{

      wyjscie = processAgilent(obiekt)

    }


  ### Affymetrix

  }else if(platform == 'Affymetrix'){

    if(class(obiekt) == 'list'){

      wyjscie = lapply(obiekt, processAffy)

    }else{

      wyjscie = processAffy(obiekt)

    }
  }

  return(wyjscie)
}
